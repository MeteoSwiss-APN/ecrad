! radiation_mcica_lw.F90 - Monte-Carlo Independent Column Approximation longtwave solver
!
! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-07-12  R. Hogan  Call fast adding method if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_mcica_lw

  public

contains

  !---------------------------------------------------------------------
  ! Longwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_lw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb,jprd
    use yomhook,  only           : lhook, dr_hook

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_two_stream_gammas_lw, &
         &                               calc_two_stream_gammas_lw_lr, &
         &                               calc_reflectance_transmittance_lw, &
         &                               calc_reflectance_transmittance_lw_lr, &
         &                               calc_no_scattering_transmittance_lw, &
         &                               calc_no_scattering_transmittance_lw_lr, &
#ifdef FAST_EXPONENTIAL
         &                               exp_fast, &
#endif
         &                               LwDiffusivity

#ifndef FAST_EXPONENTIAL
#define exp_fast exp
#endif

    use radiation_adding_ica_lw, only  : adding_ica_lw, adding_ica_lw_lr, adding_ica_lw_cond_lr, fast_adding_ica_lw, &
         &                               fast_adding_ica_lw_lr, calc_fluxes_no_scattering_lw, &
         &                               calc_fluxes_no_scattering_lw_lr, calc_fluxes_no_scattering_lw_cond_lr
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator, only: cloud_generator_lr

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    ! cos: original (g, lev). Future demote to (col, lev)
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: ref_clear, reflectance, &
                                              trans_clear, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    ! cos: original (g, lev). Future demote to (col, lev)
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: source_up_clear, source_up, &
                                              source_dn_clear, source_dn

    ! Fluxes per g point
    ! cos: temporarily promote to jcol
    real(jprb), dimension(config%n_g_lw, nlev+1,istartcol:iendcol) :: flux_up, flux_dn
    ! cos: original (g, lev). Future demote to (col, lev)
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up_clear, flux_dn_clear

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    ! cos: original (ng). Future demote to jcol
    ! cos: temporarily add nlev for loop splitting
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_total, ssa_total, g_total

    ! Combined scattering optical depth
    real(jprb) :: scat_od, scat_od_total(config%n_g_lw)

    ! Two-stream coefficients
    ! cos: original (g). Future demote to (col)
    real(jprb), dimension(istartcol:iendcol) :: gamma1, gamma2 

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    ! cos; original (g, nlev), Future demote to (jcol, nlev)
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    ! cos: temporarily added jcol & nlev for loop splitting
    real(jprb), dimension(config%n_g_lw, nlev,istartcol:iendcol) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    ! cos: original (scalar). Future demote to (scalar) again
    real(jprb), dimension(istartcol:iendcol) :: total_cloud_cover

    ! Identify clear-sky layers
    !cos : temporarily added jcol
    logical :: is_clear_sky_layer(nlev,istartcol:iendcol)

    ! Index of the highest cloudy layer
    ! cos : temporarily added jcol
    integer :: i_cloud_top(istartcol:iendcol)

    ! Number of g points
    integer :: ng

    ! Loop indices for level, column and g point
    integer :: jlev, jcol, jg

    ! cos: inlining of functions
    real(jprb) :: factor
    real(jprd) :: k_exponent, reftrans_factor
    real(jprd) :: exponential  ! = exp(-k_exponent*od)
    real(jprd) :: exponential2 ! = exp(-2*k_exponent*od)

    real(jprd) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot
    ! end cos

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: longwave McICA requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

    ng = config%n_g_lw

    ! Loop through columns
    do jg = 1, ng

      ! Clear-sky calculation
      if (config%do_lw_aerosol_scattering) then
        ! Scattering case: first compute clear-sky reflectance,
        ! transmittance etc at each model level
        do jlev = 1,nlev
          ssa_total = ssa(:,:,:)
          g_total   = g(:,:,:)
          call calc_two_stream_gammas_lw_lr(istartcol, iendcol, ssa_total(jg,jlev,:), g_total(jg,jlev,:), &
               &  gamma1, gamma2)
          call calc_reflectance_transmittance_lw_lr(istartcol, iendcol, &
               &  od(jg,jlev,:), gamma1, gamma2, &
               &  planck_hl(jg,jlev,:), planck_hl(jg,jlev+1,:), &
               &  ref_clear(jg,jlev,:), trans_clear(jg,jlev,:), &
               &  source_up_clear(jg,jlev,:), source_dn_clear(jg,jlev,:))
        end do
        ! Then use adding method to compute fluxes
        call adding_ica_lw_lr(istartcol, iendcol, nlev, &
             &  ref_clear(jg,:,:), trans_clear(jg,:,:), source_up_clear(jg,:,:), source_dn_clear(jg,:,:), &
             &  emission(jg,:), albedo(jg,:), &
             &  flux_up_clear(jg,:,:), flux_dn_clear(jg,:,:))
        
      else
        ! ! Non-scattering case: use simpler functions for
        ! ! transmission and emission
        do jlev = 1,nlev
          call calc_no_scattering_transmittance_lw_lr(istartcol, iendcol, od(jg,jlev,:), &
               &  planck_hl(jg,jlev,:), planck_hl(jg,jlev+1,:), &
               &  trans_clear(jg,jlev,:), source_up_clear(jg,jlev,:), source_dn_clear(jg,jlev,:))
        end do
        ! ! ! Simpler down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw_lr(istartcol, iendcol, nlev, &
             &  trans_clear(jg,:,:), source_up_clear(jg,:,:), source_dn_clear(jg,:,:), &
             &  emission(jg,:), albedo(jg,:), &
             &  flux_up_clear(jg,:,:), flux_dn_clear(jg,:,:))
        
        ! Ensure that clear-sky reflectance is zero since it may be
        ! used in cloudy-sky case
        ref_clear = 0.0_jprb
      end if
    ! cos: todo once all fields are promoted to 3D
    enddo
    ! cos: array syntax once data layouts are compatible
    do jcol = istartcol,iendcol

      ! Sum over g-points to compute broadband fluxes
      flux%lw_up_clear(jcol,:) = sum(flux_up_clear(:,:,jcol),1)
      flux%lw_dn_clear(jcol,:) = sum(flux_dn_clear(:,:,jcol),1)
      flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(:,nlev+1,jcol)
    enddo

      ! Do cloudy-sky calculation; add a prime number to the seed in
      ! the longwave
    call cloud_generator_lr(ng, istartcol, iendcol, nlev, config%i_overlap_scheme, &
           &  single_level%iseed, &
           &  config%cloud_fraction_threshold, &
           &  cloud%fraction, cloud%overlap_param, &
           &  config%cloud_inhom_decorr_scaling, cloud%fractional_std, &
           &  config%pdf_sampler, od_scaling, total_cloud_cover, &
           &  is_beta_overlap=config%use_beta_overlap)
      
    do jcol = istartcol, iendcol
      ! Store total cloud cover
      flux%cloud_cover_lw(jcol) = total_cloud_cover(jcol)      
    enddo
    do jcol = istartcol,iendcol
      if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
        ! Total-sky calculation

        is_clear_sky_layer(:,jcol) = .true.
        i_cloud_top(jcol) = nlev+1
        do jlev = 1,nlev
          ! Compute combined gas+aerosol+cloud optical properties
          if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
            is_clear_sky_layer(jlev,jcol) = .false.
            ! Get index to the first cloudy layer from the top
            if (i_cloud_top(jcol) > jlev) then
              i_cloud_top(jcol) = jlev
            end if
          endif
        enddo
      endif
    enddo

    do jg=1,ng
      do jlev = 1,nlev
        do jcol = istartcol,iendcol
        ! Compute combined gas+aerosol+cloud optical properties
          if ((total_cloud_cover(jcol) >= config%cloud_fraction_threshold) .and. &
&               (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold)) then
            od_cloud_new(jg,jlev,jcol) = od_scaling(jg,jlev, jcol) &
                &  * od_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol)
            od_total(jg,jlev,jcol) = od(jg,jlev,jcol) + od_cloud_new(jg,jlev,jcol)
            ssa_total(jg,jlev,jcol) = 0.0_jprb
            g_total(jg,jlev,jcol)   = 0.0_jprb
          endif
        enddo
      enddo

      do jlev = 1,nlev
        if (config%do_lw_cloud_scattering) then
          ! Scattering case: calculate reflectance and
          ! transmittance at each model level
          do jcol = istartcol,iendcol
            if ((total_cloud_cover(jcol) >= config%cloud_fraction_threshold) .and. &
&               (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold)) then

              if (config%do_lw_aerosol_scattering) then
                ! In single precision we need to protect against the
                ! case that od_total > 0.0 and ssa_total > 0.0 but
                ! od_total*ssa_total == 0 due to underflow
                scat_od_total(jg) = ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                    &     + ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                    &     *  od_cloud_new(jg,jlev,jcol)
                if (scat_od_total(jg) > 0.0_jprb) then
                   g_total(jg,jlev,jcol) = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                       &     +   g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                       &     *  od_cloud_new(jg,jlev,jcol)) &
                       &     / scat_od_total(jg)
                endif                
                if (od_total(jg,jlev,jcol) > 0.0_jprb) then
                   ssa_total(jg,jlev,jcol) = scat_od_total(jg) / od_total(jg,jlev,jcol)
                endif
              else
  !                  do jg = 1,ng
                  if (od_total(jg,jlev,jcol) > 0.0_jprb) then
                    scat_od = ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                        &     * od_cloud_new(jg,jlev,jcol)
                    ssa_total(jg,jlev,jcol) = scat_od / od_total(jg,jlev,jcol)
                    if (scat_od > 0.0_jprb) then
                      g_total(jg,jlev,jcol) = g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                          &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                          &     *  od_cloud_new(jg,jlev,jcol) / scat_od
                    end if
                  end if
                !end do
              end if
            endif
          enddo
        endif
      enddo

      do jlev = 1,nlev
        if (config%do_lw_cloud_scattering) then
      
          ! Compute cloudy-sky reflectance, transmittance etc at
          ! each model level
! cos: inlining the function due to the conditionals on the cloud cover
!              call calc_two_stream_gammas_lw(ng, ssa_total(:,jcol), g_total(:,jcol), &
!                   &  gamma1(:,jcol), gamma2(:,jcol))
          do jcol = istartcol,iendcol
            if ((total_cloud_cover(jcol) >= config%cloud_fraction_threshold) .and. &
&               (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold)) then

              ! Fu et al. (1997), Eq 2.9 and 2.10:
              !      gamma1(jg) = LwDiffusivity * (1.0_jprb - 0.5_jprb*ssa(jg) &
              !           &                    * (1.0_jprb + g(jg)))
              !      gamma2(jg) = LwDiffusivity * 0.5_jprb * ssa(jg) &
              !           &                    * (1.0_jprb - g(jg))
              ! Reduce number of multiplications
              factor = (LwDiffusivity * 0.5_jprb) * ssa_total(jg,jlev,jcol)
              gamma1(jcol) = LwDiffusivity - factor*(1.0_jprb + g_total(jg,jlev,jcol))
              gamma2(jcol) = factor * (1.0_jprb - g_total(jg,jlev,jcol))
            endif

! cos: inlining the function due to the conditionals on the cloud cover
!              call calc_reflectance_transmittance_lw(ng, &
!                   &  od_total(:,jcol), gamma1(:,jcol), gamma2(:,jcol), &
!                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
!                   &  reflectance(:,jlev,jcol), transmittance(:,jlev,jcol), &
!                   & source_up(:,jlev,jcol), source_dn(:,jlev,jcol))

            if ((total_cloud_cover(jcol) >= config%cloud_fraction_threshold) .and. &
&               (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold)) then

              if (od_total(jg,jlev,jcol) > 1.0e-3_jprd) then
                k_exponent = sqrt(max((gamma1(jcol) - gamma2(jcol)) * (gamma1(jcol) + gamma2(jcol)), &
                      1.E-12_jprd)) ! Eq 18 of Meador & Weaver (1980)
                exponential = exp_fast(-k_exponent*od_total(jg,jlev,jcol))
                exponential2 = exponential*exponential
                reftrans_factor = 1.0 / (k_exponent + gamma1(jcol) + (k_exponent - gamma1(jcol))*exponential2)
                ! Meador & Weaver (1980) Eq. 25
                reflectance(jg,jlev,jcol) = gamma2(jcol) * (1.0_jprd - exponential2) * reftrans_factor
                ! Meador & Weaver (1980) Eq. 26
                transmittance(jg,jlev,jcol) = 2.0_jprd * k_exponent * exponential * reftrans_factor
              
                ! Compute upward and downward emission assuming the Planck
                ! function to vary linearly with optical depth within the layer
                ! (e.g. Wiscombe , JQSRT 1976).
        
                ! Stackhouse and Stephens (JAS 1991) Eqs 5 & 12
                coeff = (planck_hl(jg,jlev+1,jcol)-planck_hl(jg,jlev,jcol)) / & 
                &       (od_total(jg,jlev,jcol)*(gamma1(jcol)+gamma2(jcol)))
                coeff_up_top  =  coeff + planck_hl(jg,jlev,jcol)
                coeff_up_bot  =  coeff + planck_hl(jg,jlev+1,jcol)
                coeff_dn_top  = -coeff + planck_hl(jg,jlev,jcol)
                coeff_dn_bot  = -coeff + planck_hl(jg,jlev+1,jcol)
                source_up(jg,jlev,jcol) =  coeff_up_top - reflectance(jg,jlev,jcol) * coeff_dn_top - & 
                &                     transmittance(jg,jlev,jcol) * coeff_up_bot
                source_dn(jg,jlev,jcol) =  coeff_dn_bot - reflectance(jg,jlev,jcol) * coeff_up_bot - &
                &                     transmittance(jg,jlev,jcol) * coeff_dn_top
              else
                k_exponent = sqrt(max((gamma1(jcol) - gamma2(jcol)) * (gamma1(jcol) + gamma2(jcol)), &
                      1.E-12_jprd)) ! Eq 18 of Meador & Weaver (1980)
                reflectance(jg,jlev,jcol) = gamma2(jcol) * od_total(jg,jlev,jcol)
                transmittance(jg,jlev,jcol) = (1.0_jprb - k_exponent*od_total(jg,jlev,jcol)) / (1.0_jprb + &
                &                             od_total(jg,jlev,jcol)*(gamma1(jcol)-k_exponent))
                source_up(jg,jlev,jcol) = (1.0_jprb - reflectance(jg,jlev,jcol) - transmittance(jg,jlev,jcol)) &
                      &       * 0.5 * (planck_hl(jg,jlev,jcol) + planck_hl(jg,jlev+1,jcol))
                source_dn(jg,jlev,jcol) = source_up(jg,jlev,jcol)
              end if
            endif
          end do
        endif
      enddo

      do jlev = 1,nlev
        if (config%do_lw_cloud_scattering) then
      
        else
! cos: inlining the function due to the conditionals on the cloud cover

          ! No-scattering case: use simpler functions for
          ! transmission and emission
!              call calc_no_scattering_transmittance_lw(ng, od_total(:,jcol), &
!                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
!                   &  transmittance(:,jlev,jcol), source_up(:,jlev,jcol), source_dn(:,jlev,jcol))

          do jcol = istartcol,iendcol
            if ((total_cloud_cover(jcol) >= config%cloud_fraction_threshold) .and. & 
&             (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold)) then

              ! Compute upward and downward emission assuming the Planck
              ! function to vary linearly with optical depth within the layer
              ! (e.g. Wiscombe , JQSRT 1976).
              if (od_total(jg,jlev,jcol) > 1.0e-3) then
                ! Simplified from calc_reflectance_transmittance_lw above
                coeff = LwDiffusivity*od_total(jg,jlev,jcol)
                transmittance(jg,jlev,jcol) = exp_fast(-coeff)
                coeff = (planck_hl(jg,jlev+1,jcol)-planck_hl(jg,jlev,jcol)) / coeff
                coeff_up_top  =  coeff + planck_hl(jg,jlev,jcol)
                coeff_up_bot  =  coeff + planck_hl(jg,jlev+1,jcol)
                coeff_dn_top  = -coeff + planck_hl(jg,jlev,jcol)
                coeff_dn_bot  = -coeff + planck_hl(jg,jlev+1,jcol)
                source_up(jg,jlev,jcol) =  coeff_up_top - transmittance(jg,jlev,jcol) * coeff_up_bot
                source_dn(jg,jlev,jcol) =  coeff_dn_bot - transmittance(jg,jlev,jcol) * coeff_dn_top
              else
                ! Linear limit at low optical depth
                coeff = LwDiffusivity*od_total(jg,jlev,jcol)
                transmittance(jg,jlev,jcol) = 1.0_jprb - coeff
                source_up(jg,jlev,jcol) = coeff * 0.5_jprb * (planck_hl(jg,jlev,jcol)+planck_hl(jg,jlev+1,jcol))
                source_dn(jg,jlev,jcol) = source_up(jg,jlev,jcol)
              end if
            endif
          end do      
        end if

        do jcol = istartcol,iendcol

          if ((total_cloud_cover(jcol) >= config%cloud_fraction_threshold) .and. &
&           (cloud%fraction(jcol,jlev) < config%cloud_fraction_threshold)) then

            ! Clear-sky layer: copy over clear-sky values
            reflectance(jg,jlev,jcol) = ref_clear(jg,jlev, jcol)
            transmittance(jg,jlev,jcol) = trans_clear(jg,jlev,jcol)
            source_up(jg,jlev,jcol) = source_up_clear(jg,jlev,jcol)
            source_dn(jg,jlev,jcol) = source_dn_clear(jg,jlev,jcol)
          endif
        enddo
      end do
        if (config%do_lw_aerosol_scattering) then
          ! Use adding method to compute fluxes for an overcast sky,
          ! allowing for scattering in all layers
          call adding_ica_lw_cond_lr(istartcol, iendcol, nlev, total_cloud_cover, config%cloud_fraction_threshold, &
&          reflectance(jg,:,:), transmittance(jg,:,:), source_up(jg,:,:), &
&          source_dn(jg,:,:), emission(jg,:), albedo(jg,:), flux_up(jg,:,:), flux_dn(jg,:,:))
        endif
!cos end split
        if (config%do_lw_aerosol_scattering) then
        else if (config%do_lw_cloud_scattering) then
          ! Use adding method to compute fluxes but optimize for the
          ! presence of clear-sky layers
!                      call fast_adding_ica_lw(ng, nlev, &
! &               reflectance(:,:,jcol), transmittance(:,:,jcol), source_up(:,:,jcol), &
!                 & source_dn(:,:,jcol), emission(:,jcol), albedo(:,jcol), is_clear_sky_layer(:,jcol), i_cloud_top(jcol), &
!                 & flux_dn_clear(:,:,jcol), flux_up(:,:,jcol), flux_dn(:,:,jcol))

           call fast_adding_ica_lw_lr(istartcol,iendcol, nlev, total_cloud_cover, config%cloud_fraction_threshold, &
&               reflectance(jg,:,:), transmittance(jg,:,:), source_up(jg,:,:), &
                & source_dn(jg,:,:), emission(jg,:), albedo(jg,:), is_clear_sky_layer(:,:), i_cloud_top, &
                & flux_dn_clear(jg,:,:), flux_up(jg,:,:), flux_dn(jg,:,:))
        endif
        if (config%do_lw_aerosol_scattering) then
        else if (config%do_lw_cloud_scattering) then

        else
          ! ! Simpler down-then-up method to compute fluxes
          ! call calc_fluxes_no_scattering_lw(ng, nlev, &
          !      &  transmittance(:,:,jcol), source_up(:,:,jcol), source_dn(:,:,jcol), emission(:,jcol), albedo(:,jcol), &
          !      &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
                         ! Simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw_cond_lr(istartcol,iendcol, nlev, total_cloud_cover, config%cloud_fraction_threshold, &
          &  transmittance(jg,:,:), source_up(jg,:,:), source_dn(jg,:,:), emission(jg,:), albedo(jg,:), &
          &  flux_up(jg,:,:), flux_dn(jg,:,:))

        end if
    enddo

    ! cos: here there are reductions on ng. Therefore we need to break the ng loop,
    ! and the storages can only be ng independent if the accumulation is performed in previous loops

    do jcol = istartcol,iendcol
      if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
        
        ! Store overcast broadband fluxes
        flux%lw_up(jcol,:) = sum(flux_up(:,:,jcol),1)
        flux%lw_dn(jcol,:) = sum(flux_dn(:,:,jcol),1)

        ! Cloudy flux profiles currently assume completely overcast
        ! skies; perform weighted average with clear-sky profile
        flux%lw_up(jcol,:) =  total_cloud_cover(jcol) *flux%lw_up(jcol,:) &
             &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) =  total_cloud_cover(jcol) *flux%lw_dn(jcol,:) &
             &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_dn_clear(jcol,:)
        ! Store surface spectral downwelling fluxes
        flux%lw_dn_surf_g(:,jcol) = total_cloud_cover(jcol)*flux_dn(:,nlev+1,jcol) &
             &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_dn_surf_clear_g(:,jcol)

        ! Compute the longwave derivatives needed by Hogan and Bozzo
        ! (2015) approximate radiation update scheme
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance(:,:,jcol), flux_up(:,nlev+1,jcol), &
               &                       flux%lw_derivatives)

          if (total_cloud_cover(jcol) < 1.0_jprb - config%cloud_fraction_threshold) then
            ! Modify the existing derivative with the contribution from the clear sky
            call modify_lw_derivatives_ica(ng, nlev, jcol, trans_clear(:,:,jcol), flux_up_clear(:,nlev+1,jcol), &
                 &                         1.0_jprb-total_cloud_cover(jcol), flux%lw_derivatives)
          end if
        end if
      else
        ! No cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
        flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
        flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, trans_clear(:,:,jcol), flux_up_clear(:,nlev+1,jcol), &
               &                       flux%lw_derivatives)
 
        end if
      end if ! Cloud is present in profile
    end do

    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',1,hook_handle)
    
  end subroutine solver_mcica_lw

end module radiation_mcica_lw
