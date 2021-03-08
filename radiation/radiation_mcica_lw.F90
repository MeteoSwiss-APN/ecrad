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

  implicit none
  private

  public :: solver_mcica_lw

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

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook
    USE YOERRTM  , ONLY : JPGPT

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_two_stream_gammas_lw, &
         &                               calc_reflectance_transmittance_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, fast_adding_ica_lw, &
         &                               calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator, only: cloud_generator

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

    ! Number of g points
#ifndef _OPENACC
    integer :: ng
#else
    integer, parameter :: ng = JPGPT;
    integer, parameter :: WARP_SIZE = 32
    integer, parameter :: NUM_WARPS = merge(ng / WARP_SIZE, ng / WARP_SIZE + 1, mod(ng, WARP_SIZE) == 0)
#endif

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_lw,nlev, istartcol:iendcol) :: od_scaling

    ! Total cloud cover output from the cloud generator
    real(jprb), dimension(istartcol:iendcol) :: total_cloud_cover

    ! Loop indices for level, column and g point
    integer :: jcol, ncol

    real(jprb) :: hook_handle
#ifdef DR_HOOK
    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',0,hook_handle)
#endif
    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: longwave McICA requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

#ifndef _OPENACC
    ng = config%n_g_lw
#else
    if (ng /= config%n_g_lw) then
        stop 'Mismatch between config and assumed n_g_lw'
    end if
#endif

    do jcol = istartcol,iendcol
      ! Do cloudy-sky calculation; add a prime number to the seed in
      ! the longwave
      call cloud_generator(ng, nlev, config%i_overlap_scheme, &
           &  single_level%iseed(jcol) + 997, &
           &  config%cloud_fraction_threshold, &
           &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
           &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
           &  config%pdf_sampler, od_scaling(:, :, jcol), total_cloud_cover(jcol), &
           &  is_beta_overlap=config%use_beta_overlap)
    end do
    ! Store total cloud cover
    flux%cloud_cover_lw(istartcol:iendcol) = total_cloud_cover
#ifdef _OPENACC
    call init_gpu()
#endif
    !$acc data copyin(od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, emission, albedo,&
    !$acc& flux, flux%lw_up_clear, flux%lw_dn_clear, flux%lw_dn_surf_clear_g, flux%lw_up, flux%lw_dn, flux%lw_dn_surf_g,&
    !$acc& config, config%i_band_from_reordered_g_lw, cloud, cloud%fraction, od_cloud, ssa_cloud, g_cloud,  od_scaling, total_cloud_cover)
    ncol = iendcol - istartcol + 1

    !$acc parallel default(none) num_gangs(ncol) num_workers(NUM_WARPS) vector_length(WARP_SIZE)
    ! Loop through columns
    !$acc loop gang
    do jcol = istartcol,iendcol
        call solver_mcica_lw_col(jcol, nlev, &
           &  config, cloud, &
           &  od(:, :, jcol), ssa(:, :, jcol), g(:, :, jcol), od_cloud(:, :, jcol), ssa_cloud(:, :, jcol), &
           & g_cloud(:, :, jcol), planck_hl(:, :, jcol), emission(:, jcol), albedo(:, jcol), &
           &  flux, od_scaling(:, :, jcol), total_cloud_cover(jcol))
    end do
    !$acc end parallel


    !$acc update host(flux%lw_up_clear(istartcol:iendcol,:), flux%lw_dn_clear(istartcol:iendcol,:), flux%lw_dn_surf_clear_g(:,istartcol:iendcol))
    !$acc update host(flux%lw_up(istartcol:iendcol,:), flux%lw_dn(istartcol:iendcol,:), flux%lw_dn_surf_g(istartcol:iendcol,:))
    !$acc end data
#ifdef DR_HOOK
    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',1,hook_handle)
#endif
  end subroutine solver_mcica_lw

  pure subroutine solver_mcica_lw_col(jcol, nlev, &
       &  config, cloud, &
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux, od_scaling, total_cloud_cover)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook
    USE YOERRTM  , ONLY : JPGPT

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_two_stream_gammas_lw, &
         &                               calc_reflectance_transmittance_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, fast_adding_ica_lw, &
         &                               calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator, only: cloud_generator

    implicit none

    ! Inputs
    integer, intent(in) :: jcol
    integer, intent(in) :: nlev               ! number of model levels
    type(config_type),        intent(in) :: config
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev) ::   od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev) ::  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev)   ::   od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, nlev) :: &
        & ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1) :: &
         &  planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: ref_clear, trans_clear, reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: source_up_clear, source_dn_clear, source_up, source_dn

    ! Fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up, flux_dn
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up_clear, flux_dn_clear

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_lw) :: od_total, ssa_total, g_total

    ! Combined scattering optical depth
    real(jprb) :: scat_od, scat_od_total(config%n_g_lw)

    ! Two-stream coefficients
    real(jprb), dimension(config%n_g_lw) :: gamma1, gamma2

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_lw,nlev), intent(in) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_lw) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    real(jprb), intent(in) :: total_cloud_cover

    ! Identify clear-sky layers
    logical :: is_clear_sky_layer(nlev)

    ! Index of the highest cloudy layer
    integer :: i_cloud_top

    ! Number of g points
#ifndef _OPENACC
    integer :: ng
#else
    integer, parameter :: ng = JPGPT;
#endif

    ! Loop indices for level, column and g point
    integer :: jlev, jg
    real(jprb) :: acc

    real(jprb) :: hook_handle
    !$acc routine worker
#ifdef DR_HOOK
    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw_col',0,hook_handle)
#endif

#ifndef _OPENACC
    ng = config%n_g_lw
#endif

  ! Clear-sky calculation
  if (config%do_lw_aerosol_scattering) then
    ! Scattering case: first compute clear-sky reflectance,
    ! transmittance etc at each model level
#ifndef _OPENACC
    do jlev = 1,nlev
      ssa_total = ssa(:,jlev)
      g_total   = g(:,jlev)
      call calc_two_stream_gammas_lw(ng, ssa_total, g_total, &
           &  gamma1, gamma2)
      call calc_reflectance_transmittance_lw(ng, &
           &  od(:,jlev), gamma1, gamma2, &
           &  planck_hl(:,jlev), planck_hl(:,jlev+1), &
           &  ref_clear(:,jlev), trans_clear(:,jlev), &
           &  source_up_clear(:,jlev), source_dn_clear(:,jlev))
    end do
    ! Then use adding method to compute fluxes
    call adding_ica_lw(ng, nlev, &
         &  ref_clear, trans_clear, source_up_clear, source_dn_clear, &
         &  emission, albedo, &
         &  flux_up_clear, flux_dn_clear)
#endif
  else
    ! Non-scattering case: use simpler functions for
    ! transmission and emission
    !$acc loop seq
    do jlev = 1,nlev
      call calc_no_scattering_transmittance_lw(ng, od(:,jlev), &
           &  planck_hl(:,jlev), planck_hl(:,jlev+1), &
           &  trans_clear(:,jlev), source_up_clear(:,jlev), source_dn_clear(:,jlev))
    end do
    ! Simpler down-then-up method to compute fluxes
    call calc_fluxes_no_scattering_lw(ng, nlev, &
         &  trans_clear, source_up_clear, source_dn_clear, &
         &  emission, albedo, &
         &  flux_up_clear, flux_dn_clear)

    ! Ensure that clear-sky reflectance is zero since it may be
    ! used in cloudy-sky case
    ref_clear = 0.0_jprb
  end if

  ! Sum over g-points to compute broadband fluxes
  !$acc loop seq
  do jlev = 1,nlev+1
    flux%lw_up_clear(jcol,jlev) = sum_reduction(ng, flux_up_clear(:, jlev))
    flux%lw_dn_clear(jcol,jlev) = sum_reduction(ng, flux_dn_clear(:, jlev))
  end do

  ! Store surface spectral downwelling fluxes
  flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(:,nlev+1)

#define _REMOVE

  if (total_cloud_cover >= config%cloud_fraction_threshold) then
    ! Total-sky calculation

    is_clear_sky_layer = .true.
    i_cloud_top = nlev+1
    !$acc loop seq
    do jlev = 1,nlev
      ! Compute combined gas+aerosol+cloud optical properties
      if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
        is_clear_sky_layer(jlev) = .false.
        ! Get index to the first cloudy layer from the top
        if (i_cloud_top > jlev) then
          i_cloud_top = jlev
        end if

        od_cloud_new = od_scaling(:,jlev) &
             &  * od_cloud(config%i_band_from_reordered_g_lw,jlev)
        od_total = od(:,jlev) + od_cloud_new
        ssa_total = 0.0_jprb
        g_total   = 0.0_jprb
        if (config%do_lw_cloud_scattering) then
          ! Scattering case: calculate reflectance and
          ! transmittance at each model level

          if (config%do_lw_aerosol_scattering) then
            ! In single precision we need to protect against the
            ! case that od_total > 0.0 and ssa_total > 0.0 but
            ! od_total*ssa_total == 0 due to underflow
#ifndef _OPENACC
            scat_od_total = ssa(:,jlev)*od(:,jlev) &
                 &     + ssa_cloud(config%i_band_from_reordered_g_lw,jlev) &
                 &     *  od_cloud_new
            where (scat_od_total > 0.0_jprb)
              g_total = (g(:,jlev)*ssa(:,jlev)*od(:,jlev) &
                   &     +   g_cloud(config%i_band_from_reordered_g_lw,jlev) &
                   &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev) &
                   &     *  od_cloud_new) &
                   &     / scat_od_total
            end where
            where (od_total > 0.0_jprb)
              ssa_total = scat_od_total / od_total
            end where
#endif
          else
            !$acc loop independent worker private(scat_od)
            do jg = 1,ng
              if (od_total(jg) > 0.0_jprb) then
                scat_od = ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev) &
                     &     * od_cloud_new(jg)
                ssa_total(jg) = scat_od / od_total(jg)
                if (scat_od > 0.0_jprb) then
                  g_total(jg) = g_cloud(config%i_band_from_reordered_g_lw(jg),jlev) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev) &
                       &     *  od_cloud_new(jg) / scat_od
                end if
              end if
            end do
          end if
          ! Compute cloudy-sky reflectance, transmittance etc at
          ! each model level
          call calc_two_stream_gammas_lw(ng, ssa_total, g_total, &
               &  gamma1, gamma2)
          call calc_reflectance_transmittance_lw(ng, &
               &  od_total, gamma1, gamma2, &
               &  planck_hl(:,jlev), planck_hl(:,jlev+1), &
               &  reflectance(:,jlev), transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))
        else
#ifndef _OPENACC
          ! No-scattering case: use simpler functions for
          ! transmission and emission
          call calc_no_scattering_transmittance_lw(ng, od_total, &
               &  planck_hl(:,jlev), planck_hl(:,jlev+1), &
               &  transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))
#endif
        end if
      else
        ! Clear-sky layer: copy over clear-sky values
        reflectance(:,jlev) = ref_clear(:,jlev)
        transmittance(:,jlev) = trans_clear(:,jlev)
        source_up(:,jlev) = source_up_clear(:,jlev)
        source_dn(:,jlev) = source_dn_clear(:,jlev)
      end if
    end do


    if (config%do_lw_aerosol_scattering) then
      ! Use adding method to compute fluxes for an overcast sky,
      ! allowing for scattering in all layers
#ifndef _OPENACC
      call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
           &  emission, albedo, &
           &  flux_up, flux_dn)
#endif
    else if (config%do_lw_cloud_scattering) then
      ! Use adding method to compute fluxes but optimize for the
      ! presence of clear-sky layers
!          call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
!               &  emission(:,jcol), albedo(:,jcol), &
!               &  flux_up, flux_dn)
      call fast_adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
           &  emission, albedo, &
           &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
           &  flux_up, flux_dn)
    else
#ifndef _OPENACC
      ! Simpler down-then-up method to compute fluxes
      call calc_fluxes_no_scattering_lw(ng, nlev, &
           &  transmittance, source_up, source_dn, emission, albedo, &
           &  flux_up, flux_dn)
#endif
    end if

    ! Store overcast broadband fluxes

    do jlev = 1,nlev+1
        flux%lw_up(jcol,jlev) = sum_reduction(ng, flux_up(:, jlev))
        flux%lw_dn(jcol,jlev) = sum_reduction(ng, flux_dn(:, jlev))
    end do
    ! Cloudy flux profiles currently assume completely overcast
    ! skies; perform weighted average with clear-sky profile
    flux%lw_up(jcol,:) =  total_cloud_cover *flux%lw_up(jcol,:) &
         &  + (1.0_jprb - total_cloud_cover)*flux%lw_up_clear(jcol,:)
    flux%lw_dn(jcol,:) =  total_cloud_cover *flux%lw_dn(jcol,:) &
         &  + (1.0_jprb - total_cloud_cover)*flux%lw_dn_clear(jcol,:)
    ! Store surface spectral downwelling fluxes
    flux%lw_dn_surf_g(:,jcol) = total_cloud_cover*flux_dn(:,nlev+1) &
         &  + (1.0_jprb - total_cloud_cover)*flux%lw_dn_surf_clear_g(:,jcol)

#ifndef _OPENACC
    ! Compute the longwave derivatives needed by Hogan and Bozzo
    ! (2015) approximate radiation update scheme

    if (config%do_lw_derivatives) then
      call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance, flux_up(:,nlev+1), &
           &                       flux%lw_derivatives)
      if (total_cloud_cover < 1.0_jprb - config%cloud_fraction_threshold) then
        ! Modify the existing derivative with the contribution from the clear sky
        call modify_lw_derivatives_ica(ng, nlev, jcol, trans_clear, flux_up_clear(:,nlev+1), &
             &                         1.0_jprb-total_cloud_cover, flux%lw_derivatives)
      end if
    end if
#endif
  else
    ! No cloud in profile and clear-sky fluxes already
    ! calculated: copy them over
    flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
    flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
    flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)
#ifndef _OPENACC
    if (config%do_lw_derivatives) then
      call calc_lw_derivatives_ica(ng, nlev, jcol, trans_clear, flux_up_clear(:,nlev+1), &
           &                       flux%lw_derivatives)
    end if
#endif
  end if ! Cloud is present in profile

#ifdef DR_HOOK
    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw_col',1,hook_handle)
#endif
  end subroutine solver_mcica_lw_col

  pure function sum_reduction(ng, arr) result(res)
    use parkind1, only  : jprb
    USE YOERRTM  , ONLY : JPGPT
    integer, intent(in) :: ng
    real(jprb), intent(in) :: arr(1:ng)
    real(jprb) :: res
    !$acc routine worker

    res = sum(arr)

  end function

#ifdef _OPENACC
  subroutine init_gpu()
    use :: iso_c_binding
    use cudafor
    use openacc
    logical, save :: device_initialized = .False.
    integer, parameter :: KB = 1024
    integer, parameter :: MB = 1024 * KB
    integer, parameter :: DEV_NUM = 0 ! Always pick up first found device
    integer :: num_dev, istat
    character*100 :: dev_name
    integer(kind=cuda_count_kind) :: total_gpu_mem, heap_limit
    integer :: cuda_status
    integer(kind=c_int64_t) :: gpu_mem_size, gpu_mem_size_mb

    !$omp critical
    if (.not. device_initialized) then
      num_dev = acc_get_num_devices(ACC_DEVICE_NVIDIA)
      if (num_dev .lt. 1) then
        print *, 'Error: There are no Nvidia GPU devices available on this host'
        stop
      endif
      write(*,"(' Number of Nvidia GPU devices: ',i0)") num_dev
      call acc_get_property_string(DEV_NUM, ACC_DEVICE_NVIDIA, ACC_PROPERTY_NAME, dev_name)
      write(*,"(' Selecting device:             ',i0,' ""',a,'""')") DEV_NUM, dev_name
      if (acc_get_device_num(ACC_DEVICE_NVIDIA) /= DEV_NUM) then
        call acc_set_device_num(DEV_NUM, acc_device_nvidia)
      end if
      gpu_mem_size = acc_get_property(DEV_NUM, ACC_DEVICE_NVIDIA, ACC_PROPERTY_MEMORY)
      gpu_mem_size_mb = gpu_mem_size / (MB)
      write(*,"(' Device memory size:           ',i0,'mb')") gpu_mem_size_mb
      device_initialized = .True.
      istat = cudaDeviceGetLimit(heap_limit, cudaLimitMallocHeapSize)
      write(*,"(' Default CUDA heap size limit: ',i0,'mb')") (heap_limit / MB)
      istat = cudaDeviceSetLimit(cudaLimitMallocHeapSize, gpu_mem_size / 2)
      istat = cudaDeviceGetLimit(heap_limit, cudaLimitMallocHeapSize)
      write(*,"(' New CUDA heap size limit:     ',i0,'mb')") (heap_limit / MB)
    else
      if (acc_get_device_num(ACC_DEVICE_NVIDIA) /= DEV_NUM) then
        call acc_set_device_num(DEV_NUM, ACC_DEVICE_NVIDIA)
      end if
    end if
    !$omp end critical
  end subroutine
#endif

end module radiation_mcica_lw
