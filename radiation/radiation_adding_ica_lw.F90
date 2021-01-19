! radiation_adding_ica_lw.F90 - Longwave adding method in independent column approximation
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
!   2017-07-12  R. Hogan  Fast adding method for if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_adding_ica_lw

  public

contains

  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering, by successively adding the contribution of
  ! layers starting from the surface to compute the total albedo and
  ! total upward emission of the increasingly larger block of
  ! atmospheric layers.
  subroutine adding_ica_lw(ncol, nlev, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  flux_up, flux_dn)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(ncol, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to emission below
    ! that half-level (W m-2)
    real(jprb), dimension(ncol, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(ncol, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jcol

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw',0,hook_handle)

    albedo(:,nlev+1) = albedo_surf

    ! At the surface, the source is thermal emission
    source(:,nlev+1) = emission_surf

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    do jlev = nlev,1,-1
      ! Next loop over columns. We could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  Rather, we do it with an explicit
      ! loop.
      do jcol = 1,ncol
        ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
        inv_denominator(jcol,jlev) = 1.0_jprb &
             &  / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
        ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
        albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
             &  * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 11:
        source(jcol,jlev) = source_up(jcol,jlev) &
             &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
             &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev)) &
             &                   * inv_denominator(jcol,jlev)
      end do
    end do

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(:,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to emission
    ! below that level
    flux_up(:,1) = source(:,1)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
    do jlev = 1,nlev
      do jcol = 1,ncol
        ! Shonk & Hogan (2008) Eq 14 (after simplification):
        flux_dn(jcol,jlev+1) &
             &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
             &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
             &     + source_dn(jcol,jlev)) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 12:
        flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
             &            + source(jcol,jlev+1)
      end do
    end do

    if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw',1,hook_handle)

  end subroutine adding_ica_lw

    !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering, by successively adding the contribution of
  ! layers starting from the surface to compute the total albedo and
  ! total upward emission of the increasingly larger block of
  ! atmospheric layers.
  subroutine adding_ica_lw_lr(istartcol, iendcol, nlev, &
    &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
    &  flux_up, flux_dn)

 use parkind1, only           : jprb
 use yomhook,  only           : lhook, dr_hook

 implicit none

 ! Inputs
 integer, intent(in) :: istartcol, iendcol ! number of columns (may be spectral intervals)
 integer, intent(in) :: nlev ! number of levels

 ! Surface emission (W m-2) and albedo
 real(jprb), intent(in),  dimension(istartcol:iendcol) :: emission_surf, albedo_surf

 ! Diffuse reflectance and transmittance of each layer
 real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: reflectance, transmittance

 ! Emission from each layer in an upward and downward direction
 real(jprb), intent(in),  dimension(istartcol:iendcol, nlev)   :: source_up, source_dn

 ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
 ! downwelling
 real(jprb), intent(out), dimension(istartcol:iendcol,nlev+1) :: flux_up, flux_dn
 
 ! Albedo of the entire earth/atmosphere system below each half
 ! level
 real(jprb), dimension(nlev+1, istartcol:iendcol) :: albedo

 ! Upwelling radiation at each half-level due to emission below
 ! that half-level (W m-2)
 real(jprb), dimension(nlev+1, istartcol:iendcol) :: source

 ! Equal to 1/(1-albedo*reflectance)
 real(jprb), dimension(nlev, istartcol:iendcol)   :: inv_denominator

 ! Loop index for model level and column
 integer :: jlev, jcol

 real(jprb) :: hook_handle

 if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw_lr',0,hook_handle)

 albedo(nlev+1,:) = albedo_surf

 ! At the surface, the source is thermal emission
 source(nlev+1,:) = emission_surf

 ! Work back up through the atmosphere and compute the albedo of
 ! the entire earth/atmosphere system below that half-level, and
 ! also the "source", which is the upwelling flux due to emission
 ! below that level
 do jlev = nlev,1,-1
   ! Next loop over columns. We could do this by indexing the
   ! entire inner dimension as follows, e.g. for the first line:
   !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
   ! and similarly for subsequent lines, but this slows down the
   ! routine by a factor of 2!  Rather, we do it with an explicit
   ! loop.
   do jcol = istartcol,iendcol
     ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
     inv_denominator(jlev,jcol) = 1.0_jprb &
          &  / (1.0_jprb-albedo(jlev+1,jcol)*reflectance(jcol,jlev))
     ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
     albedo(jlev,jcol) = reflectance(jcol,jlev) + transmittance(jlev,jcol)*transmittance(jlev,jcol) &
          &  * albedo(jlev+1,jcol) * inv_denominator(jlev,jcol)
     ! Shonk & Hogan (2008) Eq 11:
     source(jlev,jcol) = source_up(jcol,jlev) &
          &  + transmittance(jlev,jcol) * (source(jlev+1,jcol) &
          &                    + albedo(jlev+1,jcol)*source_dn(jcol,jlev)) &
          &                   * inv_denominator(jlev,jcol)
   end do
 end do

 ! At top-of-atmosphere there is no diffuse downwelling radiation
 flux_dn(:,1) = 0.0_jprb

 ! At top-of-atmosphere, all upwelling radiation is due to emission
 ! below that level
 flux_up(:,1) = source(1,:)

 ! Work back down through the atmosphere computing the fluxes at
 ! each half-level
 do jlev = 1,nlev
   do jcol = istartcol,iendcol
     ! Shonk & Hogan (2008) Eq 14 (after simplification):
     flux_dn(jcol,jlev+1) &
          &  = (transmittance(jlev,jcol)*flux_dn(jcol,jlev) &
          &     + reflectance(jcol,jlev)*source(jlev+1,jcol) &
          &     + source_dn(jcol,jlev)) * inv_denominator(jlev,jcol)
     ! Shonk & Hogan (2008) Eq 12:
     flux_up(jcol,jlev+1) = albedo(jlev+1,jcol)*flux_dn(jcol,jlev+1) &
          &            + source(jlev+1,jcol)
   end do
 end do

 if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw_lr',1,hook_handle)

end subroutine adding_ica_lw_lr


subroutine adding_ica_lw_cond_lr(istartcol, iendcol, nlev, total_cloud_cover, cloud_fraction_threshold, &
  &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
  &  flux_up, flux_dn)

use parkind1, only           : jprb
use yomhook,  only           : lhook, dr_hook

implicit none

! Inputs
integer, intent(in) :: istartcol, iendcol ! number of columns (may be spectral intervals)
integer, intent(in) :: nlev ! number of levels
real(jprb), intent(in), dimension(istartcol:iendcol) :: total_cloud_cover
real(jprb), intent(in) :: cloud_fraction_threshold

! Surface emission (W m-2) and albedo
real(jprb), intent(in),  dimension(istartcol:iendcol) :: emission_surf, albedo_surf

! Diffuse reflectance and transmittance of each layer
real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: reflectance, transmittance

! Emission from each layer in an upward and downward direction
real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: source_up, source_dn

! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
! downwelling
real(jprb), intent(out), dimension(istartcol:iendcol, nlev+1) :: flux_up, flux_dn

! Albedo of the entire earth/atmosphere system below each half
! level
real(jprb), dimension(nlev+1, istartcol:iendcol) :: albedo

! Upwelling radiation at each half-level due to emission below
! that half-level (W m-2)
real(jprb), dimension(nlev+1, istartcol:iendcol) :: source

! Equal to 1/(1-albedo*reflectance)
real(jprb), dimension(nlev, istartcol:iendcol)   :: inv_denominator

! Loop index for model level and column
integer :: jlev, jcol

real(jprb) :: hook_handle

if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw_cond_lr',0,hook_handle)

do jcol = istartcol,iendcol
  if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then
    albedo(nlev+1,jcol) = albedo_surf(jcol)

    ! At the surface, the source is thermal emission
    source(nlev+1,jcol) = emission_surf(jcol)
  endif
enddo

! Work back up through the atmosphere and compute the albedo of
! the entire earth/atmosphere system below that half-level, and
! also the "source", which is the upwelling flux due to emission
! below that level
do jlev = nlev,1,-1
 ! Next loop over columns. We could do this by indexing the
 ! entire inner dimension as follows, e.g. for the first line:
 !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
 ! and similarly for subsequent lines, but this slows down the
 ! routine by a factor of 2!  Rather, we do it with an explicit
 ! loop.
 do jcol = istartcol,iendcol
  if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then
   ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
   inv_denominator(jlev,jcol) = 1.0_jprb &
        &  / (1.0_jprb-albedo(jlev+1,jcol)*reflectance(jcol,jlev))
   ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
   albedo(jlev,jcol) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
        &  * albedo(jlev+1,jcol) * inv_denominator(jlev,jcol)
   ! Shonk & Hogan (2008) Eq 11:
   source(jlev,jcol) = source_up(jcol,jlev) &
        &  + transmittance(jcol,jlev) * (source(jlev+1,jcol) &
        &                    + albedo(jlev+1,jcol)*source_dn(jcol,jlev)) &
        &                   * inv_denominator(jlev,jcol)

  endif
 end do
end do

do jcol = istartcol,iendcol
  if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then

  ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(jcol,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to emission
    ! below that level
    flux_up(jcol,1) = source(1,jcol)
  endif
enddo

! Work back down through the atmosphere computing the fluxes at
! each half-level
do jlev = 1,nlev
 do jcol = istartcol,iendcol
  if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then
     ! Shonk & Hogan (2008) Eq 14 (after simplification):
     flux_dn(jcol,jlev+1) &
        &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
        &     + reflectance(jcol,jlev)*source(jlev+1,jcol) &
        &     + source_dn(jcol,jlev)) * inv_denominator(jlev,jcol)
     ! Shonk & Hogan (2008) Eq 12:
     flux_up(jcol,jlev+1) = albedo(jlev+1,jcol)*flux_dn(jcol,jlev+1) &
        &            + source(jlev+1,jcol)
  endif
 end do
end do

if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw_cond_lr',1,hook_handle)

end subroutine adding_ica_lw_cond_lr


  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering in cloudy layers only.
  subroutine fast_adding_ica_lw(ncol, nlev, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
       &  flux_up, flux_dn)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! Determine which layers are cloud-free
    logical, intent(in) :: is_clear_sky_layer(nlev)

    ! Index to highest cloudy layer
    integer, intent(in) :: i_cloud_top

    ! Pre-computed clear-sky downwelling fluxes (W m-2) at half-levels
    real(jprb), intent(in), dimension(ncol, nlev+1)  :: flux_dn_clear

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(ncol, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to emission below
    ! that half-level (W m-2)
    real(jprb), dimension(ncol, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(ncol, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jcol

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_lw:fast_adding_ica_lw',0,hook_handle)

    ! Copy over downwelling fluxes above cloud from clear sky
    flux_dn(:,1:i_cloud_top) = flux_dn_clear(:,1:i_cloud_top)

    albedo(:,nlev+1) = albedo_surf
    
    ! At the surface, the source is thermal emission
    source(:,nlev+1) = emission_surf

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    do jlev = nlev,i_cloud_top,-1
      if (is_clear_sky_layer(jlev)) then
        ! Reflectance of this layer is zero, simplifying the expression
        do jcol = 1,ncol
          albedo(jcol,jlev) = transmittance(jcol,jlev)*transmittance(jcol,jlev)*albedo(jcol,jlev+1)
          source(jcol,jlev) = source_up(jcol,jlev) &
               &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
               &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev))
        end do
      else
        ! Loop over columns; explicit loop seems to be faster
        do jcol = 1,ncol
          ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
          inv_denominator(jcol,jlev) = 1.0_jprb &
               &  / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
          ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
          albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
               &  * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
          ! Shonk & Hogan (2008) Eq 11:
          source(jcol,jlev) = source_up(jcol,jlev) &
               &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
               &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev)) &
               &                   * inv_denominator(jcol,jlev)
        end do
      end if
    end do

    ! Compute the fluxes above the highest cloud
    flux_up(:,i_cloud_top) = source(:,i_cloud_top) &
         &                 + albedo(:,i_cloud_top)*flux_dn(:,i_cloud_top)
    do jlev = i_cloud_top-1,1,-1
      flux_up(:,jlev) = transmittance(:,jlev)*flux_up(:,jlev+1) + source_up(:,jlev)
    end do

    ! Work back down through the atmosphere from cloud top computing
    ! the fluxes at each half-level
    do jlev = i_cloud_top,nlev
      if (is_clear_sky_layer(jlev)) then
        do jcol = 1,ncol
          flux_dn(jcol,jlev+1) = transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
               &               + source_dn(jcol,jlev)
          flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
               &               + source(jcol,jlev+1)
        end do
      else
        do jcol = 1,ncol
          ! Shonk & Hogan (2008) Eq 14 (after simplification):
          flux_dn(jcol,jlev+1) &
               &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
               &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
               &     + source_dn(jcol,jlev)) * inv_denominator(jcol,jlev)
          ! Shonk & Hogan (2008) Eq 12:
          flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
               &               + source(jcol,jlev+1)
        end do
      end if
    end do

    if (lhook) call dr_hook('radiation_adding_ica_lw:fast_adding_ica_lw',1,hook_handle)

  end subroutine fast_adding_ica_lw

    !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering in cloudy layers only.
  subroutine fast_adding_ica_lw_lr(istartcol,iendcol, nlev, total_cloud_cover, cloud_fraction_threshold, &
    &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
    &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
    &  flux_up, flux_dn)

 use parkind1, only           : jprb
 use yomhook,  only           : lhook, dr_hook

 implicit none

 ! Inputs
 integer, intent(in) :: istartcol,iendcol ! number of columns (may be spectral intervals)
 integer, intent(in) :: nlev ! number of levels
 real(jprb), intent(in), dimension(istartcol:iendcol) :: total_cloud_cover
 real(jprb), intent(in) :: cloud_fraction_threshold

 ! Surface emission (W m-2) and albedo
 real(jprb), intent(in),  dimension(istartcol:iendcol) :: emission_surf, albedo_surf

 ! Diffuse reflectance and transmittance of each layer
 real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: reflectance, transmittance

 ! Emission from each layer in an upward and downward direction
 real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: source_up, source_dn

 ! Determine which layers are cloud-free
 logical, intent(in) :: is_clear_sky_layer(nlev,istartcol:iendcol)

 ! Index to highest cloudy layer
 integer, intent(in), dimension(istartcol:iendcol) :: i_cloud_top

 ! Pre-computed clear-sky downwelling fluxes (W m-2) at half-levels
 real(jprb), intent(in), dimension(istartcol:iendcol,nlev+1)  :: flux_dn_clear

 ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
 ! downwelling
 real(jprb), intent(out), dimension(istartcol:iendcol,nlev+1) :: flux_up, flux_dn
 
 ! Albedo of the entire earth/atmosphere system below each half
 ! level
 real(jprb), dimension(nlev+1,istartcol:iendcol) :: albedo

 ! Upwelling radiation at each half-level due to emission below
 ! that half-level (W m-2)
 real(jprb), dimension(nlev+1,istartcol:iendcol) :: source

 ! Equal to 1/(1-albedo*reflectance)
 real(jprb), dimension(nlev,istartcol:iendcol)   :: inv_denominator

 ! Loop index for model level and column
 integer :: jlev, jcol

 real(jprb) :: hook_handle

 if (lhook) call dr_hook('radiation_adding_ica_lw:fast_adding_ica_lw_lr',0,hook_handle)

 do jcol=istartcol,iendcol
  if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then
   ! Copy over downwelling fluxes above cloud from clear sky
     flux_dn(jcol,1:i_cloud_top(jcol)) = flux_dn_clear(jcol,1:i_cloud_top(jcol))

     albedo(nlev+1,jcol) = albedo_surf(jcol)
 
     ! At the surface, the source is thermal emission
     source(nlev+1,jcol) = emission_surf(jcol)
    
  endif
 enddo

 ! Work back up through the atmosphere and compute the albedo of
 ! the entire earth/atmosphere system below that half-level, and
 ! also the "source", which is the upwelling flux due to emission
 ! below that level
! do jlev = nlev,i_cloud_top,-1
  do jlev = nlev,1,-1
   do jcol = istartcol,iendcol
    if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then

     if(jlev >= i_cloud_top(jcol)) then
     
        if (is_clear_sky_layer(jlev,jcol)) then
        ! ! Reflectance of this layer is zero, simplifying the expression

          albedo(jlev,jcol) = transmittance(jcol,jlev)*transmittance(jcol,jlev)*albedo(jlev+1,jcol)
          source(jlev,jcol) = source_up(jcol,jlev) &
                &  + transmittance(jcol,jlev) * (source(jlev+1,jcol) &
                &                    + albedo(jlev+1,jcol)*source_dn(jcol,jlev))
        else
          ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
          inv_denominator(jlev,jcol) = 1.0_jprb &
                &  / (1.0_jprb-albedo(jlev+1,jcol)*reflectance(jcol,jlev))
          ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
          albedo(jlev,jcol) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
                &  * albedo(jlev+1,jcol) * inv_denominator(jlev,jcol)
          ! Shonk & Hogan (2008) Eq 11:
          source(jlev,jcol) = source_up(jcol,jlev) &
                &  + transmittance(jcol,jlev) * (source(jlev+1,jcol) &
                &                    + albedo(jlev+1,jcol)*source_dn(jcol,jlev)) &
                &                   * inv_denominator(jlev,jcol)
        endif
      endif
    endif
  end do
 end do
 
 do jcol=istartcol,iendcol
  if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then
   ! Compute the fluxes above the highest cloud
   flux_up(jcol,i_cloud_top(jcol)) = source(i_cloud_top(jcol),jcol) &
      &                 + albedo(i_cloud_top(jcol),jcol)*flux_dn(jcol,i_cloud_top(jcol))
  endif
 enddo
 !do jlev = i_cloud_top(jcol)-1,1,-1
 ! cos: would it make sense to revert the jcol and jlev loops here
 ! in this pattern to avoid the jcol conditional exit? performance wise need to test?
 do jlev = nlev-1,1,-1
   do jcol=istartcol,iendcol
    if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then
     if(jlev < i_cloud_top(jcol)) then
       flux_up(jcol,jlev) = transmittance(jcol,jlev)*flux_up(jcol,jlev+1) + source_up(jcol,jlev)
     endif
    endif
   end do
 enddo

 ! Work back down through the atmosphere from cloud top computing
 ! the fluxes at each half-level

 !do jlev = i_cloud_top,nlev
 do jlev = 1,nlev
  do jcol=istartcol,iendcol
    if ((total_cloud_cover(jcol) >= cloud_fraction_threshold) .and. jlev >= i_cloud_top(jcol)) then

   if (is_clear_sky_layer(jlev,jcol)) then
       flux_dn(jcol,jlev+1) = transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
            &               + source_dn(jcol,jlev)
       flux_up(jcol,jlev+1) = albedo(jlev+1,jcol)*flux_dn(jcol,jlev+1) &
            &               + source(jlev+1,jcol)
   else
       ! Shonk & Hogan (2008) Eq 14 (after simplification):
       flux_dn(jcol,jlev+1) &
            &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
            &     + reflectance(jcol,jlev)*source(jlev+1,jcol) &
            &     + source_dn(jcol,jlev)) * inv_denominator(jlev,jcol)
       ! Shonk & Hogan (2008) Eq 12:
       flux_up(jcol,jlev+1) = albedo(jlev+1,jcol)*flux_dn(jcol,jlev+1) &
            &               + source(jlev+1,jcol)
   end if
  endif
  enddo
 end do

 if (lhook) call dr_hook('radiation_adding_ica_lw:fast_adding_ica_lw_lr',1,hook_handle)

end subroutine fast_adding_ica_lw_lr


  !---------------------------------------------------------------------
  ! If there is no scattering then fluxes may be computed simply by
  ! passing down through the atmosphere computing the downwelling
  ! fluxes from the transmission and emission of each layer, and then
  ! passing back up through the atmosphere to compute the upwelling
  ! fluxes in the same way.
  subroutine calc_fluxes_no_scattering_lw(ncol, nlev, &
       &  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! Loop index for model level
    integer :: jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw',0,hook_handle)

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(:,1) = 0.0_jprb

    ! Work down through the atmosphere computing the downward fluxes
    ! at each half-level
    do jlev = 1,nlev
      flux_dn(:,jlev+1) = transmittance(:,jlev)*flux_dn(:,jlev) + source_dn(:,jlev)
    end do

    ! Surface reflection and emission
    flux_up(:,nlev+1) = emission_surf + albedo_surf * flux_dn(:,nlev+1)

    ! Work back up through the atmosphere computing the upward fluxes
    ! at each half-level
    do jlev = nlev,1,-1
      flux_up(:,jlev) = transmittance(:,jlev)*flux_up(:,jlev+1) + source_up(:,jlev)
    end do
    
    if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw',1,hook_handle)

  end subroutine calc_fluxes_no_scattering_lw

    !---------------------------------------------------------------------
  ! If there is no scattering then fluxes may be computed simply by
  ! passing down through the atmosphere computing the downwelling
  ! fluxes from the transmission and emission of each layer, and then
  ! passing back up through the atmosphere to compute the upwelling
  ! fluxes in the same way.
  subroutine calc_fluxes_no_scattering_lw_lr(istartcol, iendcol, nlev, &
    &  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)

 use parkind1, only           : jprb
 use yomhook,  only           : lhook, dr_hook

 implicit none

 ! Inputs
 integer, intent(in) :: istartcol, iendcol ! number of columns (may be spectral intervals)
 integer, intent(in) :: nlev ! number of levels

 ! Surface emission (W m-2) and albedo
 real(jprb), intent(in),  dimension(istartcol:iendcol) :: emission_surf, albedo_surf

 ! Diffuse reflectance and transmittance of each layer
 real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: transmittance

 ! Emission from each layer in an upward and downward direction
 real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: source_up, source_dn

 ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
 ! downwelling
 real(jprb), intent(out), dimension(istartcol:iendcol,nlev+1) :: flux_up, flux_dn
 
 ! Loop index for model level
 integer :: jlev

 real(jprb) :: hook_handle

 if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw_lr',0,hook_handle)

 ! At top-of-atmosphere there is no diffuse downwelling radiation
 flux_dn(:,1) = 0.0_jprb

 ! Work down through the atmosphere computing the downward fluxes
 ! at each half-level
 do jlev = 1,nlev
   flux_dn(:,jlev+1) = transmittance(:,jlev)*flux_dn(:,jlev) + source_dn(:,jlev)
 end do

 ! Surface reflection and emission
 flux_up(:,nlev+1) = emission_surf + albedo_surf * flux_dn(:,nlev+1)

 ! Work back up through the atmosphere computing the upward fluxes
 ! at each half-level
 do jlev = nlev,1,-1
   flux_up(:,jlev) = transmittance(:,jlev)*flux_up(:,jlev+1) + source_up(:,jlev)
 end do
 
 if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw_lr',1,hook_handle)

end subroutine calc_fluxes_no_scattering_lw_lr

subroutine calc_fluxes_no_scattering_lw_cond_lr(istartcol, iendcol, nlev, total_cloud_cover, cloud_fraction_threshold, &
  &  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)

use parkind1, only           : jprb
use yomhook,  only           : lhook, dr_hook

implicit none

! Inputs
integer, intent(in) :: istartcol, iendcol ! number of columns (may be spectral intervals)
integer, intent(in) :: nlev ! number of levels
real(jprb), intent(in), dimension(istartcol:iendcol) :: total_cloud_cover
real(jprb), intent(in) :: cloud_fraction_threshold

! Surface emission (W m-2) and albedo
real(jprb), intent(in),  dimension(istartcol:iendcol) :: emission_surf, albedo_surf

! Diffuse reflectance and transmittance of each layer
real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: transmittance

! Emission from each layer in an upward and downward direction
real(jprb), intent(in),  dimension(istartcol:iendcol,nlev)   :: source_up, source_dn

! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
! downwelling
real(jprb), intent(out), dimension(nlev+1,istartcol:iendcol) :: flux_up, flux_dn

! Loop index for model level
integer :: jlev, jcol

real(jprb) :: hook_handle

if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw_cond_lr',0,hook_handle)

do jcol=istartcol,iendcol
  if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then
! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(1,jcol) = 0.0_jprb
  endif
enddo

! Work down through the atmosphere computing the downward fluxes
! at each half-level
do jlev = 1,nlev
  do jcol=istartcol,iendcol
    if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then  
      flux_dn(jlev+1,jcol) = transmittance(jcol,jlev)*flux_dn(jlev,jcol) + source_dn(jcol,jlev)
    endif
  enddo
end do

do jcol=istartcol,iendcol
  if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then
! Surface reflection and emission
    flux_up(nlev+1,jcol) = emission_surf(jcol) + albedo_surf(jcol) * flux_dn(nlev+1,jcol)
  endif
enddo

! Work back up through the atmosphere computing the upward fluxes
! at each half-level
do jlev = nlev,1,-1
  do jcol=istartcol,iendcol
    if (total_cloud_cover(jcol) >= cloud_fraction_threshold) then 
      flux_up(jlev,jcol) = transmittance(jcol,jlev)*flux_up(jlev+1,jcol) + source_up(jcol,jlev)
    endif
  enddo
end do

if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw_cond_lr',1,hook_handle)

end subroutine calc_fluxes_no_scattering_lw_cond_lr


end module radiation_adding_ica_lw
