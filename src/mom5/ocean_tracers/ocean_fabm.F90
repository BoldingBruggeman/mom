!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include <fms_platform.h>

!
!
!<CONTACT EMAIL="jorn@bolding-bruggeman.com"> Jorn Bruggeman
!</CONTACT>
!
!<REVIEWER EMAIL="TODO"> TODO
!</REVIEWER>
!
!<OVERVIEW>
! MOM5 driver for the Framework for Aquatic Biogeochemical Models (FABM)
!</OVERVIEW>
!
!<DESCRIPTION>
!       TODO
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! TODO
! </REFERENCE>
!
! </INFO>
!
! $Id: ocean_fabm.F90 69 2010-09-07 15:00:18Z jornbr $
!

!
!------------------------------------------------------------------
!
!       Module ocean_fabm_mod
!
!       TODO
!
!------------------------------------------------------------------
!

module  ocean_fabm_mod  !{

!
!------------------------------------------------------------------
!
!       Global definitions
!
!------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Modules
!
!----------------------------------------------------------------------
!

use time_manager_mod,         only: time_type, get_date, set_date, operator(-), get_time
use diag_manager_mod,         only: send_data
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len, fm_get_index
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value
use fms_mod,                  only: write_data
use fms_io_mod,               only: restart_file_type, register_restart_field, restore_state, save_restart
use mpp_mod,                  only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_interp_external_mod, only: time_interp_external
use mpp_domains_mod,          only: domain2d, mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use data_override_mod,only: data_override

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer, otpm_set_diag_tracer
use fm_util_mod, only: fm_util_check_for_bad_fields, fm_util_set_value, fm_util_set_value_string_array
use fm_util_mod, only: fm_util_get_string, fm_util_get_integer_array, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod, only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod, only: fm_util_start_namelist, fm_util_end_namelist
!use fm_util_mod, only: grid, dtts
!use fm_util_mod, only: taum1, tau, taup1
!use fm_util_mod, only: t_prog, t_diag
use coupler_types_mod,  only: coupler_2d_bc_type ! ind_alpha, ind_csurf
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_diag_tracer_type, ocean_density_type

!use ocean_sbc_mod,      only: use_waterflux
!use ocean_types_mod,    only: ocean_thickness_type,ocean_density_type

use fabm
use fabm_types
use fabm_config
use fabm_driver

!
!----------------------------------------------------------------------
!
!       force all variables to be "typed"
!
!----------------------------------------------------------------------
!

implicit none

!
!----------------------------------------------------------------------
!
!       Make all routines and variables private by default
!
!----------------------------------------------------------------------
!

private

!
!----------------------------------------------------------------------
!
!       Public routines
!
!----------------------------------------------------------------------
!

public  :: ocean_fabm_bbc
public  :: ocean_fabm_end
public  :: ocean_fabm_init
public  :: ocean_fabm_sbc
public  :: ocean_fabm_source
public  :: ocean_fabm_start
public  :: ocean_fabm_tracer
public  :: ocean_fabm_avg_sfc
public  :: ocean_fabm_sum_sfc
public  :: ocean_fabm_zero_sfc
public  :: ocean_fabm_flux_init
public  :: ocean_fabm_sfc_end
public  :: ocean_fabm_init_sfc

!
!----------------------------------------------------------------------
!
!       Private parameters
!
!----------------------------------------------------------------------
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_fabm'
character(len=48), parameter                    :: mod_name = 'ocean_fabm_mod'
character(len=fm_string_len), parameter         :: default_restart_file = 'ocean_fabm.res.nc'
character(len=fm_string_len), parameter         :: default_2d_restart_file = 'ocean_fabm_2d.res.nc'
character(len=fm_string_len), parameter         :: default_ice_restart_file = 'ice_ocean_fabm.res.nc'
character(len=fm_string_len), parameter         :: default_ocean_restart_file = 'ocean_fabm_airsea_flux.res.nc'

!
!----------------------------------------------------------------------
!
!       Private types
!
!----------------------------------------------------------------------
!

type,extends(type_base_driver) :: type_mom_driver
contains
   procedure :: mom_driver_fatal_error
   procedure :: mom_driver_log_message
end type

type type_forcing_2d
   character(len=128)             :: name =''
   real, _ALLOCATABLE             :: data(:,:) _NULL
   type(type_forcing_2d), pointer :: next => null()
end type

type type_forcing_3d
   character(len=128)             :: name =''
   real, _ALLOCATABLE             :: data(:,:,:) _NULL
   type(type_forcing_3d), pointer :: next => null()
end type

!
!----------------------------------------------------------------------
!
!       Public variables
!
!----------------------------------------------------------------------
!

logical, public :: do_ocean_fabm

!
!----------------------------------------------------------------------
!
!       Private variables
!
!----------------------------------------------------------------------
!

type (type_model) :: model

! MOM-FABM namelist parameters
logical :: inplace_repair
logical :: zero_river_concentration
logical :: disable_sources
logical :: disable_vertical_movement

! Identifiers for environmental depenendeices provided by MOM
type (type_bulk_variable_id) :: id_temp
type (type_bulk_variable_id) :: id_salt
type (type_bulk_variable_id) :: id_pres
type (type_bulk_variable_id) :: id_par
type (type_bulk_variable_id) :: id_dens

! Indices of prognostic variables (i.e., index in T_prog for each interior state variable)
integer,_ALLOCATABLE,dimension(:) :: interior_state_indices

! Indices of diagnostic outputs
integer,_ALLOCATABLE,dimension(:) :: interior_diagnostic_indices
integer,_ALLOCATABLE,dimension(:) :: horizontal_diagnostic_indices
integer,_ALLOCATABLE,dimension(:) :: bottom_state_indices
integer,_ALLOCATABLE,dimension(:) :: surface_state_indices
integer,_ALLOCATABLE,dimension(:) :: inds_clip
integer,_ALLOCATABLE,dimension(:) :: inds_cons_tot
integer,_ALLOCATABLE,dimension(:) :: inds_cons_ave

! Array for surface/bottom attached state variables [not stored by MOM]
real, _ALLOCATABLE :: bottom_state(:,:,:) _NULL
real, _ALLOCATABLE :: surface_state(:,:,:) _NULL

! Arrays to hold information on externally provided fields (specified by the user in data_table)
type(type_forcing_2d), pointer,save :: first_forcing_2d => null()
type(type_forcing_3d), pointer,save :: first_forcing_3d => null()

! Work arrays
real,_ALLOCATABLE,dimension(:,:)   :: work_cons _NULL
real,_ALLOCATABLE,dimension(:,:,:) :: w _NULL
real,_ALLOCATABLE,dimension(:,:)   :: adv _NULL,work_dy _NULL,work_dy_sf _NULL,work_dy_bt _NULL
real, _ALLOCATABLE :: surface_fluxes(:,:,:) _NULL
real, _ALLOCATABLE :: bottom_fluxes(:,:,:) _NULL
real,_ALLOCATABLE,dimension(:,:,:) :: clipped _NULL

! Environmental dependencies computed by the MOM-FABM coupler
real, target :: days_since_start_of_the_year

! Information on additonal restarts (surface/bottom attached state)
type(restart_file_type) :: restart_2d

integer      :: package_index
logical,save :: module_initialized = .false.
logical,save :: initialization_complete = .false.

! Indices of variable maintained by MOM
integer :: indsal
integer :: indtemp
integer :: index_irr
integer :: index_chl

! Indices of additional timers
integer :: id_clock_fabm_setup_env
integer :: id_clock_fabm_repair
integer :: id_clock_fabm_conservation
integer :: id_clock_fabm_interior
integer :: id_clock_fabm_surface
integer :: id_clock_fabm_bottom
integer :: id_clock_fabm_vert_mov

!
!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------
!

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>

subroutine ocean_fabm_flux_init  !{

   use fms_mod, only : close_file,file_exist
   use mpp_io_mod, only: mpp_open,MPP_RDONLY

   character(len=64), parameter    :: sub_name = 'ocean_fabm_flux_init'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   character(len=256) :: configuration_file
   integer            :: yaml_unit
   character(len=256) :: caller_str

   ! First, perform some initialization if this module has not been
   ! initialized because the normal initialization routine will
   ! not have been called as part of the normal ocean model
   ! initialization if this is an Atmosphere pe of a coupled
   ! model running in concurrent mode
   if (.not. module_initialized) then

      ! Connect custom logging/error reporting object to FABM.
      allocate(type_mom_driver::driver)

      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

      package_index = otpm_set_tracer_package(package_name,            &
         restart_file = default_restart_file,                        &
         caller = caller_str)

      call fm_util_start_namelist('', package_name, caller = caller_str, no_overwrite = .true., &
            check = .true.)
      call fm_util_set_value('use', .false.)
      call fm_util_set_value('configuration_file', 'fabm.yaml')

   !  call fm_util_start_namelist(package_name, name, caller = caller_str)
      do_ocean_fabm = fm_util_get_logical('use', caller = caller_str, scalar = .true.)
      if(.not. do_ocean_fabm) return 

      ! Read FABM configuration from fabm.yaml.
      configuration_file = fm_util_get_string ('configuration_file', caller = caller_str, scalar = .true.)
      call mpp_open(yaml_unit, configuration_file, action=MPP_RDONLY )
      call fabm_create_model_from_yaml_file(model,unit=yaml_unit)
      call close_file (yaml_unit)
   endif

end subroutine  ocean_fabm_flux_init  !}
!</SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocean_fabm_init  !{

   use fms_mod, only : close_file,file_exist
   use mpp_io_mod, only: mpp_open,MPP_RDONLY
   use diag_manager_mod, only: register_diag_field

   character(len=64), parameter    :: sub_name = 'ocean_fabm_init'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   character(len=256)                                      :: caller_str
   character(len=fm_string_len), pointer, dimension(:)     :: good_list

   integer            :: yaml_unit
   integer            :: ivar
   character(len=256) :: configuration_file

   if (module_initialized) return 

   caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

   ! Initialize the FABM package
   package_index = otpm_set_tracer_package(package_name,            &
      restart_file = default_restart_file,     &
      caller = caller_str)

   ! Connect custom logging/error reporting object to FABM.
   allocate(type_mom_driver::driver)

   ! Add the package name to the list of good namelists, to be used
   ! later for a consistency check
   if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) <= 0) then  !{
     call mpp_error(FATAL, trim(error_header) //                           &
          ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
   endif  !}

   call fm_util_start_namelist('', package_name, caller = caller_str, no_overwrite = .true., &
        check = .true.)

   call fm_util_set_value('use', .false.)
   call fm_util_set_value('configuration_file', 'fabm.yaml')
   call fm_util_set_value('inplace_repair', .false.)
   call fm_util_set_value('zero_river_concentration', .false.)
   call fm_util_set_value('disable_sources', .false.)
   call fm_util_set_value('disable_vertical_movement', .false.)

   do_ocean_fabm = fm_util_get_logical('use', caller = caller_str, scalar = .true.)

   if (.not.do_ocean_fabm) return

   ! Read FABM configuration from fabm.yaml.
   configuration_file = fm_util_get_string('configuration_file', caller = caller_str, scalar = .true.)
   call mpp_open(yaml_unit, configuration_file, action=MPP_RDONLY )
   call fabm_create_model_from_yaml_file(model,unit=yaml_unit)
   call close_file(yaml_unit)

   ! Register interior biogeochemical state variables
   allocate(interior_state_indices(size(model%state_variables)))
   do ivar=1,size(model%state_variables)
      interior_state_indices(ivar) = otpm_set_prog_tracer(           &
         trim(model%state_variables(ivar)%name),                     &
         package_name,                                               &
         longname = trim(model%state_variables(ivar)%long_name),     &
         units = trim(model%state_variables(ivar)%units)//' m3 kg-1',&
         flux_units = trim(model%state_variables(ivar)%units)//'/s', &
         caller = trim(mod_name)//'('//trim(sub_name)//')',          &
         min_tracer_limit = model%state_variables(ivar)%minimum,     &
         max_tracer_limit = model%state_variables(ivar)%maximum,     &
         !min_tracer = model%state_variables(ivar)%minimum,          &
         !max_tracer = model%state_variables(ivar)%maximum,          &
         const_init_tracer = .true.,                                 &
         const_init_value = model%state_variables(ivar)%initial_value)
   end do

!  call fm_util_end_namelist(package_name, '*global*', caller = caller_str)
   call fm_util_end_namelist('', package_name, check = .true., caller = caller_str)

   ! Obtain ids of environmental variables that MOM can provide
   id_temp = fabm_get_bulk_variable_id(model,standard_variables%temperature)
   id_salt = fabm_get_bulk_variable_id(model,standard_variables%practical_salinity)
   id_pres = fabm_get_bulk_variable_id(model,standard_variables%pressure)
   id_dens = fabm_get_bulk_variable_id(model,standard_variables%density)
   id_par  = fabm_get_bulk_variable_id(model,standard_variables%downwelling_photosynthetic_radiative_flux)

!
!       Check for any errors in the number of fields in the namelists for this package
!

   !good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
   !     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
   !if (associated(good_list)) then  !{
   !  call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
   !       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
   !  deallocate(good_list)
   !else  !}{
   !  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
   !end if  !}

   ! Create a diagnostic for chlorophyll
   ! This will be used to compute shortwave atentuation if read_chl=.false. in the shortwave namelist
   index_chl = otpm_set_diag_tracer('chl',                                           &
    caller=trim(mod_name)//'('//trim(sub_name)//')',                                 &
    longname='chlorophyll', units='mg/m^3',                                          &
    conversion=1.0, offset=0.0, min_tracer=0.0, max_tracer=1.e20,                    &
    min_range=-10.0, max_range=100.0, const_init_tracer=.true.,const_init_value=0.0, &
    restart_file='ocean_chl.res.nc')
   if (index_chl==-1) call mpp_error(FATAL,trim(error_header) // 'index of chlorophyll equals -1.')

   id_clock_fabm_setup_env    = mpp_clock_id('(Ocean FABM: setup environment) ',  grain=CLOCK_ROUTINE)
   id_clock_fabm_repair       = mpp_clock_id('(Ocean FABM: repair state) ',       grain=CLOCK_ROUTINE)
   id_clock_fabm_conservation = mpp_clock_id('(Ocean FABM: conserved quantities)',grain=CLOCK_ROUTINE)
   id_clock_fabm_interior     = mpp_clock_id('(Ocean FABM: interior) ',           grain=CLOCK_ROUTINE)
   id_clock_fabm_surface      = mpp_clock_id('(Ocean FABM: surface) ',            grain=CLOCK_ROUTINE)
   id_clock_fabm_bottom       = mpp_clock_id('(Ocean FABM: bottom) ',             grain=CLOCK_ROUTINE)
   id_clock_fabm_vert_mov     = mpp_clock_id('(Ocean FABM: vertical movement) ',  grain=CLOCK_ROUTINE)

   module_initialized = .true.

end subroutine ocean_fabm_init  !}
! </SUBROUTINE> NAME="ocean_fabm_init"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_fabm_zero_sfc(Ocean_fields)  !{

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

end subroutine ocean_fabm_zero_sfc  !}
! </SUBROUTINE> NAME="ocean_fabm_zero_sfc"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_fabm_sum_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     Ocean_fields, T_prog, rho, taum1, model_time, grid_tmask)  !{

   integer, intent(in)                                     :: isc
   integer, intent(in)                                     :: iec
   integer, intent(in)                                     :: jsc
   integer, intent(in)                                     :: jec
   integer, intent(in)                                     :: nk
   integer, intent(in)                                     :: isd
   integer, intent(in)                                     :: ied
   integer, intent(in)                                     :: jsd
   integer, intent(in)                                     :: jed
   integer, intent(in)                                     :: isc_bnd
   integer, intent(in)                                     :: iec_bnd
   integer, intent(in)                                     :: jsc_bnd
   integer, intent(in)                                     :: jec_bnd
   type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
   type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
   real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
   integer, intent(in)                                     :: taum1
   type(time_type), intent(in)                             :: model_time
   real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

end subroutine ocean_fabm_sum_sfc  !}
! </SUBROUTINE> NAME="ocean_fabm_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocean_fabm_sbc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     T_prog, taum1, model_time, grid_tmask, ice_ocean_boundary_fluxes)  !{

   use coupler_types_mod, only       : coupler_2d_bc_type !, ind_flux
   use mpp_mod, only                 : mpp_sum

   integer, intent(in)                                             :: isc
   integer, intent(in)                                             :: iec
   integer, intent(in)                                             :: jsc
   integer, intent(in)                                             :: jec
   integer, intent(in)                                             :: nk
   integer, intent(in)                                             :: isd
   integer, intent(in)                                             :: ied
   integer, intent(in)                                             :: jsd
   integer, intent(in)                                             :: jed
   integer, intent(in)                                             :: isc_bnd
   integer, intent(in)                                             :: iec_bnd
   integer, intent(in)                                             :: jsc_bnd
   integer, intent(in)                                             :: jec_bnd
   type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
   integer, intent(in)                                             :: taum1
   type(time_type), intent(in)                                     :: model_time
   real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask
   type(coupler_2d_bc_type), intent(in)                            :: ice_ocean_boundary_fluxes

   integer :: ivar

   do ivar=1,size(model%state_variables)
      if (model%state_variables(ivar)%no_precipitation_dilution) &
         T_prog(interior_state_indices(ivar))%tpme  (isc:iec,jsc:jec) = T_prog(interior_state_indices(ivar))%field(isc:iec,jsc:jec,1,taum1)
      if ((.not.zero_river_concentration).and.model%state_variables(ivar)%no_river_dilution) &
         T_prog(interior_state_indices(ivar))%triver(isc:iec,jsc:jec) = T_prog(interior_state_indices(ivar))%field(isc:iec,jsc:jec,1,taum1)

      ! stf is defined as bottom tracer flux [rho*m/sec*tracer concen] in ocean_types_mod.
      ! Since MOM defines tracer concentration as tracer quantity per seawater mass [not volume!],
      ! multiplication with rho (kg m-3) ensures it becomes tracer quantity per volume, as FABM already uses internally.
      ! Therefore, we can directly copy the FABM's flux values.
      t_prog(interior_state_indices(ivar))%stf(isc:iec,jsc:jec) = surface_fluxes(isc:iec,jsc:jec,ivar)
   end do

end subroutine  ocean_fabm_sbc  !}
! </SUBROUTINE> NAME="ocean_fabm_sbc"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!

subroutine ocean_fabm_start(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,      &
     T_prog, T_diag, taup1, model_time, grid_dat, grid_tmask, grid_kmt,                 &
     grid_xt, grid_yt, grid_zt, grid_zw, grid_ht, grid_dzt, grid_name, grid_tracer_axes, &
     mpp_domain2d, rho_dzt, dzt, swflx, swflx_vis, current_wave_stress)  !{

   use time_interp_external_mod, only: init_external_field
   use diag_manager_mod, only        : register_diag_field, diag_axis_init
   use field_manager_mod, only       : fm_get_index

   integer, intent(in)                                     :: isc
   integer, intent(in)                                     :: iec
   integer, intent(in)                                     :: jsc
   integer, intent(in)                                     :: jec
   integer, intent(in)                                     :: nk
   integer, intent(in)                                     :: isd
   integer, intent(in)                                     :: ied
   integer, intent(in)                                     :: jsd
   integer, intent(in)                                     :: jed
   type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
   type(ocean_diag_tracer_type), intent(in), dimension(:)  :: T_diag
   integer, intent(in)                                     :: taup1
   type(time_type), intent(in)                             :: model_time
   real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
   real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
   integer, dimension(isd:,jsd:), intent(in)               :: grid_kmt
   real, dimension(isd:,jsd:), intent(in), target          :: grid_xt
   real, dimension(isd:,jsd:), intent(in), target          :: grid_yt
   real, dimension(:), intent(in)                          :: grid_zt
   real, dimension(:), intent(in)                          :: grid_zw
   real, dimension(isd:,jsd:), intent(in)                  :: grid_ht
   real, dimension(:), intent(in)                          :: grid_dzt
   character(len=*), intent(in)                            :: grid_name
   integer, dimension(3), intent(in)                       :: grid_tracer_axes
   type(domain2d), intent(in)                              :: mpp_domain2d
   real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt
   real, dimension(isd:,jsd:,:), target                    :: dzt
   real, dimension(isd:,jsd:), target                      :: swflx, swflx_vis, current_wave_stress

   character(len=64), parameter    :: sub_name = 'ocean_fabm_start'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   character(len=256) :: caller_str
   integer            :: ivar
   integer            :: id_restart

   write(stdout(),*)
   write(stdout(),*) trim(note_header),                     &
                     'Starting ', trim(package_name), ' module'

   indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
   if (indtemp <= 0) then  !{
     call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
   endif  !}

   indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
   if (indsal <= 0) then  !{
     call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
   endif  !}

   index_irr = fm_get_index('/ocean_mod/diag_tracers/irr')
   if (index_irr==-1) call mpp_error(FATAL,trim(error_header) // 'index of shortwave irradiance equals -1.')

   caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

   call fm_util_start_namelist('', package_name, caller = caller_str)

   inplace_repair            = fm_util_get_logical('inplace_repair',           scalar = .true.)
   zero_river_concentration  = fm_util_get_logical('zero_river_concentration', scalar = .true.)
   disable_sources           = fm_util_get_logical('disable_sources',          scalar = .true.)
   disable_vertical_movement = fm_util_get_logical('disable_vertical_movement',scalar = .true.)

   call fm_util_end_namelist('', package_name, caller = caller_str)

   ! Register diagnostic variables defined on interior model domain.
   allocate(interior_diagnostic_indices(size(model%diagnostic_variables)))
   do ivar=1,size(model%diagnostic_variables)
      interior_diagnostic_indices(ivar) = register_diag_field('ocean_model',      &
         trim(model%diagnostic_variables(ivar)%name), grid_tracer_axes(1:3),                       &
         model_time, trim(model%diagnostic_variables(ivar)%long_name), &
         trim(model%diagnostic_variables(ivar)%units),            &
         missing_value = model%diagnostic_variables(ivar)%missing_value)
      model%diagnostic_variables(ivar)%save = interior_diagnostic_indices(ivar)/=-1
   end do

   ! Register diagnostic variables defined on horizontal slice of domain.
   allocate(horizontal_diagnostic_indices(size(model%horizontal_diagnostic_variables)))
   do ivar=1,size(model%horizontal_diagnostic_variables)
      horizontal_diagnostic_indices(ivar) = register_diag_field('ocean_model',      &
         trim(model%horizontal_diagnostic_variables(ivar)%name), grid_tracer_axes(1:2),                       &
         model_time, trim(model%horizontal_diagnostic_variables(ivar)%long_name), &
         trim(model%horizontal_diagnostic_variables(ivar)%units),            &
         missing_value = model%horizontal_diagnostic_variables(ivar)%missing_value)
      model%horizontal_diagnostic_variables(ivar)%save = horizontal_diagnostic_indices(ivar)/=-1
   end do

   ! Register clipping diagnostic (increase/time) for interior state variables.
   allocate(inds_clip(size(model%state_variables)))
   do ivar=1,size(model%state_variables)
      inds_clip(ivar) = register_diag_field('ocean_model',      &
         trim(model%state_variables(ivar)%name)//'_clip', grid_tracer_axes(1:3),                       &
         model_time, trim(model%state_variables(ivar)%long_name)//' clipping increase', &
         trim(model%state_variables(ivar)%units),            &
         missing_value = -1.0e+10)
   end do

   ! Register local values and domain-wide integrals of conserved quantities.
   allocate(inds_cons_tot(size(model%conserved_quantities)))
   allocate(inds_cons_ave(size(model%conserved_quantities)))
   do ivar=1,size(model%conserved_quantities)
      inds_cons_tot(ivar) = register_diag_field('ocean_model',                               &
         trim(model%conserved_quantities(ivar)%name)//'_global_int',                         &
         model_time, 'global integrated '//trim(model%conserved_quantities(ivar)%long_name), &
         trim(model%conserved_quantities(ivar)%units)//'*m3/1e21',                           &
         missing_value = -1.0e+10)
      inds_cons_ave(ivar) = register_diag_field('ocean_model',                            &
         trim(model%conserved_quantities(ivar)%name)//'_global_ave',                 &
         model_time, 'global mean '//trim(model%conserved_quantities(ivar)%long_name)//' in liquid seawater', &
         trim(model%conserved_quantities(ivar)%units),                     &
         missing_value = -1.0e+10)
   end do

   call fabm_set_domain(model,iec-isc+1,jec-jsc+1,nk)
   call fabm_set_mask(model,grid_tmask(isc:iec,jsc:jec,:),grid_tmask(isc:iec,jsc:jec,1))
   call model%set_surface_index(1)
   call model%set_bottom_index(grid_kmt(isc:iec,jsc:jec))

   write(stdout(),*)
   write(stdout(),*) trim(note_header), 'Tracer runs initialized'
   write(stdout(),*)

   allocate(clipped(isc:iec,jsc:jec,nk))

   allocate(work_dy(isc:iec,     size(model%state_variables)))
   allocate(work_dy_sf(isc:iec,  size(model%surface_state_variables)))
   allocate(work_dy_bt(isc:iec,  size(model%bottom_state_variables)))
   allocate(w      (isc:iec,nk+1,size(model%state_variables)))
   allocate(adv    (        nk+1,size(model%state_variables)))
   if (any(inds_cons_tot>0).or.any(inds_cons_ave>0)) then
      allocate(work_cons(isd:ied,size(model%conserved_quantities)))
      work_cons = 0
   end if

   ! Provide FABM with current biogeochemical state, as maintained by MOM.
   ! We create a copy of the biogeochemical state as maintained by MOM, in order to
   ! - convert tracer per seawater mass (MOM) to tracer per volume (FABM) by multiplying with density
   ! - allow us to repair the state (ensure that FABM sees valid values) while not affecting [clipping] the mass in the system.
   do ivar=1,size(model%state_variables)
      call model%link_interior_state_data(ivar,t_prog(interior_state_indices(ivar))%wrk1(isc:iec,jsc:jec,:))
   end do

   allocate(bottom_state(isc:iec,jsc:jec,size(model%bottom_state_variables)))
   bottom_state = 0.0
   call model%link_all_bottom_state_data(bottom_state)
   allocate(bottom_state_indices(size(model%bottom_state_variables)))

   allocate(surface_state(isc:iec,jsc:jec,size(model%surface_state_variables)))
   surface_state = 0.0
   call model%link_all_surface_state_data(surface_state)
   allocate(surface_state_indices(size(model%surface_state_variables)))

   ! Add surface/bottom attached tracers to 2D restart file.
   ! Also make them available to he output system by registering them as diagnostics.
   do ivar=1,size(model%bottom_state_variables)
      id_restart = register_restart_field(restart_2d, default_2d_restart_file, trim(model%bottom_state_variables(ivar)%name), bottom_state(:,:,ivar), domain=mpp_domain2d)
      bottom_state_indices(ivar) = register_diag_field('ocean_model',             &
         trim(model%bottom_state_variables(ivar)%name), grid_tracer_axes(1:2), &
         model_time, trim(model%bottom_state_variables(ivar)%long_name),       &
         trim(model%bottom_state_variables(ivar)%units),                       &
         missing_value = model%bottom_state_variables(ivar)%missing_value)
   end do
   do ivar=1,size(model%surface_state_variables)
      id_restart = register_restart_field(restart_2d, default_2d_restart_file, trim(model%surface_state_variables(ivar)%name), surface_state(:,:,ivar), domain=mpp_domain2d)
      surface_state_indices(ivar) = register_diag_field('ocean_model',             &
         trim(model%surface_state_variables(ivar)%name), grid_tracer_axes(1:2), &
         model_time, trim(model%surface_state_variables(ivar)%long_name),       &
         trim(model%surface_state_variables(ivar)%units),                       &
         missing_value = model%surface_state_variables(ivar)%missing_value)
   end do

   ! Load 2D state from restart file.
   call restore_state(restart_2d)

   if (fabm_variable_needs_values(model,id_par)) call model%link_interior_data(id_par, t_diag(index_irr)%field(isc:iec,jsc:jec,:))

   call model%link_horizontal_data(standard_variables%surface_downwelling_shortwave_flux,               swflx(isc:iec,jsc:jec))
   call model%link_horizontal_data(standard_variables%surface_downwelling_photosynthetic_radiative_flux,swflx_vis(isc:iec,jsc:jec))
   call model%link_horizontal_data(standard_variables%bottom_stress,                                    current_wave_stress(isc:iec,jsc:jec))
   call model%link_interior_data  (standard_variables%cell_thickness,                                   dzt(isc:iec,jsc:jec,:))
   call model%link_horizontal_data(standard_variables%longitude,                                        grid_xt(isc:iec,jsc:jec))
   call model%link_horizontal_data(standard_variables%latitude,                                         grid_yt(isc:iec,jsc:jec))
   call model%link_horizontal_data(standard_variables%bottom_depth_below_geoid,                         grid_ht(isc:iec,jsc:jec))
   call model%link_scalar         (standard_variables%number_of_days_since_start_of_the_year,           days_since_start_of_the_year)

end subroutine  ocean_fabm_start  !}
! </SUBROUTINE> NAME="ocean_fabm_start"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_source">
!
! <DESCRIPTION>
!     compute the source terms for the BIOTICs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!

subroutine ocean_fabm_source(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,     &
     T_prog, T_diag, taum1, model_time, grid_dat, grid_zw, grid_tmask, grid_kmt, rho_dzt, dzt, Dens, dtts, mpp_domain2d)  !{

   integer, intent(in)                                             :: isc
   integer, intent(in)                                             :: iec
   integer, intent(in)                                             :: jsc
   integer, intent(in)                                             :: jec
   integer, intent(in)                                             :: nk
   integer, intent(in)                                             :: isd
   integer, intent(in)                                             :: ied
   integer, intent(in)                                             :: jsd
   integer, intent(in)                                             :: jed
   type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
   type(ocean_diag_tracer_type), intent(inout), dimension(:)       :: T_diag
   integer, intent(in)                                             :: taum1
   type(time_type), intent(in)                                     :: model_time
   real, dimension(isd:,jsd:), intent(in)                          :: grid_dat
   real, dimension(nk), intent(in)                                 :: grid_zw
   real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask
   integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt
   real, dimension(isd:,jsd:,:,:), intent(in)                      :: rho_dzt
   real, dimension(isd:ied,jsd:jed,nk), intent(in)                 :: dzt
   type(ocean_density_type), intent(in)                            :: Dens
   real, intent(in)                                                :: dtts
   type(domain2d), intent(in)                                      :: mpp_domain2d

   character(len=64), parameter    :: sub_name = 'ocean_fabm_source'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   real, parameter :: SECONDS_PER_DAY  = 8.640000E+04

   integer :: i
   integer :: j
   integer :: k
   integer :: ivar
   logical :: used,valid
   real :: total_tracer(size(model%conserved_quantities)),total_volume
   integer :: years, months, days, hours, minutes, seconds
   type (type_forcing_2d), pointer :: forcing_2d
   type (type_forcing_3d), pointer :: forcing_3d

   call mpp_clock_begin(id_clock_fabm_setup_env)

   ! Set a reasonable chlorophyll concentration
   T_diag(index_chl)%field(isc:iec,jsc:jec,:) = 0.08

   ! Send environmental variables that MOM can provide.
   call model%link_interior_data(id_temp,t_prog(indtemp)%field(isc:iec,jsc:jec,:,taum1))
   call model%link_interior_data(id_salt,t_prog(indsal )%field(isc:iec,jsc:jec,:,taum1))
   call model%link_interior_data(id_dens,Dens%rho             (isc:iec,jsc:jec,:,taum1))
   call model%link_interior_data(id_pres,Dens%pressure_at_depth(isc:iec,jsc:jec,:))

   ! Compute number of days since start of the year
   call get_date(model_time, years, months, days, hours, minutes, seconds)
   call get_time(model_time-set_date(years,1,1), seconds, days)
   days_since_start_of_the_year = days + seconds/SECONDS_PER_DAY

   ! Convert tracer per seawater mass (MOM) to tracer per volume (FABM) by multiplying with density (kg m-3)
   do ivar=1,size(model%state_variables)
      t_prog(interior_state_indices(ivar))%wrk1(isc:iec,jsc:jec,:) = t_prog(interior_state_indices(ivar))%field(isc:iec,jsc:jec,:,taum1)*Dens%rho(isc:iec,jsc:jec,:,taum1)
   end do

   if (.not.initialization_complete) then
      ! This is the first time that this subroutine is called. Therefore, probe for
      ! each FABM variable whether external fields have been provided. If so, register
      ! the field name and increment the number of external fields. The total count
      ! is then used to allocate arrays for the external data, which are then coupled to
      ! FABM.
      !
      ! It would have been nice to do this during initialization (ocean_fabm_init or ocean_fabm_start),
      ! instead of on demand, but during initialization the data override module has not
      ! been told about the ocean domain yet (the relevant data_override_init call in
      ! coupler_main/ocean_solo occurs after ocean_model_init, which in turn call both init and start).
      ! Thus, data_override cannot be called from ocean_fabm_init or ocean_fabm_start.

      ! Process interior variables.
      allocate(forcing_3d)
      allocate(forcing_3d%data(isc:iec,jsc:jec,nk))
      do ivar=1,size(model%dependencies)
         forcing_3d%name = 'fabm_'//trim(model%dependencies(ivar))
         call data_override('OCN', forcing_3d%name, forcing_3d%data, model_time, override=used )
         if (used) then
            call model%link_interior_data(fabm_get_bulk_variable_id(model,model%dependencies(ivar)),forcing_3d%data,source=data_source_user)
            forcing_3d%next => first_forcing_3d
            first_forcing_3d => forcing_3d
            allocate(forcing_3d)
            allocate(forcing_3d%data(isc:iec,jsc:jec,nk))
         end if
      end do
      deallocate(forcing_3d)

      ! Process horizontal variables.
      allocate(forcing_2d)
      allocate(forcing_2d%data(isc:iec,jsc:jec))
      do ivar=1,size(model%dependencies_hz)
         forcing_2d%name = 'fabm_'//trim(model%dependencies_hz(ivar))
         call data_override('OCN', forcing_2d%name, forcing_2d%data, model_time, override=used )
         if (used) then
            call model%link_horizontal_data(fabm_get_horizontal_variable_id(model,model%dependencies_hz(ivar)),forcing_2d%data,source=data_source_user)
            forcing_2d%next => first_forcing_2d
            first_forcing_2d => forcing_2d
            allocate(forcing_2d)
            allocate(forcing_2d%data(isc:iec,jsc:jec))
         end if
      end do
      deallocate(forcing_2d)

      ! Make sure FABM has all dependencies fulfilled.
      call fabm_check_ready(model)

      initialization_complete = .true.
   end if

   ! Update values of external variables.
   forcing_3d => first_forcing_3d
   do while (associated(forcing_3d))
      call data_override('OCN',forcing_3d%name, forcing_3d%data, model_time)
      forcing_3d => forcing_3d%next
   end do
   forcing_2d => first_forcing_2d
   do while (associated(forcing_2d))
      call data_override('OCN',forcing_2d%name, forcing_2d%data, model_time)
      forcing_2d => forcing_2d%next
   end do

   call mpp_clock_end(id_clock_fabm_setup_env)

   call mpp_clock_begin(id_clock_fabm_repair)

   ! Repair biogeochemical state at the start of the time step.
   ! This modifies the data we sent to link_interior_state_data in-place!
   do j = jsc, jec  !{
      call fabm_check_surface_state(model,1,iec-isc+1,j-jsc+1,.true.,valid)
      call fabm_check_bottom_state(model,1,iec-isc+1,j-jsc+1,.true.,valid)
   end do
   do k = 1, nk  !{
      do j = jsc, jec  !{
         call fabm_check_state(model,1,iec-isc+1,j-jsc+1,k,.true.,valid)
      end do
   end do

   ! Send per-grid-point, per-variable clipping-induced change to diagnostic manager
   do ivar=1,size(model%state_variables)
      if (inds_clip(ivar)>0) then
         ! Compute value change due to clipping
         clipped(isc:iec,jsc:jec,:) = t_prog(interior_state_indices(ivar))%wrk1(isc:iec,jsc:jec,:) - t_prog(interior_state_indices(ivar))%field(isc:iec,jsc:jec,:,taum1)*Dens%rho(isc:iec,jsc:jec,:,taum1)
         used = send_data(inds_clip(ivar),clipped(isc:iec,jsc:jec,:),model_time,rmask=grid_tmask(isc:iec,jsc:jec,:))
      end if
   end do

   call mpp_clock_end(id_clock_fabm_repair)

   call mpp_clock_begin(id_clock_fabm_conservation)

   if (any(inds_cons_tot>0).or.any(inds_cons_ave>0)) then
      ! Get conserved quantities from FABM
      total_tracer = 0.0
      do k = 1, nk  !{
         do j = jsc, jec  !{
            call fabm_get_conserved_quantities(model,1,iec-isc+1,j-jsc+1,k,work_cons(isc:iec,:))
            do ivar=1,size(model%conserved_quantities)
               if (inds_cons_tot(ivar)>0.or.inds_cons_ave(ivar)>0) &
                  total_tracer(ivar) = total_tracer(ivar) + sum(grid_tmask(:,j,k)*grid_dat(:,j)*work_cons(:,ivar)*dzt(:,j,k))
            end do
         end do
      end do
      do j = jsc, jec  !{
         call fabm_get_horizontal_conserved_quantities(model,1,iec-isc+1,j-jsc+1,work_cons(isc:iec,:))
         do ivar=1,size(model%conserved_quantities)
            if (inds_cons_tot(ivar)>0.or.inds_cons_ave(ivar)>0) &
               total_tracer(ivar) = total_tracer(ivar) + sum(grid_tmask(:,j,1)*grid_dat(:,j)*work_cons(:,ivar))
         end do
      end do

      ! If we need the global average for any variable, calculate the global integral of mass here.
      ! (we use mass rather than volume because biogeochemical variables ar expresed in per mass rather than per volume)
      if (any(inds_cons_ave>0)) then
         total_volume = 0.0
         do k=1,nk
            total_volume = total_volume + sum(grid_tmask(:,:,k)*grid_dat(:,:)*dzt(:,:,k))
            !if(have_obc) tracer_k(:,:) = tracer_k(:,:)*Grd%obc_tmask(:,:)
            call mpp_sum(total_volume)
         end do
      end if

      do ivar=1,size(model%conserved_quantities)
         if (inds_cons_tot(ivar)>0.or.inds_cons_ave(ivar)>0) then
            ! Integrate tracer across full domain.
            !if(have_obc) tracer_k(:,:) = tracer_k(:,:)*Grd%obc_tmask(:,:)
            call mpp_sum(total_tracer(ivar))

            ! Send global integral and mean to diagnostic manager.
            if (inds_cons_tot(ivar)>0) used = send_data(inds_cons_tot(ivar),total_tracer(ivar)*1e-21,model_time)
            if (inds_cons_ave(ivar)>0) used = send_data(inds_cons_ave(ivar),total_tracer(ivar)/total_volume,model_time)
         end if
      end do
   end if

   call mpp_clock_end(id_clock_fabm_conservation)

   do k = 1, nk  !{
      do j = jsc, jec  !{
         call fabm_get_light_extinction(model,1,iec-isc+1,j-jsc+1,k,work_dy(isc:iec,1))
      end do
   end do

   do j = jsc, jec  !{
      do i = isc, iec  !{
         call fabm_get_light(model,1,nk,i-isc+1,j-jsc+1)
      end do
   end do

   call mpp_clock_begin(id_clock_fabm_bottom)

   bottom_fluxes = 0.0
   do j = jsc, jec  !{
      work_dy_bt = 0.0
      call fabm_do_bottom(model,1,iec-isc+1,j-jsc+1,bottom_fluxes(isc:iec,j,:),work_dy_bt(isc:iec,:))
      if (any(isnan(bottom_fluxes(isc:iec,j,:)))) then
         call mpp_error(FATAL,trim(error_header) // ' NaN in FABM bottom fluxes.')
      end if
      if (any(isnan(work_dy_bt(isc:iec,:)))) then
         call mpp_error(FATAL,trim(error_header) // ' NaN in FABM bottom sink/source terms.')
      end if
      bottom_state(isc:iec,j,:) = bottom_state(isc:iec,j,:) + dtts*work_dy_bt(isc:iec,:)
   end do  !} j

   call mpp_clock_end(id_clock_fabm_bottom)

   call mpp_clock_begin(id_clock_fabm_surface)

   surface_fluxes = 0.0
   do j = jsc, jec  !{
      work_dy_sf = 0.0
      call fabm_do_surface(model,1,iec-isc+1,j-jsc+1,surface_fluxes(isc:iec,j,:),work_dy_sf(isc:iec,:))
      if (any(isnan(surface_fluxes(isc:iec,j,:)))) then
         call mpp_error(FATAL,trim(error_header) // ' NaN in FABM surface fluxes.')
      end if
      if (any(isnan(work_dy_sf(isc:iec,:)))) then
         call mpp_error(FATAL,trim(error_header) // ' NaN in FABM surface sink/source terms.')
      end if
      surface_state(isc:iec,j,:) = surface_state(isc:iec,j,:) + dtts*work_dy_sf(isc:iec,:)
   end do  !} j

   call mpp_clock_end(id_clock_fabm_surface)

   call mpp_clock_begin(id_clock_fabm_interior)

   do k = 1, nk  !{
      do j = jsc, jec  !{

         ! Initialize derivatives to zero, because FABM will increment/decrement values rather than set them.
         work_dy = 0.0

         call fabm_do(model,1,iec-isc+1,j-jsc+1,k,work_dy(isc:iec,:))
         if (any(isnan(work_dy(isc:iec,:)))) then
            call mpp_error(FATAL,trim(error_header) // ' NaN in FABM sink/source terms.')
         end if

         if (.not.disable_sources) then
            ! Update tendencies with current sink and source terms.
            do ivar=1,size(model%state_variables)
               if (inplace_repair) then
                  ! We need to repair the actual state [clipping = non-conservative!], not only the strate as seen by the biogeochemical models.
                  ! Add clipping difference divided by time step as tracer source term
                  work_dy(isc:iec,ivar) = work_dy(isc:iec,ivar) + &
                                          (t_prog(interior_state_indices(ivar))%wrk1(isc:iec,j,k) - t_prog(interior_state_indices(ivar))%field(isc:iec,j,k,taum1)*Dens%rho(isc:iec,j,k,taum1))/dtts
               endif

               ! Increment tracer tendency
               ! This is the change per time in (tracer per seawater mass), multiplied by density and layer thickness.
               ! Note that (tracer per seawater mass)*density is tracer per volume, as FABM already uses internally.
               t_prog(interior_state_indices(ivar))%th_tendency(isc:iec,j,k) = t_prog(interior_state_indices(ivar))%th_tendency(isc:iec,j,k) &
                  + work_dy(isc:iec,ivar)*dzt(isc:iec,j,k)*grid_tmask(isc:iec,j,k)
            end do
         end if

      end do  !} j
   end do  !} k

   call mpp_clock_end(id_clock_fabm_interior)

   ! Save interior diagnostic variables.
   do ivar=1,size(model%diagnostic_variables)
      if (interior_diagnostic_indices(ivar) > 0) then
         used = send_data(interior_diagnostic_indices(ivar), fabm_get_interior_diagnostic_data(model,ivar), &
            model_time, rmask = grid_tmask(isc:iec,jsc:jec,:))
      end if
   end do

   ! Save horizontal diagnostic variables.
   do ivar=1,size(model%horizontal_diagnostic_variables)
      if (horizontal_diagnostic_indices(ivar) > 0) then
         used = send_data(horizontal_diagnostic_indices(ivar), fabm_get_horizontal_diagnostic_data(model,ivar), &
            model_time, rmask = grid_tmask(isc:iec,jsc:jec,1))
      end if
   end do

   ! Save 2D state (surface- and bottom-attached state variables)
   do ivar=1,size(model%bottom_state_variables)
      if (bottom_state_indices(ivar) > 0) &
         used = send_data(bottom_state_indices(ivar), bottom_state(:,:,ivar), model_time, rmask = grid_tmask(isc:iec,jsc:jec,1))
   end do
   do ivar=1,size(model%surface_state_variables)
      if (surface_state_indices(ivar) > 0) &
         used = send_data(surface_state_indices(ivar), surface_state(:,:,ivar), model_time, rmask = grid_tmask(isc:iec,jsc:jec,1))
   end do

   call mpp_clock_begin(id_clock_fabm_vert_mov)

   ! Vertical movement is applied with a first-order upwind scheme.
   do j = jsc, jec
      ! For every i: get sinking speed in cell centers, over entire column, for all state variables.
      do k=1,nk
         call fabm_get_vertical_movement(model,1,iec-isc+1,j-jsc+1,k,w(:,k,:))
      end do

      do i = isc, iec
         if (grid_tmask(i,j,1)/=1.) cycle

         ! Interpolate to sinking speed (m s-1) at interfaces
         do ivar=1,size(model%state_variables)
            w(i,2:grid_kmt(i,j),ivar) = (w(i,2:grid_kmt(i,j),ivar) + w(i,1:grid_kmt(i,j)-1,ivar))/2
         end do
         w(i,1,              :) = 0.0   ! Surface boundary condition
         w(i,grid_kmt(i,j)+1,:) = 0.0   ! Bottom boundary condition

         adv = 0.0

         do ivar=1,size(model%state_variables)

            ! Get upstream-biased tracer flux at all interfaces.
            do k=2,grid_kmt(i,j)
               if (w(i,k,ivar)>0.) then
                  ! floating
                  adv(k,ivar) = w(i,k,ivar)*t_prog(interior_state_indices(ivar))%field(i,j,k,taum1)*Dens%rho(i,j,k,taum1)
               elseif (w(i,k,ivar)<0.) then
                  ! sinking
                  adv(k,ivar) = w(i,k,ivar)*t_prog(interior_state_indices(ivar))%field(i,j,k-1,taum1)*Dens%rho(i,j,k-1,taum1)
               end if
            end do

            if (.not.disable_vertical_movement) then
               ! Add transport to tracer tendencies (conservative formulation)
               ! Note: normally the transport terms should be divided by the layer thickness.
               ! However, as MOM needs thickness (and density)-weighted tendencies (multiplication by thickness)
               ! that can be skipped here.
               do k=1,grid_kmt(i,j)
                  t_prog(interior_state_indices(ivar))%th_tendency(i,j,k) = t_prog(interior_state_indices(ivar))%th_tendency(i,j,k) + &
                     grid_tmask(i,j,k)*(adv(k+1,ivar)-adv(k,ivar)) !*Dens%rho(i,j,k,taum1)
               end do  !} k
            end if

         end do  !} ivar
      end do  !} i
   end do  !} j

   call mpp_clock_end(id_clock_fabm_vert_mov)

end subroutine  ocean_fabm_source  !}
! </SUBROUTINE> NAME="ocean_fabm_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM basecode
! </DESCRIPTION>
!

subroutine ocean_fabm_tracer  !{

   use mpp_mod, only : mpp_sum

   character(len=64), parameter    :: sub_name = 'ocean_fabm_tracer'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

end subroutine  ocean_fabm_tracer  !}
! </SUBROUTINE> NAME="ocean_fabm_tracer"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_init_sfc">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocean_fabm_start
! </DESCRIPTION>

subroutine ocean_fabm_init_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,   &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     Ocean_fields, T_prog, rho, taum1, model_time, grid_tmask)  !{

   integer, intent(in)                                     :: isc
   integer, intent(in)                                     :: iec
   integer, intent(in)                                     :: jsc
   integer, intent(in)                                     :: jec
   integer, intent(in)                                     :: nk
   integer, intent(in)                                     :: isd
   integer, intent(in)                                     :: ied
   integer, intent(in)                                     :: jsd
   integer, intent(in)                                     :: jed
   integer, intent(in)                                     :: isc_bnd
   integer, intent(in)                                     :: iec_bnd
   integer, intent(in)                                     :: jsc_bnd
   integer, intent(in)                                     :: jec_bnd
   type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
   type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
   real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
   integer, intent(in)                                     :: taum1
   type(time_type), intent(in)                             :: model_time
   real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

   character(len=64), parameter    :: sub_name = 'ocean_fabm_init_sfc'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

  allocate(surface_fluxes(isc:iec,jsc:jec,1:size(model%state_variables)))
  surface_fluxes = 0.0

  allocate(bottom_fluxes(isc:iec,jsc:jec,1:size(model%state_variables)))
  bottom_fluxes = 0.0

end subroutine ocean_fabm_init_sfc  !}
! </SUBROUTINE> NAME="ocean_fabm_init_sfc"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_fabm_avg_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd, Ocean_fields, Ocean_avg_kount, grid_tmask)  !{

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: iec_bnd
integer, intent(in)                                     :: jsc_bnd
integer, intent(in)                                     :: jec_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
integer                                                 :: Ocean_avg_kount
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

end subroutine ocean_fabm_avg_sfc  !}
! </SUBROUTINE> NAME="ocean_fabm_avg_sfc"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_sfc_end">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_fabm_sfc_end  !{

end subroutine ocean_fabm_sfc_end  !}
! </SUBROUTINE> NAME="ocean_fabm_sfc_end"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocean_fabm_bbc(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, grid_kmt)  !{

   integer, intent(in)                                             :: isc
   integer, intent(in)                                             :: iec
   integer, intent(in)                                             :: jsc
   integer, intent(in)                                             :: jec
   integer, intent(in)                                             :: isd
   integer, intent(in)                                             :: ied
   integer, intent(in)                                             :: jsd
   integer, intent(in)                                             :: jed
   type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
   integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt

   character(len=64), parameter    :: sub_name = 'ocean_fabm_bbc'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   integer :: ivar

   do ivar=1,size(model%state_variables)
      ! btf is defined as bottom tracer flux [rho*m/sec*tracer concen] in ocean_types_mod.
      ! Since MOM defines tracer concentration as tracer quantity per seawater mass [not volume!],
      ! multiplication with rho (kg m-3) ensures it becomes tracer quantity per volume, as FABM already uses internally.
      ! Note: btf<0 means tracer entering the ocean (see e.g. geothermal heating implementation in ocean_bbc_mod),
      ! but FABM has bottom fluxes > 0 when entering the ocean. Hence the minus sign.
      t_prog(interior_state_indices(ivar))%btf(isc:iec,jsc:jec) = -bottom_fluxes(isc:iec,jsc:jec,ivar)
   end do

end subroutine  ocean_fabm_bbc  !}
! </SUBROUTINE> NAME="ocean_fabm_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_end">
!
! <DESCRIPTION>
!     Clean up various BIOTIC quantities for this run.
! </DESCRIPTION>
!

subroutine ocean_fabm_end(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
     T_prog, grid_dat, grid_tmask, mpp_domain2d, rho_dzt, taup1)  !{

   integer, intent(in)                                     :: isc
   integer, intent(in)                                     :: iec
   integer, intent(in)                                     :: jsc
   integer, intent(in)                                     :: jec
   integer, intent(in)                                     :: nk
   integer, intent(in)                                     :: isd
   integer, intent(in)                                     :: ied
   integer, intent(in)                                     :: jsd
   integer, intent(in)                                     :: jed
   type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
   integer, intent(in)                                     :: taup1
   real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
   real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
   type(domain2d), intent(in)                              :: mpp_domain2d
   real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

   character(len=64), parameter    :: sub_name = 'ocean_fabm_end'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   ! Save additional restart file with 2d state variables
  call save_restart(restart_2d)
end subroutine  ocean_fabm_end  !}
! </SUBROUTINE> NAME="ocean_fabm_end"

subroutine mom_driver_fatal_error(self,location,message)
   class (type_mom_driver), intent(inout) :: self
   character(len=*),        intent(in)    :: location,message

   call mpp_error(FATAL, trim(location)//': '//trim(message))
end subroutine

subroutine mom_driver_log_message(self,message)
   class (type_mom_driver), intent(inout) :: self
   character(len=*),        intent(in)    :: message

   write (stdout(),*) trim(message)
end subroutine

end module  ocean_fabm_mod  !}
