module CCPP_data

!! \section arg_table_CCPP_data Argument Table
!! \htmlinclude CCPP_data.html
!!

    use CCPP_typedefs, only: GFS_interstitial_type,  &
                             GFDL_interstitial_type
    use GFS_typedefs,  only: GFS_control_type,       &
                             GFS_statein_type,       &
                             GFS_stateout_type,      &
                             GFS_grid_type,          &
                             GFS_tbd_type,           &
                             GFS_cldprop_type,       &
                             GFS_sfcprop_type,       &
                             GFS_radtend_type,       &
                             GFS_coupling_type,      &
                             GFS_diag_type
    use machine,       only: kind_phys
    use iso_c_binding, only: c_int, c_int64_t

    implicit none

    private

    public GFDL_interstitial,      &
           GFS_control,            &
           GFS_statein,            &
           GFS_stateout,           &
           GFS_grid,               &
           GFS_tbd,                &
           GFS_cldprop,            &
           GFS_sfcprop,            &
           GFS_radtend,            &
           GFS_coupling,           &
           GFS_intdiag,            &
           GFS_interstitial

    public catchem_context_id,             &
           catchem_context_generation,     &
           catchem_interval_id,            &
           catchem_interval_token,         &
           catchem_year,                   &
           catchem_month,                  &
           catchem_day,                    &
           catchem_hour,                   &
           catchem_minute,                 &
           catchem_second,                 &
           catchem_day_of_year,            &
           catchem_duration_seconds,       &
           catchem_process_subcycle,       &
           catchem_coupling_site

    !-------------------------------------------------------!
    !  GFS data containers;                                 !
    !  GFS_Interstitial has dimension nthreads              !
    !-------------------------------------------------------!
    type(GFS_control_type),                                    save, target :: GFS_control
    type(GFS_statein_type),                                    save, target :: GFS_statein
    type(GFS_stateout_type),                                   save, target :: GFS_stateout
    type(GFS_grid_type),                                       save, target :: GFS_grid
    type(GFS_tbd_type),                                        save, target :: GFS_tbd
    type(GFS_cldprop_type),                                    save, target :: GFS_cldprop
    type(GFS_sfcprop_type),                                    save, target :: GFS_sfcprop
    type(GFS_radtend_type),                                    save, target :: GFS_radtend
    type(GFS_coupling_type),                                   save, target :: GFS_coupling
    type(GFS_diag_type),                                       save, target :: GFS_intdiag
    type(GFS_interstitial_type),  dimension(:),   allocatable, save, target :: GFS_interstitial

    !------------------------------------------------------!
    !  CCPP data containers for dynamics (fast physics)    !
    !------------------------------------------------------!
    type(GFDL_interstitial_type),                              save, target :: GFDL_interstitial

    !------------------------------------------------------!
    !  CATChem CCPP host contract state.                   !
    !  The context/interval identities are established by  !
    !  the host (CCPP_driver) before ccpp_register; the    !
    !  clock fields are refreshed each timestep.           !
    !------------------------------------------------------!
    integer(c_int),                        save :: catchem_context_id = 0_c_int
    integer(c_int64_t),                    save :: catchem_context_generation = 0_c_int64_t
    integer(c_int64_t),                    save :: catchem_interval_id = 0_c_int64_t
    integer(c_int64_t),                    save :: catchem_interval_token = 0_c_int64_t
    integer,                               save :: catchem_year = 2000
    integer,                               save :: catchem_month = 1
    integer,                               save :: catchem_day = 1
    integer,                               save :: catchem_hour = 0
    integer,                               save :: catchem_minute = 0
    integer,                               save :: catchem_second = 0
    integer,                               save :: catchem_day_of_year = 1
    real(kind_phys),                       save :: catchem_duration_seconds = 3600.0_kind_phys
    integer(c_int),                        save :: catchem_process_subcycle = 1_c_int
    character(len=256),                    save :: catchem_coupling_site = ''

end module CCPP_data
