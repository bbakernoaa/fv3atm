module CCPP_driver

  use ufs_ccpp_cap,       only: ccpp_register,               &
                                ccpp_init,                   &
                                ccpp_physics_init,           &
                                ccpp_physics_timestep_init,  &
                                ccpp_physics_run,            &
                                ccpp_physics_timestep_final, &
                                ccpp_physics_final,          &
                                ccpp_final
  use CCPP_data,          only: GFS_control,                 &
                                GFS_Intdiag,                 &
                                GFS_Interstitial
  use iso_fortran_env,    only: error_unit

#ifdef CATCHEM_CCPP
  ! CATChem CCPP process-suite integration: the host owns the CATChem context
  ! lifecycle (create/bind/initialize/destroy), advances the host clock and
  ! interval token each timestep, runs the catchem1/catchem2 groups after the
  ! GFS physics, and exchanges advected tracers with the framework-owned CCPP
  ! constituent array. Mirrors the CCPP-SCM host contract (CATChem
  ! tests/ccpp-scm patches 0002/0003).
  use iso_c_binding,      only: c_int, c_int64_t
  use machine,            only: kind_phys
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use CATChem_API,        only: CATChem_Manifest, CATChem_Model
  use CATChem_Context_Registry, only: ccpp_context_handle, ccpp_context_create, &
                                ccpp_context_bind_manifest, ccpp_context_initialize_model, &
                                ccpp_context_destroy, CCPP_CTX_OK
  use catchem_ccpp_context, only: catchem_ccpp_context_model
  use catchem_bridge_precision, only: fp
  use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t, ccpp_constituent_prop_ptr_t
  use ccpp_host_constituents, only: ccpp_register_constituents, ccpp_initialize_constituents, &
                                ccpp_deallocate_dynamic_constituents, ccpp_constituents_array, &
                                ccpp_model_const_properties
  use CCPP_data,          only: GFS_statein,                 &
                                GFS_stateout,                &
                                GFS_coupling,                &
                                GFS_sfcprop,                 &
                                catchem_context_id,          &
                                catchem_context_generation,  &
                                catchem_interval_id,         &
                                catchem_interval_token,      &
                                catchem_year,                &
                                catchem_month,               &
                                catchem_day,                 &
                                catchem_hour,                &
                                catchem_minute,              &
                                catchem_second,              &
                                catchem_day_of_year,         &
                                catchem_duration_seconds,    &
                                catchem_process_subcycle,    &
                                catchem_coupling_site
#endif

  implicit none

  !--------------------------------------------------------!
  ! CCPP control (mandatory) data.
  !--------------------------------------------------------!
  integer :: mythread
  integer :: nthreads
  integer :: nphys_threads
  integer :: lb
  integer :: ub
  integer :: errflg
  character(len=512) :: errmsg
  character(len=256) :: ccpp_suite='undefined'
  character(len=256) :: group_name='undefined'

!--------------------------------------------------------!
!  Number of OpenMP threads                             !
!--------------------------------------------------------!
  integer :: nthrds

#ifdef CATCHEM_CCPP
  ! CATChem CCPP host state: enabled only when the active suite is a CATChem
  ! process suite (suite name contains 'CATChem'). The context handle is
  ! created in the 'register' step (before ccpp_register, whose scheme
  ! register phases require a live context) and destroyed in 'physics_final'.
  logical :: catchem_enabled = .false.
  type(ccpp_context_handle) :: catchem_handle
#endif

  !--------------------------------------------------------!
  ! Public Entities
  !--------------------------------------------------------!
  public CCPP_step

CONTAINS
  !--------------------------------------------------------!
  !  CCPP step
  !--------------------------------------------------------!
  subroutine CCPP_step (step, nblks, ierr, dycore)

#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    character(len=*),         intent(in)  :: step
    integer,                  intent(in)  :: nblks
    integer,                  intent(out) :: ierr
    character(len=*),         intent(in)  :: dycore
    ! Local variables
    integer :: nb, nt
    integer :: kdt_iau
    logical :: iauwindow_center
    ! DH* 20210104 - remove kdt_rad when code to clear diagnostic buckets is removed
    integer :: kdt_rad

    ierr = 0

    ! CCPP register (same for all dynamical cores)
    if (trim(step)=="register") then
#ifdef CATCHEM_CCPP
       ! Create the CATChem context before ccpp_register: the suite's scheme
       ! register phases query the context species manifest.
       call catchem_host_context_init(ierr)
       if (ierr/=0) return
#endif
       call ccpp_register(ccpp_suite=trim(ccpp_suite), errmsg=errmsg, errflg=errflg)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_register: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if
#ifdef CATCHEM_CCPP
       ! Merge the suite-registered CATChem constituents into the framework
       ! registry and allocate the constituent array (after ccpp_register).
       call catchem_host_constituent_init(ierr)
       if (ierr/=0) return
#endif
       
    ! CCPP Framework init (same for all dynamical cores)
    else if (trim(step)=="init") then

       ! Get and set number of OpenMP threads (module
       ! variable) that are available to run physics
#ifdef _OPENMP
       nthrds = omp_get_max_threads()
#else
       nthrds = 1
#endif
      
       call ccpp_init(ccpp_suite=trim(ccpp_suite), errmsg=errmsg, errflg=errflg)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_init: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if

    ! Physics init (same for all dynamical cores)
    else if (trim(step)=="physics_init") then

       call ccpp_physics_init( ccpp_suite=trim(ccpp_suite), group_name='all', &
            errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,         &
            mythread=1, nthreads=1, nphys_threads=nthrds)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_init: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if

    ! Timestep init = time_vary (dycore specific)
    else if (trim(step)=="timestep_init") then

#ifdef CATCHEM_CCPP
       ! Advance the host clock/interval contract before the CATChem
       ! timestep_init scheme opens the shared physics interval.
       call catchem_host_timestep_init(ierr)
       if (ierr/=0) return
#endif

       call ccpp_physics_timestep_init(ccpp_suite=trim(ccpp_suite), group_name='all', &
            errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,         &
            mythread=1, nthreads=1, nphys_threads=nthrds)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_timestep_init for group time_vary: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! DH* 20210104 - this block of code will be removed once the CCPP framework    !
       ! fully supports handling diagnostics through its metadata, work in progress   !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !--- determine if radiation diagnostics buckets need to be cleared
       if (nint(GFS_control%fhzero*3600) >= nint(max(GFS_control%fhswr,GFS_control%fhlwr))) then
          if (mod(GFS_control%kdt,GFS_control%nszero) == 1) then
             call GFS_Intdiag%rad_zero(GFS_control)
          endif
       else
          kdt_rad = nint(min(GFS_control%fhswr,GFS_control%fhlwr)/GFS_control%dtp)
          if (mod(GFS_control%kdt,kdt_rad) == 1) then
             call GFS_Intdiag%rad_zero(GFS_control)
          endif
       endif

       !--- determine if physics diagnostics buckets need to be cleared
       iauwindow_center = .false.
       if (GFS_control%iau_offset > 0) then
          kdt_iau = nint(GFS_control%iau_offset*3600./GFS_control%dtp)
          if (GFS_control%kdt-1 == kdt_iau) then
             iauwindow_center = .true.
             if( GFS_control%me == 0)print *,'in ccpp step vary, iauwindow_center=',iauwindow_center,&
                  'kdt=',GFS_control%kdt,'dtp=',GFS_control%dtp,'iau_offset=',GFS_control%iau_offset
          endif
      endif
      if ((mod(GFS_control%kdt-1,GFS_control%nszero)) == 0) then
         call GFS_Intdiag%phys_zero(GFS_control, iauwindow_center=iauwindow_center)
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! *DH 20210104                                                                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Radiation, physics and and stochastic physics - threaded regions using blocked data structures
    else if (trim(step)=="radiation" .or. trim(step)=="physics" .or. trim(step)=="stochastics" .or. trim(step)=="microphysics") then

!$OMP parallel num_threads (nthrds)                        &
!$OMP          default (none)                              &
!$OMP          shared (nblks, ccpp_suite,                  &
!$OMP                  step, GFS_Control, GFS_Interstitial,&
!$OMP                  dycore, nthrds)                     &
!$OMP          private (nb, nt, errmsg, errflg)            &
!$OMP          reduction (+:ierr)
#ifdef _OPENMP
       nt = omp_get_thread_num()+1
#else
       nt = 1
#endif
!$OMP do schedule (dynamic,1)
       do nb = 1,nblks
          ! Allocate physics interstitals for current thread
          call GFS_Interstitial(nt)%create(ixs=GFS_control%chunk_begin(nb), ixe=GFS_control%chunk_end(nb), model=GFS_control)
          !--- Call CCPP radiation/physics/stochastics group
          if (trim(step)=="physics") then
             if (trim(dycore)=="fv3") then
                ! Reset GFS_Interstitial DDT fields for this thread
                call GFS_Interstitial(nt)%reset(GFS_control)
                ! Process-split physics
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="phys_ps", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group phys_ps: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                end if
                ! Time-split physics
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="phys_ts", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group phys_ts: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                end if
             endif
             if (trim(dycore)=="mpas") then
                ! Physics
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="physics", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group physics: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                endif
             endif
          else
             ! DH* WHY WAS THIS FOR RADIATION ONLY??? PERFORMANCE?
             ! Reset GFS_Interstitial DDT fields for this thread
             call GFS_Interstitial(nt)%reset(GFS_control)
             ! *DH
             ! Radiation
             if (trim(step)=="radiation") then
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="radiation", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group radiation: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                end if
             ! Microphysics (MPAS only)
             else if (trim(step)=="microphysics") then
                if (trim(dycore)=="mpas") then
                   call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="microphysics", &
                        errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                        mythread=nt, nthreads=nthrds, nphys_threads=1)
                   if (errflg/=0) then
                      write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group microphysics: ' // trim(errmsg) // '. Exiting...'
                      ierr = ierr + errflg
                   end if
                else
                   write(error_unit,'(a)') "An error occurred in ccpp_physics_run for group microphysics. Group microphysics only valid with MPAS dycore."
                   ierr = ierr + errflg
                endif
             ! Stochastic physics
             else if (trim(step)=="stochastics") then
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="stochastics", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group stochastics: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                end if
             ! Catchall
             else
                write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run: unknown group ' // trim(step) // '. Exiting...'
                ierr = ierr + 1
             end if
          endif
          call GFS_Interstitial(nt)%destroy(GFS_control)
       end do
!$OMP end do

!$OMP end parallel
       if (ierr/=0) return

#ifdef CATCHEM_CCPP
       ! CATChem process groups run after the GFS physics over all columns:
       ! the CATChem state contract requires every met binding to span the
       ! full column count, so the groups cannot run per block. A temporary
       ! full-domain interstitial supplies the interstitial-derived inputs.
       if (trim(step)=="physics" .and. trim(dycore)=="fv3") then
          call catchem_host_run(ierr)
          if (ierr/=0) return
       end if
#endif

    ! Timestep final = time_vary (same for all dynamical cores)
    else if (trim(step)=="timestep_final") then

       call ccpp_physics_timestep_final(ccpp_suite=trim(ccpp_suite), group_name="all", &
            errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,         &
            mythread=1, nthreads=1, nphys_threads=nthrds)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_timestep_final group time_vary: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if

    ! Physics final (same for all dynamical cores)
    else if (trim(step)=="physics_final") then
       call ccpp_physics_final(ccpp_suite=trim(ccpp_suite), group_name='all', &
            errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,         &
            mythread=1, nthreads=1, nphys_threads=nthrds)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_final: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if
#ifdef CATCHEM_CCPP
       ! Release the framework constituent array and destroy the CATChem
       ! context after the scheme final phases have completed.
       call catchem_host_final(ierr)
       if (ierr/=0) return
#endif

    ! Frameowrk final (same for all dynamical cores)
    else if (trim(step)=="final") then
       call ccpp_final(ccpp_suite=trim(ccpp_suite), errmsg=errmsg, errflg=errflg)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_final: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if
    ! Undefined CCPP step
    else
       write(error_unit,'(2a)') 'Error, undefined CCPP step ', trim(step)
       ierr=errflg
       return
    end if

  end subroutine CCPP_step

#ifdef CATCHEM_CCPP

  !> \brief Create and initialize the CATChem context for a CATChem suite.
  !! Reads the species manifest from the runtime YAML (CATCHEM_CONFIG env var
  !! or ./CATChem_new_config.yml), creates the registry context, binds the
  !! manifest and initializes the CATChem model. No-op for non-CATChem suites.
  subroutine catchem_host_context_init(ierr)
    integer, intent(out) :: ierr

    type(CATChem_Manifest) :: manifest
    character(len=512) :: config_path
    integer :: manifest_rc, species_count, env_status
    integer(c_int) :: status

    ierr = 0
    catchem_enabled = .false.
    if (index(trim(ccpp_suite), 'CATChem') == 0) return

    config_path = ''
    call get_environment_variable('CATCHEM_CONFIG', config_path, status=env_status)
    if (env_status /= 0 .or. len_trim(config_path) == 0) config_path = 'CATChem_new_config.yml'

    call manifest%open(trim(config_path), manifest_rc)
    if (manifest_rc /= 0) then
       write(error_unit,'(a)') 'CATChem: unable to read manifest '//trim(config_path)//': '//trim(manifest%last_error)
       ierr = 1
       return
    end if
    call manifest%species_count(species_count, manifest_rc)
    if (manifest_rc /= 0 .or. species_count <= 0) then
       write(error_unit,'(a)') 'CATChem: configuration has no valid species manifest: '//trim(config_path)
       ierr = 1
       return
    end if

    call ccpp_context_create(GFS_control%ncols, GFS_control%levs, species_count, catchem_handle, status)
    if (status == CCPP_CTX_OK) then
       call ccpp_context_bind_manifest(catchem_handle, trim(config_path), status)
    end if
    if (status == CCPP_CTX_OK) then
       call ccpp_context_initialize_model(catchem_handle, GFS_control%ncols, GFS_control%levs, status)
    end if
    if (status /= CCPP_CTX_OK) then
       write(error_unit,'(a,i0)') 'CATChem: unable to initialize the CCPP context, status=', status
       ierr = 1
       return
    end if
    catchem_context_id = catchem_handle%id
    catchem_context_generation = catchem_handle%generation
    call manifest%close(manifest_rc)
    catchem_enabled = .true.
    write(error_unit,'(a,a,a,i0,a,i0,a)') 'CATChem: CCPP context initialized for suite ', trim(ccpp_suite), &
         ' (', species_count, ' species, context ', catchem_context_id, ')'
  end subroutine catchem_host_context_init

  !> \brief Merge suite-registered constituents and allocate the array.
  !! Must run after ccpp_register (the scheme register phases have populated
  !! the suite constituent buffers) and before ccpp_physics_init.
  subroutine catchem_host_constituent_init(ierr)
    integer, intent(out) :: ierr

    type(ccpp_constituent_properties_t), allocatable :: host_consts(:)
    integer :: const_rc
    character(len=512) :: const_msg

    ierr = 0
    if (.not. catchem_enabled) return
    allocate(host_consts(0))
    call ccpp_register_constituents(host_consts, const_rc, const_msg)
    if (const_rc /= 0) then
       write(error_unit,'(a)') 'CATChem: constituent registration failed: '//trim(const_msg)
       ierr = max(1, const_rc)
       return
    end if
    call ccpp_initialize_constituents(GFS_control%ncols, GFS_control%levs, const_rc, const_msg)
    if (const_rc /= 0) then
       write(error_unit,'(a)') 'CATChem: constituent allocation failed: '//trim(const_msg)
       ierr = max(1, const_rc)
       return
    end if
    deallocate(host_consts)
  end subroutine catchem_host_constituent_init

  !> \brief Advance the host clock and interval token for this timestep.
  subroutine catchem_host_timestep_init(ierr)
    integer, intent(out) :: ierr

    ierr = 0
    if (.not. catchem_enabled) return
    if (catchem_interval_id == huge(catchem_interval_id)) then
       write(error_unit,'(a)') 'CATChem: interval identifier exhausted'
       ierr = 1
       return
    end if
    catchem_interval_id = catchem_interval_id + 1_c_int64_t
    catchem_interval_token = catchem_interval_id
    catchem_process_subcycle = 1_c_int
    catchem_coupling_site = 'catchem-processes'
    catchem_year   = GFS_control%jdat(1)
    catchem_month  = GFS_control%jdat(2)
    catchem_day    = GFS_control%jdat(3)
    catchem_hour   = GFS_control%jdat(5)
    catchem_minute = GFS_control%jdat(6)
    catchem_second = GFS_control%jdat(7)
    catchem_day_of_year = catchem_ordinal_day(catchem_year, catchem_month, catchem_day)
    catchem_duration_seconds = real(GFS_control%dtp, kind_phys)
  end subroutine catchem_host_timestep_init

  !> \brief Run the CATChem process groups after the GFS physics.
  !! Refreshes per-call host inputs, recreates a full-domain interstitial
  !! (the blocked physics loop destroys the per-thread interstitials, but the
  !! dry-deposition driver reads snow cover from GFS_Interstitial(1)), imports
  !! host tracers into the CCPP constituent array, runs catchem1
  !! (emissions/removal) and catchem2 (chemistry) over all columns, exports
  !! the updated constituents back to the host tracer state, and releases the
  !! temporary interstitial.
  subroutine catchem_host_run(ierr)
    integer, intent(out) :: ierr

    real(kind_phys) :: aero_colmass_before(GFS_control%ncols)

    ierr = 0
    if (.not. catchem_enabled) return

    ! Column-integrated aerosol mass before the CATChem groups for the net
    ! column-change diagnostic.
    aero_colmass_before = catchem_aerosol_colmass()

    ! The wet-deposition driver expects per-call rain amounts (m) in the
    ! coupling rain fields, with convective rain no greater than total rain.
    ! The GFS host otherwise accumulates these over the forecast, and
    ! rain_cpl tracks liquid-only rain while rainc_cpl includes convective
    ! snow, so refresh them from the per-timestep diagnostics (matching the
    ! NUOPC inst_rainfall_amount semantics).
    GFS_coupling%rain_cpl(:) = max(0.0_kind_phys, GFS_Intdiag%rain(:))
    GFS_coupling%rainc_cpl(:) = min(max(0.0_kind_phys, GFS_Intdiag%rainc(:)), GFS_coupling%rain_cpl(:))

    ! NoahMP leaves fill values (9.99e20) in xlaixy over non-land columns;
    ! the CATChem import contract rejects fill-valued host fields. Replace
    ! them with zero (no vegetation), which is physically correct over
    ! ocean/ice and harmless for the land scheme (ocean columns are skipped).
    where (.not. ieee_is_finite(GFS_sfcprop%xlaixy) .or. &
           abs(GFS_sfcprop%xlaixy) >= 1.0e20_kind_phys)
       GFS_sfcprop%xlaixy = 0.0_kind_phys
    end where

    ! Recreate the thread-1 interstitial for the full domain. The only
    ! interstitial field consumed by the CATChem drivers is the snow cover
    ! fraction (dry deposition); restore it from the persistent surface
    ! snow-cover field updated by the land scheme this timestep.
    call GFS_Interstitial(1)%create(ixs=1, ixe=GFS_control%ncols, model=GFS_control)
    call GFS_Interstitial(1)%reset(GFS_control)
    GFS_Interstitial(1)%snowc(:) = GFS_sfcprop%sncovr(:)

    call catchem_host_sync(.true., ierr)
    if (ierr /= 0) then
       call GFS_Interstitial(1)%destroy(GFS_control)
       return
    end if

    call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="catchem1", &
         errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,            &
         mythread=1, nthreads=1, nphys_threads=nthrds)
    if (errflg /= 0) then
       write(error_unit,'(a)') 'An error occurred in ccpp_physics_run for group catchem1: '//trim(errmsg)//'. Exiting...'
       ierr = errflg
       call GFS_Interstitial(1)%destroy(GFS_control)
       return
    end if

    call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="catchem2", &
         errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,            &
         mythread=1, nthreads=1, nphys_threads=nthrds)
    if (errflg /= 0) then
       write(error_unit,'(a)') 'An error occurred in ccpp_physics_run for group catchem2: '//trim(errmsg)//'. Exiting...'
       ierr = errflg
       call GFS_Interstitial(1)%destroy(GFS_control)
       return
    end if

    call catchem_host_sync(.false., ierr)
    call catchem_host_diagnostics(aero_colmass_before)
    call GFS_Interstitial(1)%destroy(GFS_control)
  end subroutine catchem_host_run

  !> \brief Names of the aerosol species carried as host tracers (lower case).
  !! Used for the aggregate wet-removal and net column-change diagnostics.
  pure logical function catchem_is_aerosol_tracer(name)
    character(len=*), intent(in) :: name
    character(len=8), parameter :: aerosols(15) = &
         ['seas1   ', 'seas2   ', 'seas3   ', 'seas4   ', 'seas5   ', &
          'dust1   ', 'dust2   ', 'dust3   ', 'dust4   ', 'dust5   ', &
          'bc1     ', 'bc2     ', 'oc1     ', 'oc2     ', 'so4     ']
    integer :: i
    catchem_is_aerosol_tracer = .false.
    do i = 1, size(aerosols)
       if (catchem_same_species(trim(aerosols(i)), name)) then
          catchem_is_aerosol_tracer = .true.
          return
       end if
    end do
  end function catchem_is_aerosol_tracer

  !> \brief Column-integrated aerosol mass (kg m-2) of the post-physics tracer
  !! state, summed over all aerosol bins and levels using layer pressure
  !! thickness (mass = sum(q * dp) / g).
  function catchem_aerosol_colmass() result(colmass)
    real(kind_phys) :: colmass(GFS_control%ncols)
    real(kind_phys), parameter :: grav = 9.80665_kind_phys
    integer :: itrc, k
    colmass = 0.0_kind_phys
    do itrc = 1, size(GFS_control%tracer_names)
       if (.not. catchem_is_aerosol_tracer(GFS_control%tracer_names(itrc))) cycle
       do k = 1, GFS_control%levs
          colmass(:) = colmass(:) + GFS_stateout%gq0(:,k,itrc) * &
               abs(GFS_statein%prsi(:,k+1) - GFS_statein%prsi(:,k)) / grav
       end do
    end do
  end function catchem_aerosol_colmass

  !> \brief Populate the CATChem host output diagnostics (GFS_Intdiag fields
  !! registered in GFS_diagnostics and selected through diag_table):
  !!  - catchem_seasemis/seasnum: sea-salt column emission fluxes from the
  !!    CATChem diagnostic manager;
  !!  - catchem_wetdep: total aerosol column wet removal (sum of the
  !!    per-species per-level removed mass divided by the physics timestep);
  !!  - catchem_aerodq: net aerosol column mass change rate across all CATChem
  !!    processes, from the before/after host tracer state.
  subroutine catchem_host_diagnostics(aero_colmass_before)
    real(kind_phys), intent(in) :: aero_colmass_before(:)

    type(CATChem_Model), pointer :: model
    real(fp), allocatable :: ddata(:,:,:)
    real(kind_phys) :: wetcol(GFS_control%ncols)
    integer(c_int) :: status
    integer :: rc, i
    character(len=8), parameter :: aerosols(15) = &
         ['seas1   ', 'seas2   ', 'seas3   ', 'seas4   ', 'seas5   ', &
          'dust1   ', 'dust2   ', 'dust3   ', 'dust4   ', 'dust5   ', &
          'bc1     ', 'bc2     ', 'oc1     ', 'oc2     ', 'so4     ']

    GFS_Intdiag%catchem_aerodq(:) = &
         (catchem_aerosol_colmass() - aero_colmass_before(:)) / GFS_control%dtp

    call catchem_ccpp_context_model(catchem_handle, model, status)
    if (status /= CCPP_CTX_OK .or. .not. associated(model)) return

    call model%get_diagnostic('seasalt_mass_emission_total', ddata, rc)
    if (rc == 0 .and. allocated(ddata)) then
       if (size(ddata,1)*size(ddata,2) == GFS_control%ncols) then
          GFS_Intdiag%catchem_seasemis(:) = reshape(sum(sum(ddata, dim=3), dim=2), [GFS_control%ncols])
       end if
       deallocate(ddata)
    end if
    call model%get_diagnostic('seasalt_number_emission_total', ddata, rc)
    if (rc == 0 .and. allocated(ddata)) then
       if (size(ddata,1)*size(ddata,2) == GFS_control%ncols) then
          GFS_Intdiag%catchem_seasnum(:) = reshape(sum(sum(ddata, dim=3), dim=2), [GFS_control%ncols])
       end if
       deallocate(ddata)
    end if

    wetcol = 0.0_kind_phys
    do i = 1, size(aerosols)
       call model%get_diagnostic('wetdep_mass_'//trim(aerosols(i)), ddata, rc)
       if (rc == 0 .and. allocated(ddata)) then
          if (size(ddata,1)*size(ddata,2) == GFS_control%ncols) then
             wetcol(:) = wetcol(:) + reshape(sum(sum(ddata, dim=3), dim=2), [GFS_control%ncols])
          end if
          deallocate(ddata)
       end if
    end do
    GFS_Intdiag%catchem_wetdep(:) = wetcol(:) / GFS_control%dtp
  end subroutine catchem_host_diagnostics

  !> \brief Exchange advected tracers with the CCPP constituent array.
  !! Species are matched by constituent diagnostic name against the host
  !! tracer names (case-insensitive; O3/o3mr/ozone_mixing_ratio are aliases).
  !! Constituents without a matching host tracer keep their framework values.
  subroutine catchem_host_sync(to_constituents, ierr)
    logical, intent(in)  :: to_constituents
    integer, intent(out) :: ierr

    type(ccpp_constituent_prop_ptr_t), pointer :: properties(:)
    real(kind_phys), pointer :: constituents(:,:,:)
    character(len=64) :: species_name, prop_msg
    integer :: prop_rc, const_idx, itrc, iprop

    ierr = 0
    constituents => ccpp_constituents_array()
    properties => ccpp_model_const_properties()
    if (.not. associated(constituents) .or. .not. associated(properties)) then
       write(error_unit,'(a)') 'CATChem: CCPP constituent storage or property registry is not initialized.'
       ierr = 1
       return
    end if
    if (size(GFS_stateout%gq0,1) /= size(constituents,1) .or. &
        size(GFS_stateout%gq0,2) /= size(constituents,2)) then
       write(error_unit,'(a)') 'CATChem: host tracer and CCPP constituent array extents do not agree.'
       ierr = 1
       return
    end if

    do iprop = 1, size(properties)
       call properties(iprop)%diagnostic_name(species_name, prop_rc, prop_msg)
       if (prop_rc /= 0) then
          write(error_unit,'(a)') 'CATChem: unable to read a constituent diagnostic name: '//trim(prop_msg)
          ierr = max(1, prop_rc)
          return
       end if
       call properties(iprop)%const_index(const_idx, prop_rc, prop_msg)
       if (prop_rc /= 0 .or. const_idx < 1 .or. const_idx > size(constituents,3)) then
          write(error_unit,'(a)') 'CATChem: constituent has an invalid framework index: '//trim(species_name)
          ierr = max(1, prop_rc)
          return
       end if
       do itrc = 1, size(GFS_control%tracer_names)
          if (.not. catchem_same_species(species_name, GFS_control%tracer_names(itrc))) cycle
          ! Exchange with the POST-physics tracer state (Stateout%gq0): that is
          ! what the FV3 dycore consumes after physics; Statein%qgrs holds the
          ! pre-physics state and is refreshed from the dycore next step.
          if (to_constituents) then
             ! FV3 tracer transport can leave small negative mass ratios in the
             ! post-physics state before the dycore filler; CATChem's conversion
             ! contract requires nonnegative inputs. Clip on import and write the
             ! cleaned value back to the host state so any constituent gather
             ! sourced from gq0 also sees nonnegative values.
             constituents(:,:,const_idx) = max(0.0_kind_phys, GFS_stateout%gq0(:,:,itrc))
             GFS_stateout%gq0(:,:,itrc) = constituents(:,:,const_idx)
          else
             GFS_stateout%gq0(:,:,itrc) = constituents(:,:,const_idx)
             ! Guard the host state: report the first tracer that leaves the
             ! physical range so a runaway process is attributable.
             if (any(.not. ieee_is_finite(GFS_stateout%gq0(:,:,itrc))) .or. &
                 maxval(abs(GFS_stateout%gq0(:,:,itrc))) > 1.0e-2_kind_phys) then
                write(error_unit,'(a,a,a,es14.6,a,es14.6,a,i0)') 'CATChem: tracer ', &
                     trim(GFS_control%tracer_names(itrc)), ' out of range after export: min=', &
                     minval(GFS_stateout%gq0(:,:,itrc)), ' max=', maxval(GFS_stateout%gq0(:,:,itrc)), &
                     ' kdt=', GFS_control%kdt
             end if
          end if
          exit
       end do
    end do
  end subroutine catchem_host_sync

  !> \brief Release constituent storage and destroy the CATChem context.
  subroutine catchem_host_final(ierr)
    integer, intent(out) :: ierr

    integer(c_int) :: status

    ierr = 0
    if (.not. catchem_enabled) return
    call ccpp_deallocate_dynamic_constituents()
    call ccpp_context_destroy(catchem_handle, status)
    if (status /= CCPP_CTX_OK) then
       write(error_unit,'(a,i0)') 'CATChem: unable to destroy the CCPP context, status=', status
       ierr = 1
    end if
    catchem_context_id = 0_c_int
    catchem_context_generation = 0_c_int64_t
    catchem_enabled = .false.
  end subroutine catchem_host_final

  !> \brief Gregorian ordinal day for the host clock.
  pure integer function catchem_ordinal_day(year, month, day)
    integer, intent(in) :: year, month, day
    integer :: month_index
    integer :: month_lengths(12)
    month_lengths = [31,28,31,30,31,30,31,31,30,31,30,31]
    if (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) month_lengths(2) = 29
    catchem_ordinal_day = day
    do month_index = 1, month - 1
       if (month_index >= 1 .and. month_index <= 12) catchem_ordinal_day = catchem_ordinal_day + month_lengths(month_index)
    end do
  end function catchem_ordinal_day

  !> \brief Case-insensitive species name match with ozone aliases.
  pure logical function catchem_same_species(left, right)
    character(len=*), intent(in) :: left, right
    integer :: i, n, lc, rc
    if (catchem_is_ozone(left) .and. catchem_is_ozone(right)) then
       catchem_same_species = .true.
       return
    end if
    catchem_same_species = .false.
    n = len_trim(left)
    if (n /= len_trim(right)) return
    do i = 1, n
       lc = iachar(left(i:i))
       rc = iachar(right(i:i))
       if (lc >= iachar('a') .and. lc <= iachar('z')) lc = lc - 32
       if (rc >= iachar('a') .and. rc <= iachar('z')) rc = rc - 32
       if (lc /= rc) return
    end do
    catchem_same_species = .true.
  end function catchem_same_species

  !> \brief Recognize CATChem O3 and GFS ozone tracer identifiers as aliases.
  pure logical function catchem_is_ozone(name)
    character(len=*), intent(in) :: name
    character(len=32) :: key
    integer :: i, code
    key = adjustl(name)
    do i = 1, len_trim(key)
       code = iachar(key(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) key(i:i) = achar(code+32)
    end do
    catchem_is_ozone = trim(key)=='o3' .or. trim(key)=='o3mr' .or. trim(key)=='ozone_mixing_ratio'
  end function catchem_is_ozone

#endif

end module CCPP_driver
