program main
USE MOD_MESH
USE MOD_PARTIT
USE MOD_TRACER
USE MOD_DYN
USE MOD_ICE
USE MOD_PARSUP
USE par_support_interfaces
USE restart_derivedtype_module
USE fortran_utils
USE mod_dderived_type
implicit none

character(len=500)                            :: arg
character(len=500)                            :: dtype
character(len=500)                            :: dpath1
character(len=500)                            :: dpath2
character(len=500)                            :: tol
character(LEN=500)                            :: npepath1
character(LEN=500)                            :: npepath2
logical                                       :: dir_exist
type(t_mesh),       target, save              :: mesh
type(t_partit),     target, save              :: partit
type(t_tracer),     target, save              :: tracers1
type(t_tracer),     target, save              :: tracers2
type(t_dyn),        target, save              :: dyn1
type(t_dyn),        target, save              :: dyn2
type(t_ice),        target, save              :: ice1
type(t_ice),        target, save              :: ice2
integer                                       :: ierr, ix
real(kind=WP)                                 :: wtol=1.0E-15_WP

call MPI_INIT(ierr)
call par_init(partit)

! Parse command line arguments
ix = 1
do while (ix <= command_argument_count()) 
    call get_command_argument(ix, arg)
    select case (arg)
        case("-h","--help")
            if (partit%mype == 0) call print_help()
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        case("-t","--types")
            ix = ix + 1
            call get_command_argument(ix, dtype)
            ix = ix + 1
        case("-d","--directories")
            ix = ix + 1
            call get_command_argument(ix, dpath1)
            ix = ix + 1
            call get_command_argument(ix, dpath2)
            ix = ix + 1
        case("-r", "--tolerance")
            ix = ix + 1
            call get_command_argument(ix,tol)
            read(tol, *) wtol
            ix = ix + 1
        case default
            if (partit%mype == 0) print *, "wrong argument"
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    end select
end do

if (partit%mype == 0) then
    print *, "Command line arguments:"
    print *, "type       : ", trim(dtype)
    print *, "directory 1: ", trim(dpath1)
    print *, "directory 2: ", trim(dpath2)
    print *, "tolerance  : ", wtol
end if

! Check the directories exist
npepath1 = trim(dpath1)//"/np"//int_to_txt(partit%npes)
#if defined(__PGI)
INQUIRE(file=trim(npepath1), EXIST=dir_exist)
#else
INQUIRE(directory=trim(npepath1), EXIST=dir_exist)
#endif
if (.not. dir_exist) then
    if (partit%mype==0) print *, achar(27)//'[1;31m'//' -ERROR-> could not find:'//trim(npepath1)//achar(27)//'[0m'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    stop
end if

npepath2 = trim(dpath2)//"/np"//int_to_txt(partit%npes)
#if defined(__PGI)
INQUIRE(file=trim(npepath2), EXIST=dir_exist)
#else
INQUIRE(directory=trim(npepath2), EXIST=dir_exist)
#endif
if (.not. dir_exist) then
    if (partit%mype==0) print *, achar(27)//'[1;31m'//' -ERROR-> could not find:'//trim(npepath2)//achar(27)//'[0m'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    stop
end if

! Compare files
select case(trim(dtype))

    case("tracers")
        call read_all_bin_restarts(npepath1, partit=partit, mesh=mesh, tracers=tracers1)
        call read_all_bin_restarts(npepath2, partit=partit, mesh=mesh, tracers=tracers2)
        if (size(tracers1%data) .ne. size(tracers2%data)) then
            if (partit%mype==0) print *, "Inconsistent number of tracers"
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        end if
        call dderived_type(partit, tracers1, tracers2, wtol)

    case("dynamics")
        call read_all_bin_restarts(npepath1, partit=partit, mesh=mesh, dynamics=dyn1)
        call read_all_bin_restarts(npepath2, partit=partit, mesh=mesh, dynamics=dyn2)
        call dderived_type(partit, dyn1, dyn2, wtol)

    case("ice")
        call read_all_bin_restarts(npepath1, partit=partit, mesh=mesh, ice=ice1)
        call read_all_bin_restarts(npepath2, partit=partit, mesh=mesh, ice=ice2)
        call dderived_type(partit, ice1, ice2, wtol)

    case default
        if (partit%mype==0) print *, "Wrong type argument"
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)

end select

call par_ex(partit%MPI_COMM_FESOM, partit%mype)

contains

subroutine print_help()

    print '(a, /)', 'command-line options:'
    print '(a)',    '  -t, --type            derived type comparison (tracers, dynamics, ice)'
    print '(a)',    '  -d, --directories     two directories with results to compare'
    print '(a)',    '  -r, --tolerance       relative difference tolerance'
    print '(a, /)', '  -h, --help            print usage information and exit'

end subroutine print_help

end program main
