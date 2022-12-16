MODULE mod_dderived_type
USE o_PARAM
USE MOD_PARTIT
USE MOD_TRACER
USE MOD_DYN
USE MOD_ICE
USE fortran_utils
USE mod_darray
IMPLICIT NONE
PUBLIC

INTERFACE dderived_type
  MODULE PROCEDURE dderived_type_tracers
  MODULE PROCEDURE dderived_type_dynamics
  MODULE PROCEDURE dderived_type_ice
END INTERFACE dderived_type

contains

subroutine dderived_type_tracers(partit, tracers1, tracers2, wtol)
    type(t_partit), intent(in) :: partit
    type(t_tracer), intent(in) :: tracers1
    type(t_tracer), intent(in) :: tracers2
    real(kind=WP)              :: wtol
    integer                    :: i, fileunit

    call darray_set_tolerance(wtol)

    fileunit = partit%mype+300
    open(newunit = fileunit, &
         file     = "output/d_tracer."//mpirank_to_txt(partit%MPI_COMM_FESOM), &
         status   = 'replace')

    do i=1,size(tracers1%data)
        call diff_array(fileunit,                                  &
                        tracers1%data(i)%values,                   &
                        tracers2%data(i)%values,                   &
                        "tracers_data_"//int_to_txt(i)//"_values"  )
        call diff_array(fileunit,                                  &
                        tracers1%data(i)%valuesAB,                 &
                        tracers2%data(i)%valuesAB,                 &
                        "tracers_data_"//int_to_txt(i)//"_valuesAB")
    end do

    call diff_array(fileunit,                       &
                    tracers1%work%del_ttf,          &
                    tracers2%work%del_ttf,          &
                    "tracers_work_del_ttf"          )

    call diff_array(fileunit,                       &
                    tracers1%work%del_ttf_advhoriz, &
                    tracers2%work%del_ttf_advhoriz, &
                    "tracers_work_del_ttf_advhoriz" )

    call diff_array(fileunit,                       &
                    tracers1%work%del_ttf_advvert,  &
                    tracers2%work%del_ttf_advvert,  &
                    "tracers_work_del_ttf_advvert"  )

    call diff_array(fileunit,                       &
                    tracers1%work%tr_dvd_horiz,     &
                    tracers2%work%tr_dvd_horiz,     &
                    "tracers_work_tr_dvd_horiz"     )

    call diff_array(fileunit,                       &
                    tracers1%work%tr_dvd_vert,      &
                    tracers2%work%tr_dvd_vert,      &
                    "tracers_work_tr_dvd_vert"      )

    call diff_array(fileunit,                       &
                    tracers1%work%fct_LO,           &
                    tracers2%work%fct_LO,           &
                    "tracers_work_fct_LO"           )

    call diff_array(fileunit,                       &
                    tracers1%work%adv_flux_hor,     &
                    tracers2%work%adv_flux_hor,     &
                    "tracers_work_adv_flux_hor"     )

    call diff_array(fileunit,                       &
                    tracers1%work%adv_flux_ver,     &
                    tracers2%work%adv_flux_ver,     &
                    "tracers_work_adv_flux_ver"     )

    call diff_array(fileunit,                       &
                    tracers1%work%fct_ttf_max,      &
                    tracers2%work%fct_ttf_max,      &
                    "tracers_work_fct_ttf_max"      )

    call diff_array(fileunit,                       &
                    tracers1%work%fct_ttf_min,      &
                    tracers2%work%fct_ttf_min,      &
                    "tracers_work_fct_ttf_min"      )

    call diff_array(fileunit,                       &
                    tracers1%work%fct_plus,         &
                    tracers2%work%fct_plus,         &
                    "tracers_work_fct_plus"         )

    call diff_array(fileunit,                       &
                    tracers1%work%fct_minus,        &
                    tracers2%work%fct_minus,        &
                    "tracers_work_fct_minus")

    call diff_array(fileunit,                       &
                    tracers1%work%edge_up_dn_grad,  &
                    tracers2%work%edge_up_dn_grad,  &
                    "tracers_work_edge_up_dn_grad")

    close(fileunit)

end subroutine dderived_type_tracers


subroutine dderived_type_dynamics(partit, dyn1, dyn2, wtol)
    type(t_partit), intent(in) :: partit
    type(t_dyn), intent(in)    :: dyn1
    type(t_dyn), intent(in)    :: dyn2
    real(kind=WP)              :: wtol
    integer                    :: i, fileunit

    call darray_set_tolerance(wtol)

    fileunit = partit%mype+300
    open(newunit = fileunit, &
         file     = "output/d_dyn."//mpirank_to_txt(partit%MPI_COMM_FESOM), &
         status   = 'new')


    close(fileunit)

end subroutine dderived_type_dynamics


subroutine dderived_type_ice(partit, ice1, ice2, wtol)
    type(t_partit), intent(in) :: partit
    type(t_ice), intent(in)    :: ice1
    type(t_ice), intent(in)    :: ice2
    real(kind=WP)              :: wtol
    integer                    :: i, fileunit

    call darray_set_tolerance(wtol)

    fileunit = partit%mype+300
    open(newunit = fileunit, &
         file     = "output/d_ice."//mpirank_to_txt(partit%MPI_COMM_FESOM), &
         status   = 'new')


    close(fileunit)

end subroutine dderived_type_ice

END MODULE mod_dderived_type
