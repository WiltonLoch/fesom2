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
                    tracers1%work%nboundary_lay,    &
                    tracers2%work%nboundary_lay,    &
                    "tracers_work_nboundary_lay")

    call diff_array(fileunit,                       &
                    tracers1%work%edge_up_dn_tri,   &
                    tracers2%work%edge_up_dn_tri,   &
                    "tracers_work_edge_up_dn_tri")

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
         status   = 'replace')


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
         status   = 'replace')

    call diff_array(fileunit,                       &
                    ice1%uice,  &
                    ice2%uice,  &
                    "ice_uice")

    call diff_array(fileunit,                       &
                    ice1%uice_rhs,  &
                    ice2%uice_rhs,  &
                    "ice_uice_rhs")

    call diff_array(fileunit,                       &
                    ice1%uice_old,  &
                    ice2%uice_old,  &
                    "ice_uice_old")

    call diff_array(fileunit,                       &
                    ice1%uice_aux,  &
                    ice2%uice_aux,  &
                    "ice_uice_aux")

    call diff_array(fileunit,                       &
                    ice1%vice,  &
                    ice2%vice,  &
                    "ice_vice")

    call diff_array(fileunit,                       &
                    ice1%vice_rhs,  &
                    ice2%vice_rhs,  &
                    "ice_vice_rhs")

    call diff_array(fileunit,                       &
                    ice1%vice_old,  &
                    ice2%vice_old,  &
                    "ice_vice_old")

    call diff_array(fileunit,                       &
                    ice1%vice_aux,  &
                    ice2%vice_aux,  &
                    "ice_vice_aux")

    call diff_array(fileunit,                       &
                    ice1%stress_atmice_x,  &
                    ice2%stress_atmice_x,  &
                    "ice_stress_atmice_x")

    call diff_array(fileunit,                       &
                    ice1%stress_iceoce_x,  &
                    ice2%stress_iceoce_x,  &
                    "ice_stress_iceoce_x")

    call diff_array(fileunit,                       &
                    ice1%stress_atmice_y,  &
                    ice2%stress_atmice_y,  &
                    "ice_stress_atmice_y")

    call diff_array(fileunit,                       &
                    ice1%stress_iceoce_y,  &
                    ice2%stress_iceoce_y,  &
                    "ice_stress_iceoce_y")

    call diff_array(fileunit,                       &
                    ice1%srfoce_temp,  &
                    ice2%srfoce_temp,  &
                    "ice_srfoce_temp")

    call diff_array(fileunit,                       &
                    ice1%srfoce_salt,  &
                    ice2%srfoce_salt,  &
                    "ice_srfoce_salt")

    call diff_array(fileunit,                       &
                    ice1%srfoce_ssh,  &
                    ice2%srfoce_ssh,  &
                    "ice_srfoce_ssh")

    call diff_array(fileunit,                       &
                    ice1%srfoce_u,  &
                    ice2%srfoce_u,  &
                    "ice_srfoce_u")

    call diff_array(fileunit,                       &
                    ice1%srfoce_v,  &
                    ice2%srfoce_v,  &
                    "ice_srfoce_v")

    call diff_array(fileunit,                       &
                    ice1%flx_fw,  &
                    ice2%flx_fw,  &
                    "ice_flx_fw")

    call diff_array(fileunit,                       &
                    ice1%flx_h,  &
                    ice2%flx_h,  &
                    "ice_flx_h")

    call diff_array(fileunit,                       &
                    ice1%alpha_evp_array,  &
                    ice2%alpha_evp_array,  &
                    "ice_alpha_evp_array")

    call diff_array(fileunit,                       &
                    ice1%beta_evp_array,  &
                    ice2%beta_evp_array,  &
                    "ice_beta_evp_array")

    do i=1,size(ice1%data)
           call diff_array(fileunit,                       &
                           ice1%data(i)%values,  &
                           ice2%data(i)%values,  &
                           "ice_data_"//int_to_txt(i)//"_values") 

           call diff_array(fileunit,                       &
                           ice1%data(i)%values_old,  &
                           ice2%data(i)%values_old,  &
                           "ice_data_"//int_to_txt(i)//"_values_old")

           call diff_array(fileunit,                       &
                           ice1%data(i)%values_rhs,  &
                           ice2%data(i)%values_rhs,  &
                           "ice_data_"//int_to_txt(i)//"_values_rhs")

           call diff_array(fileunit,                       &
                           ice1%data(i)%values_div_rhs,  &
                           ice2%data(i)%values_div_rhs,  &
                           "ice_data_"//int_to_txt(i)//"_values_div_rhs")

           call diff_array(fileunit,                       &
                           ice1%data(i)%dvalues,  &
                           ice2%data(i)%dvalues,  &
                           "ice_data_"//int_to_txt(i)//"_dvalues")

           call diff_array(fileunit,                       &
                           ice1%data(i)%valuesl,  &
                           ice2%data(i)%valuesl,  &
                           "ice_data_"//int_to_txt(i)//"_valuesl")
    end do

    call diff_array(fileunit,                       &
                    ice1%work%fct_tmax,  &
                    ice2%work%fct_tmax,  &
                    "ice_work_fct_tmax")

    call diff_array(fileunit,                       &
                    ice1%work%fct_tmin,  &
                    ice2%work%fct_tmin,  &
                    "ice_work_fct_tmin")

    call diff_array(fileunit,                       &
                    ice1%work%fct_plus,  &
                    ice2%work%fct_plus,  &
                    "ice_work_fct_plus")

    call diff_array(fileunit,                       &
                    ice1%work%fct_minus,  &
                    ice2%work%fct_minus,  &
                    "ice_work_fct_minus")

    call diff_array(fileunit,                       &
                    ice1%work%fct_fluxes,  &
                    ice2%work%fct_fluxes,  &
                    "ice_work_fct_fluxes")

    call diff_array(fileunit,                       &
                    ice1%work%fct_massmatrix,  &
                    ice2%work%fct_massmatrix,  &
                    "ice_work_fct_massmatrix")

    call diff_array(fileunit,                       &
                    ice1%work%sigma11,  &
                    ice2%work%sigma11,  &
                    "ice_work_sigma11")

    call diff_array(fileunit,                       &
                    ice1%work%sigma12,  &
                    ice2%work%sigma12,  &
                    "ice_work_sigma12")

    call diff_array(fileunit,                       &
                    ice1%work%sigma22,  &
                    ice2%work%sigma22,  &
                    "ice_work_sigma22")

    call diff_array(fileunit,                       &
                    ice1%work%eps11,  &
                    ice2%work%eps11,  &
                    "ice_work_eps11")

    call diff_array(fileunit,                       &
                    ice1%work%eps12,  &
                    ice2%work%eps12,  &
                    "ice_work_eps12")

    call diff_array(fileunit,                       &
                    ice1%work%eps22,  &
                    ice2%work%eps22,  &
                    "ice_work_eps22")

    call diff_array(fileunit,                       &
                    ice1%work%ice_strength,  &
                    ice2%work%ice_strength,  &
                    "ice_work_ice_strength")

    call diff_array(fileunit,                       &
                    ice1%work%inv_areamass,  &
                    ice2%work%inv_areamass,  &
                    "ice_work_inv_areamass")

    call diff_array(fileunit,                       &
                    ice1%work%inv_mass,  &
                    ice2%work%inv_mass,  &
                    "ice_work_inv_mass")

    call diff_array(fileunit,                       &
                    ice1%thermo%t_skin,  &
                    ice2%thermo%t_skin,  &
                    "ice_thermo_t_skin")

    call diff_array(fileunit,                       &
                    ice1%thermo%thdgr,  &
                    ice2%thermo%thdgr,  &
                    "ice_thermo_thdgr")

    call diff_array(fileunit,                       &
                    ice1%thermo%thdgrsn,  &
                    ice2%thermo%thdgrsn,  &
                    "ice_thermo_thdgrsn")

    call diff_array(fileunit,                       &
                    ice1%thermo%thdgr_old,  &
                    ice2%thermo%thdgr_old,  &
                    "ice_thermo_thdgr_old")

    call diff_array(fileunit,                       &
                    ice1%thermo%ustar,  &
                    ice2%thermo%ustar,  &
                    "ice_thermo_ustar")

#if defined (__oasis) || defined (__ifsinterface)
    call diff_array(fileunit,                       &
                    ice1%atmcoupl%oce_flx_h,  &
                    ice2%atmcoupl%oce_flx_h,  &
                    "ice_atmcoupl_oce_flx_h")

    call diff_array(fileunit,                       &
                    ice1%atmcoupl%ice_flx_h,  &
                    ice2%atmcoupl%ice_flx_h,  &
                    "ice_atmcoupl_ice_flx_h")

    call diff_array(fileunit,                       &
                    ice1%atmcoupl%tmpoce_flx_h,  &
                    ice2%atmcoupl%tmpoce_flx_h,  &
                    "ice_atmcoupl_tmpoce_flx_h")

    call diff_array(fileunit,                       &
                    ice1%atmcoupl%tmpice_flx_h,  &
                    ice2%atmcoupl%tmpice_flx_h,  &
                    "ice_atmcoupl_tmpice_flx_h")

#if defined (__oifs) || defined (__ifsinterface)
        call diff_array(fileunit,                       &
                    ice1%atmcoupl%ice_alb,  &
                    ice2%atmcoupl%ice_alb,  &
                    "ice_atmcoupl_ice_alb")

        call diff_array(fileunit,                       &
                    ice1%atmcoupl%enthalpyoffuse,  &
                    ice2%atmcoupl%enthalpyoffuse,  &
                    "ice_atmcoupl_enthalpyoffuse")
#endif
#endif

    close(fileunit)

end subroutine dderived_type_ice

END MODULE mod_dderived_type
