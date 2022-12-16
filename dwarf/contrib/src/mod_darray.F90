MODULE mod_darray
USE o_PARAM
IMPLICIT NONE
PRIVATE

PUBLIC :: diff_array, darray_set_tolerance

real(KIND=WP), private :: tol

INTERFACE diff_array
  MODULE PROCEDURE diff_array1D_wp
  MODULE PROCEDURE diff_array2D_wp
  MODULE PROCEDURE diff_array3D_wp
END INTERFACE diff_array

contains

subroutine darray_set_tolerance(wtol)
    real(kind=WP), intent(in) :: wtol

    tol = wtol
end subroutine darray_set_tolerance

subroutine diff_array1D_wp(fileunit, a, b, vname)
    integer, intent(in)                                       :: fileunit
    real(kind=WP), intent(in), allocatable, dimension(:)   :: a
    real(kind=WP), intent(in), allocatable, dimension(:)   :: b
    character(LEN=*), intent(in)                              :: vname

    integer, dimension(1)                                     :: sa, sb
    integer                                                   :: i

end subroutine diff_array1D_wp

subroutine diff_array2D_wp(fileunit, a, b, vname)
    integer, intent(in)                                       :: fileunit
    real(kind=WP), intent(in), allocatable, dimension(:,:) :: a
    real(kind=WP), intent(in), allocatable, dimension(:,:) :: b
    character(LEN=*), intent(in)                              :: vname

    integer, dimension(2)                                     :: sa, sb
    integer                                                   :: i,j

    if (allocated(a) .and. allocated(b)) then
        sa = SHAPE(a)
        sb = SHAPE(b)
        if (all(sa == sb)) then
            do j = 1,sa(2)
                do i = 1,sa(1)
                    if (abs((a(i,j)-b(i,j))/a(i,j)) > tol) then
                        write(fileunit, '(a,a,i8,a,i8,a,e15.7,a,e15.7,a,e15.7,a,e15.7)') &
                              trim(vname)," i=",i," j=",j," v1=",a(i,j)," v2=",b(i,j)," abs=",abs(a(i,j)-b(i,j))," rel=",abs((a(i,j)-b(i,j))/a(i,j))
                    end if
                end do
            end do
        end if
    end if

end subroutine diff_array2D_wp

subroutine diff_array3D_wp(fileunit, a, b, vname)
    integer, intent(in)                                         :: fileunit
    real(kind=WP), intent(in), allocatable, dimension(:,:,:) :: a
    real(kind=WP), intent(in), allocatable, dimension(:,:,:) :: b
    character(LEN=*), intent(in)                                :: vname

    integer, dimension(3)                                       :: sa, sb
    integer                                                     :: i,j,k


end subroutine diff_array3D_wp

END MODULE mod_darray
