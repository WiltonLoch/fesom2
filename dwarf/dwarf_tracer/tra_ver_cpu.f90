# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90"
module oce_adv_tra_ver_interfaces
  interface
! implicit 1st order upwind vertical advection with to solve for fct_LO
! updates the input tracer ttf
    subroutine adv_tra_vert_impl(dt, w, ttf, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      real(kind=WP), intent(in), target  :: dt
      type(t_partit),intent(in), target  :: partit
      type(t_mesh),  intent(in), target  :: mesh
      real(kind=WP), intent(inout)       :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)          :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    end subroutine
!===============================================================================
! 1st order upwind (explicit)
! returns flux given at vertical interfaces of scalar volumes
! IF o_init_zero=.TRUE.  : flux will be set to zero before computation
! IF o_init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_ver_upw1(w, ttf, partit, mesh, flux, o_init_zero)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_partit),intent(in), target :: partit
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
      logical, optional                 :: o_init_zero
    end subroutine
!===============================================================================
! QR (4th order centerd)
! returns flux given at vertical interfaces of scalar volumes
! IF o_init_zero=.TRUE.  : flux will be set to zero before computation
! IF o_init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_ver_qr4c(w, ttf, partit, mesh, num_ord, flux, o_init_zero)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_partit),intent(in), target :: partit
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
      logical, optional                 :: o_init_zero
    end subroutine
!===============================================================================
! Vertical advection with PPM reconstruction (5th order)
! returns flux given at vertical interfaces of scalar volumes
! IF o_init_zero=.TRUE.  : flux will be set to zero before computation
! IF o_init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
   subroutine adv_tra_vert_ppm(dt, w, ttf, partit, mesh, flux, o_init_zero)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      real(kind=WP), intent(in), target :: dt
      type(t_partit),intent(in), target :: partit
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP)                     :: tvert(mesh%nl), tv
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
      logical, optional                 :: o_init_zero
    end subroutine
! central difference reconstruction (2nd order, use only with FCT)
! returns flux given at vertical interfaces of scalar volumes
! IF o_init_zero=.TRUE.  : flux will be set to zero before computation
! IF o_init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_ver_cdiff(w, ttf, partit, mesh, flux, o_init_zero)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_partit),intent(in), target :: partit
      type(t_mesh),  intent(in), target :: mesh
      integer                           :: n, nz, nl1
      real(kind=WP)                     :: tvert(mesh%nl), tv
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
      logical, optional                 :: o_init_zero
    end subroutine
  end interface
end module
!===============================================================================
subroutine adv_tra_vert_impl(dt, w, ttf, partit, mesh)
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto

    implicit none
    real(kind=WP), intent(in) , target :: dt
    type(t_partit),intent(in), target  :: partit
    type(t_mesh),  intent(in) , target :: mesh
    real(kind=WP), intent(inout)       :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)          :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP)                      :: a(mesh%nl), b(mesh%nl), c(mesh%nl), tr(mesh%nl)
    real(kind=WP)                      :: cp(mesh%nl), tp(mesh%nl)
    real(kind=WP)                      :: zbar_n(mesh%nl), z_n(mesh%nl-1)
    integer                            :: nz, n, nzmax, nzmin
    real(kind=WP)                      :: m, zinv, dt_inv, dz
    real(kind=WP)                      :: c1, v_adv


# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_def.h" 1

  integer,          pointer     :: MPI_COMM_FESOM ! FESOM communicator (for ocean only runs if often a copy of MPI_COMM_WORLD)
  type(com_struct), pointer     :: com_nod2D
  type(com_struct), pointer     :: com_elem2D
  type(com_struct), pointer     :: com_elem2D_full
  integer                       :: ub, lb ! to work with r(s)_mpitype_elem3D(nod3D)

  integer, dimension(:),     pointer :: s_mpitype_edge2D, r_mpitype_edge2D
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D, r_mpitype_elem2D
  integer, dimension(:),     pointer :: s_mpitype_elem2D_full_i, r_mpitype_elem2D_full_i
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D_full,   r_mpitype_elem2D_full
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D,        r_mpitype_elem3D
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D_full,   r_mpitype_elem3D_full

  integer, dimension(:),     pointer :: s_mpitype_nod2D, r_mpitype_nod2D
  integer, dimension(:),     pointer :: s_mpitype_nod2D_i, r_mpitype_nod2D_i
  integer, dimension(:,:,:), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

  integer, pointer   :: MPIERR
  integer, pointer   :: npes
  integer, pointer   :: mype
  integer, pointer   :: maxPEnum

  integer, dimension(:), pointer     :: part

  ! Mesh partition
  integer, pointer                :: myDim_nod2D, eDim_nod2D
  integer, dimension(:), pointer  :: myList_nod2D
  integer, pointer                :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, dimension(:), pointer  :: myList_elem2D
  integer, pointer                :: myDim_edge2D, eDim_edge2D
  integer, dimension(:), pointer  :: myList_edge2D

  integer, pointer                :: pe_status

  integer, dimension(:), pointer  ::  remPtr_nod2D(:),  remList_nod2D(:)
  integer, dimension(:), pointer  ::  remPtr_elem2D(:), remList_elem2D(:)

  logical, pointer                :: elem_full_flag
# 111 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_def.h" 1
integer         , pointer :: nod2D
integer         , pointer :: elem2D
integer         , pointer :: edge2D
integer         , pointer :: edge2D_in
real(kind=WP)   , pointer :: ocean_area
real(kind=WP)   , pointer :: ocean_areawithcav
integer         , pointer :: nl
integer         , pointer :: nn_size
real(kind=WP), dimension(:,:), pointer :: coord_nod2D, geo_coord_nod2D
integer, dimension(:,:)      , pointer :: elem2D_nodes
integer, dimension(:,:)      , pointer :: edges
integer, dimension(:,:)      , pointer :: edge_tri
integer, dimension(:,:)      , pointer :: elem_edges
real(kind=WP), dimension(:)  , pointer :: elem_area
real(kind=WP), dimension(:,:), pointer :: edge_dxdy, edge_cross_dxdy
real(kind=WP), dimension(:)  , pointer :: elem_cos, metric_factor
integer,       dimension(:,:), pointer :: elem_neighbors
integer,       dimension(:,:), pointer :: nod_in_elem2D
real(kind=WP), dimension(:,:), pointer :: x_corners, y_corners
integer,       dimension(:)  , pointer :: nod_in_elem2D_num
real(kind=WP), dimension(:)  , pointer :: depth
real(kind=WP), dimension(:,:), pointer :: gradient_vec
real(kind=WP), dimension(:,:), pointer :: gradient_sca
integer,       dimension(:)  , pointer :: bc_index_nod2D
real(kind=WP), dimension(:)  , pointer :: zbar, Z, elem_depth
integer,       dimension(:)  , pointer :: nlevels, nlevels_nod2D, nlevels_nod2D_min
real(kind=WP), dimension(:,:), pointer :: area, area_inv, areasvol, areasvol_inv
real(kind=WP), dimension(:)  , pointer :: mesh_resolution
real(kind=WP), dimension(:)  , pointer :: lump2d_north, lump2d_south
type(sparse_matrix)          , pointer :: ssh_stiff
integer,       dimension(:)  , pointer :: cavity_flag_n, cavity_flag_e
real(kind=WP), dimension(:)  , pointer :: cavity_depth
integer,       dimension(:)  , pointer :: ulevels, ulevels_nod2D, ulevels_nod2D_max
integer,       dimension(:)  , pointer :: nn_num
integer,       dimension(:,:), pointer :: nn_pos

real(kind=WP), dimension(:,:), pointer :: hnode
real(kind=WP), dimension(:,:), pointer :: hnode_new
real(kind=WP), dimension(:,:), pointer :: zbar_3d_n
real(kind=WP), dimension(:,:), pointer :: Z_3d_n
real(kind=WP), dimension(:,:), pointer :: helem
real(kind=WP), dimension(:)  , pointer :: bottom_elem_thickness
real(kind=WP), dimension(:)  , pointer :: bottom_node_thickness
real(kind=WP), dimension(:)  , pointer :: dhe
real(kind=WP), dimension(:)  , pointer :: hbar
real(kind=WP), dimension(:)  , pointer :: hbar_old
!real(kind=WP), dimension(:)  , pointer :: zbar_n
!real(kind=WP), dimension(:)  , pointer :: Z_n
real(kind=WP), dimension(:)  , pointer :: zbar_n_bot
real(kind=WP), dimension(:)  , pointer :: zbar_e_bot
real(kind=WP), dimension(:)  , pointer :: zbar_n_srf
real(kind=WP), dimension(:)  , pointer :: zbar_e_srf
# 112 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_ass.h" 1
MPI_COMM_FESOM  => partit%MPI_COMM_FESOM
com_nod2D       => partit%com_nod2D
com_elem2D      => partit%com_elem2D
com_elem2D_full => partit%com_elem2D_full
myDim_nod2D     => partit%myDim_nod2D
eDim_nod2D      => partit%eDim_nod2D
myDim_elem2D    => partit%myDim_elem2D
eDim_elem2D     => partit%eDim_elem2D
eXDim_elem2D    => partit%eXDim_elem2D
myDim_edge2D    => partit%myDim_edge2D
eDim_edge2D     => partit%eDim_edge2D
pe_status       => partit%pe_status
elem_full_flag  => partit%elem_full_flag
MPIERR          => partit%MPIERR
npes            => partit%npes
mype            => partit%mype
maxPEnum        => partit%maxPEnum
part            => partit%part

lb=lbound(partit%s_mpitype_elem3D, 2)
ub=ubound(partit%s_mpitype_elem3D, 2)

myList_nod2D (1:myDim_nod2D +eDim_nod2D)                => partit%myList_nod2D(:)
myList_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)  => partit%myList_elem2D(:)
myList_edge2D(1:myDim_edge2D+eDim_edge2D)               => partit%myList_edge2D(:)

if (allocated(partit%remPtr_nod2D)) then
   remPtr_nod2D  (1:npes)                => partit%remPtr_nod2D(:)
   remList_nod2D (1:remPtr_nod2D(npes))  => partit%remList_nod2D(:)
end if

if (allocated(partit%remPtr_elem2D)) then
remPtr_elem2D (1:npes)                => partit%remPtr_elem2D(:)
remList_elem2D(1:remPtr_elem2D(npes)) => partit%remList_elem2D(:)
end if

s_mpitype_elem2D(1:com_elem2D%sPEnum, 1:4) => partit%s_mpitype_elem2D(:,:)
r_mpitype_elem2D(1:com_elem2D%rPEnum, 1:4) => partit%r_mpitype_elem2D(:,:)

s_mpitype_elem2D_full_i(1:com_elem2D_full%sPEnum) => partit%s_mpitype_elem2D_full_i(:)
r_mpitype_elem2D_full_i(1:com_elem2D_full%rPEnum) => partit%r_mpitype_elem2D_full_i(:)

s_mpitype_elem2D_full(1:com_elem2D_full%sPEnum, 1:4) => partit%s_mpitype_elem2D_full(:,:)
r_mpitype_elem2D_full(1:com_elem2D_full%rPEnum, 1:4) => partit%r_mpitype_elem2D_full(:,:)

s_mpitype_elem3D(1:com_elem2D%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D(:,:,:)
r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D(:,:,:)

s_mpitype_elem3D_full(1:com_elem2D_full%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D_full(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4)           => partit%r_mpitype_elem3D(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

s_mpitype_nod2D(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D(:)
r_mpitype_nod2D(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D(:)

s_mpitype_nod2D_i(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D_i(:)
r_mpitype_nod2D_i(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D_i(:)

s_mpitype_nod3D(1:com_nod2D%sPEnum, lb:ub, 1:3) => partit%s_mpitype_nod3D(:,:,:)
r_mpitype_nod3D(1:com_nod2D%rPEnum, lb:ub, 1:3) => partit%r_mpitype_nod3D(:,:,:)
# 113 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_ass.h" 1
nod2D              => mesh%nod2D
elem2D             => mesh%elem2D
edge2D             => mesh%edge2D
edge2D_in          => mesh%edge2D_in
ocean_area         => mesh%ocean_area
nl                 => mesh%nl
nn_size            => mesh%nn_size
ocean_areawithcav  => mesh%ocean_areawithcav
coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)                  => mesh%coord_nod2D(:,:)
geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)              => mesh%geo_coord_nod2D(:,:)
elem2D_nodes(1:3, 1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%elem2D_nodes(:,:)
edges(1:2,1:myDim_edge2D+eDim_edge2D)                      => mesh%edges(:,:)
edge_tri(1:2,1:myDim_edge2D+eDim_edge2D)                   => mesh%edge_tri(:,:)
elem_edges(1:3,1:myDim_elem2D)                             => mesh%elem_edges(:,:)
elem_area(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)         => mesh%elem_area(:)
edge_dxdy(1:2,1:myDim_edge2D+eDim_edge2D)                  => mesh%edge_dxdy(:,:)
edge_cross_dxdy(1:4,1:myDim_edge2D+eDim_edge2D)            => mesh%edge_cross_dxdy(:,:)
elem_cos(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)          => mesh%elem_cos(:)
metric_factor(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%metric_factor(:)
elem_neighbors(1:3,1:myDim_elem2D)                         => mesh%elem_neighbors(:,:)
nod_in_elem2D      => mesh%nod_in_elem2D   ! (maxval(rmax),myDim_nod2D+eDim_nod2D)
x_corners          => mesh%x_corners   ! (myDim_nod2D, maxval(rmax))
y_corners          => mesh%y_corners   ! (myDim_nod2D, maxval(rmax))
nod_in_elem2D_num(1:myDim_nod2D+eDim_nod2D)                => mesh%nod_in_elem2D_num(:)
depth(1:myDim_nod2D+eDim_nod2D)                            => mesh%depth(:)
gradient_vec(1:6,1:myDim_elem2D)                           => mesh%gradient_vec(:,:)
gradient_sca(1:6,1:myDim_elem2D)                           => mesh%gradient_sca(:,:)
bc_index_nod2D(1:myDim_nod2D+eDim_nod2D)                   => mesh%bc_index_nod2D(:)
zbar(1:mesh%nl)                                            => mesh%zbar(:)
Z(1:mesh%nl-1)                                             => mesh%Z(:)
elem_depth         => mesh%elem_depth      ! never used, not even allocated
nlevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%nlevels(:)
nlevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%nlevels_nod2D(:)
nlevels_nod2D_min(1:myDim_nod2D+eDim_nod2D)                => mesh%nlevels_nod2D_min(:)
area(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)                   => mesh%area(:,:)
areasvol(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%areasvol(:,:)
area_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%area_inv(:,:)
areasvol_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)           => mesh%areasvol_inv(:,:)
mesh_resolution(1:myDim_nod2d+eDim_nod2D)                  => mesh%mesh_resolution(:)
ssh_stiff                                                  => mesh%ssh_stiff
lump2d_north(1:myDim_nod2d)                                => mesh%lump2d_north(:)
lump2d_south(1:myDim_nod2d)                                => mesh%lump2d_south(:)
cavity_flag_n(1:myDim_nod2D+eDim_nod2D)                    => mesh%cavity_flag_n(:)
cavity_flag_e(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%cavity_flag_e(:)
!!$cavity_lev_nod2D(1:myDim_nod2D+eDim_nod2D)                 => mesh%cavity_lev_nod2D
!!$cavity_lev_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%cavity_lev_elem2D
cavity_depth(1:myDim_nod2D+eDim_nod2D)                     => mesh%cavity_depth(:)
ulevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%ulevels(:)
ulevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%ulevels_nod2D(:)
ulevels_nod2D_max(1:myDim_nod2D+eDim_nod2D)                => mesh%ulevels_nod2D_max(:)
nn_num(1:myDim_nod2D)                                      => mesh%nn_num(:)
nn_pos(1:mesh%nn_size, 1:myDim_nod2D)                      => mesh%nn_pos(:,:)
hnode(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)               => mesh%hnode(:,:)
hnode_new(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)           => mesh%hnode_new(:,:)
zbar_3d_n(1:mesh%nl, 1:myDim_nod2D+eDim_nod2D)             => mesh%zbar_3d_n(:,:)
Z_3d_n(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)              => mesh%Z_3d_n(:,:)
helem(1:mesh%nl-1, 1:myDim_elem2D)                         => mesh%helem(:,:)
bottom_elem_thickness(1:myDim_elem2D)                      => mesh%bottom_elem_thickness(:)
bottom_node_thickness(1:myDim_nod2D+eDim_nod2D)            => mesh%bottom_node_thickness(:)
dhe(1:myDim_elem2D)                                        => mesh%dhe(:)
hbar(1:myDim_nod2D+eDim_nod2D)                             => mesh%hbar(:)
hbar_old(1:myDim_nod2D+eDim_nod2D)                         => mesh%hbar_old(:)
!zbar_n(1:mesh%nl)                                          => mesh%zbar_n
!Z_n(1:mesh%nl-1)                                           => mesh%Z_n
zbar_n_bot(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_bot(:)
zbar_e_bot(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_bot(:)
zbar_n_srf(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_srf(:)
zbar_e_srf(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_srf(:)
# 114 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

    dt_inv=1.0_WP/dt
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a, b, c, tr, cp, tp, n, nz, nzmax, nzmin, m, zinv, dz, c1, v_adv)
!$OMP DO
    !___________________________________________________________________________
    ! loop over local nodes
    do n=1,myDim_nod2D

        ! initialise
        a  = 0.0_WP
        b  = 0.0_WP
        c  = 0.0_WP
        tr = 0.0_WP
        tp = 0.0_WP
        cp = 0.0_WP

        ! max. number of levels at node n
        nzmax=nlevels_nod2D(n)

        ! upper surface index, in case of cavity !=1
        nzmin=ulevels_nod2D(n)

        !___________________________________________________________________________
        ! Here can not exchange zbar_n & Z_n with zbar_3d_n & Z_3d_n because
        ! they be calculate from the actualized mesh with hnode_new
        ! calculate new zbar (depth of layers) and Z (mid depths of layers)
        ! depending on layer thinkness over depth at node n
        ! Be carefull here vertical operation have to be done on NEW vertical mesh !!!
        zbar_n=0.0_WP
        Z_n=0.0_WP
        zbar_n(nzmax)=zbar_n_bot(n)
        Z_n(nzmax-1) =zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP
        do nz=nzmax-1,nzmin+1,-1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
            Z_n(nz-1)  = zbar_n(nz)   + hnode_new(nz-1,n)/2.0_WP
        end do
        zbar_n(nzmin) = zbar_n(nzmin+1) + hnode_new(nzmin,n)

        !_______________________________________________________________________
        ! Regular part of coefficients: --> surface layer
        nz=nzmin

        ! 1/dz(nz)
        zinv=1.0_WP*dt    ! no .../(zbar(1)-zbar(2)) because of  ALE

        !!PS a(nz)=0.0_WP
        !!PS v_adv=zinv*areasvol(nz+1,n)/areasvol(nz,n)
        !!PS b(nz)= hnode_new(nz,n)+W(nz, n)*zinv-min(0._WP, W(nz+1, n))*v_adv
        !!PS c(nz)=-max(0._WP, W(nz+1, n))*v_adv

        a(nz)=0.0_WP
        v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
        b(nz)= hnode_new(nz,n)+W(nz, n)*v_adv

        v_adv=zinv*area(nz+1,n)/areasvol(nz,n)
        b(nz)= b(nz)-min(0._WP, W(nz+1, n))*v_adv
        c(nz)=-max(0._WP, W(nz+1, n))*v_adv

        !_______________________________________________________________________
        ! Regular part of coefficients: --> 2nd...nl-2 layer
        do nz=nzmin+1, nzmax-2
            ! update from the vertical advection
            v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
            a(nz)=min(0._WP, W(nz, n))*v_adv
            b(nz)=hnode_new(nz,n)+max(0._WP, W(nz, n))*v_adv

            v_adv=zinv*area(nz+1,n)/areasvol(nz,n)
            b(nz)=b(nz)-min(0._WP, W(nz+1, n))*v_adv
            c(nz)=     -max(0._WP, W(nz+1, n))*v_adv
        end do ! --> do nz=2, nzmax-2

        !_______________________________________________________________________
        ! Regular part of coefficients: --> nl-1 layer
        nz=nzmax-1
        ! update from the vertical advection
        !!PS a(nz)=                min(0._WP, W(nz, n))*zinv
        !!PS b(nz)=hnode_new(nz,n)+max(0._WP, W(nz, n))*zinv
        !!PS c(nz)=0.0_WP
        v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
        a(nz)=                min(0._WP, W(nz, n))*v_adv
        b(nz)=hnode_new(nz,n)+max(0._WP, W(nz, n))*v_adv
        c(nz)=0.0_WP

        !_______________________________________________________________________
        nz=nzmin
        dz=hnode_new(nz,n) ! It would be (zbar(nz)-zbar(nz+1)) if not ALE
        tr(nz)=-(b(nz)-dz)*ttf(nz,n)-c(nz)*ttf(nz+1,n)

        do nz=nzmin+1,nzmax-2
            dz=hnode_new(nz,n)
            tr(nz)=-a(nz)*ttf(nz-1,n)-(b(nz)-dz)*ttf(nz,n)-c(nz)*ttf(nz+1,n)
        end do
        nz=nzmax-1
        dz=hnode_new(nz,n)
        tr(nz)=-a(nz)*ttf(nz-1,n)-(b(nz)-dz)*ttf(nz,n)

        !_______________________________________________________________________
        nz = nzmin
        cp(nz) = c(nz)/b(nz)
        tp(nz) = tr(nz)/b(nz)

        ! solve for vectors c-prime and t, s-prime
        do nz = nzmin+1,nzmax-1
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
        end do

        !_______________________________________________________________________
        ! start with back substitution
        tr(nzmax-1) = tp(nzmax-1)

        ! solve for x from the vectors c-prime and d-prime
        do nz = nzmax-2, nzmin, -1
            tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
        end do

        !_______________________________________________________________________
        ! update tracer
        do nz=nzmin,nzmax-1
            ttf(nz,n)=ttf(nz,n)+tr(nz)
        end do
    end do ! --> do n=1,myDim_nod2D
!$OMP END DO
!$OMP END PARALLEL
end subroutine adv_tra_vert_impl
!
!
!===============================================================================
subroutine adv_tra_ver_upw1(w, ttf, partit, mesh, flux, o_init_zero)
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto

    implicit none
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP)                     :: tvert(mesh%nl)
    integer                           :: n, nz, nzmax, nzmin
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
    logical, optional                 :: o_init_zero
    logical                           :: l_init_zero

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_def.h" 1

  integer,          pointer     :: MPI_COMM_FESOM ! FESOM communicator (for ocean only runs if often a copy of MPI_COMM_WORLD)
  type(com_struct), pointer     :: com_nod2D
  type(com_struct), pointer     :: com_elem2D
  type(com_struct), pointer     :: com_elem2D_full
  integer                       :: ub, lb ! to work with r(s)_mpitype_elem3D(nod3D)

  integer, dimension(:),     pointer :: s_mpitype_edge2D, r_mpitype_edge2D
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D, r_mpitype_elem2D
  integer, dimension(:),     pointer :: s_mpitype_elem2D_full_i, r_mpitype_elem2D_full_i
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D_full,   r_mpitype_elem2D_full
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D,        r_mpitype_elem3D
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D_full,   r_mpitype_elem3D_full

  integer, dimension(:),     pointer :: s_mpitype_nod2D, r_mpitype_nod2D
  integer, dimension(:),     pointer :: s_mpitype_nod2D_i, r_mpitype_nod2D_i
  integer, dimension(:,:,:), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

  integer, pointer   :: MPIERR
  integer, pointer   :: npes
  integer, pointer   :: mype
  integer, pointer   :: maxPEnum

  integer, dimension(:), pointer     :: part

  ! Mesh partition
  integer, pointer                :: myDim_nod2D, eDim_nod2D
  integer, dimension(:), pointer  :: myList_nod2D
  integer, pointer                :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, dimension(:), pointer  :: myList_elem2D
  integer, pointer                :: myDim_edge2D, eDim_edge2D
  integer, dimension(:), pointer  :: myList_edge2D

  integer, pointer                :: pe_status

  integer, dimension(:), pointer  ::  remPtr_nod2D(:),  remList_nod2D(:)
  integer, dimension(:), pointer  ::  remPtr_elem2D(:), remList_elem2D(:)

  logical, pointer                :: elem_full_flag
# 261 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_def.h" 1
integer         , pointer :: nod2D
integer         , pointer :: elem2D
integer         , pointer :: edge2D
integer         , pointer :: edge2D_in
real(kind=WP)   , pointer :: ocean_area
real(kind=WP)   , pointer :: ocean_areawithcav
integer         , pointer :: nl
integer         , pointer :: nn_size
real(kind=WP), dimension(:,:), pointer :: coord_nod2D, geo_coord_nod2D
integer, dimension(:,:)      , pointer :: elem2D_nodes
integer, dimension(:,:)      , pointer :: edges
integer, dimension(:,:)      , pointer :: edge_tri
integer, dimension(:,:)      , pointer :: elem_edges
real(kind=WP), dimension(:)  , pointer :: elem_area
real(kind=WP), dimension(:,:), pointer :: edge_dxdy, edge_cross_dxdy
real(kind=WP), dimension(:)  , pointer :: elem_cos, metric_factor
integer,       dimension(:,:), pointer :: elem_neighbors
integer,       dimension(:,:), pointer :: nod_in_elem2D
real(kind=WP), dimension(:,:), pointer :: x_corners, y_corners
integer,       dimension(:)  , pointer :: nod_in_elem2D_num
real(kind=WP), dimension(:)  , pointer :: depth
real(kind=WP), dimension(:,:), pointer :: gradient_vec
real(kind=WP), dimension(:,:), pointer :: gradient_sca
integer,       dimension(:)  , pointer :: bc_index_nod2D
real(kind=WP), dimension(:)  , pointer :: zbar, Z, elem_depth
integer,       dimension(:)  , pointer :: nlevels, nlevels_nod2D, nlevels_nod2D_min
real(kind=WP), dimension(:,:), pointer :: area, area_inv, areasvol, areasvol_inv
real(kind=WP), dimension(:)  , pointer :: mesh_resolution
real(kind=WP), dimension(:)  , pointer :: lump2d_north, lump2d_south
type(sparse_matrix)          , pointer :: ssh_stiff
integer,       dimension(:)  , pointer :: cavity_flag_n, cavity_flag_e
real(kind=WP), dimension(:)  , pointer :: cavity_depth
integer,       dimension(:)  , pointer :: ulevels, ulevels_nod2D, ulevels_nod2D_max
integer,       dimension(:)  , pointer :: nn_num
integer,       dimension(:,:), pointer :: nn_pos

real(kind=WP), dimension(:,:), pointer :: hnode
real(kind=WP), dimension(:,:), pointer :: hnode_new
real(kind=WP), dimension(:,:), pointer :: zbar_3d_n
real(kind=WP), dimension(:,:), pointer :: Z_3d_n
real(kind=WP), dimension(:,:), pointer :: helem
real(kind=WP), dimension(:)  , pointer :: bottom_elem_thickness
real(kind=WP), dimension(:)  , pointer :: bottom_node_thickness
real(kind=WP), dimension(:)  , pointer :: dhe
real(kind=WP), dimension(:)  , pointer :: hbar
real(kind=WP), dimension(:)  , pointer :: hbar_old
!real(kind=WP), dimension(:)  , pointer :: zbar_n
!real(kind=WP), dimension(:)  , pointer :: Z_n
real(kind=WP), dimension(:)  , pointer :: zbar_n_bot
real(kind=WP), dimension(:)  , pointer :: zbar_e_bot
real(kind=WP), dimension(:)  , pointer :: zbar_n_srf
real(kind=WP), dimension(:)  , pointer :: zbar_e_srf
# 262 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_ass.h" 1
MPI_COMM_FESOM  => partit%MPI_COMM_FESOM
com_nod2D       => partit%com_nod2D
com_elem2D      => partit%com_elem2D
com_elem2D_full => partit%com_elem2D_full
myDim_nod2D     => partit%myDim_nod2D
eDim_nod2D      => partit%eDim_nod2D
myDim_elem2D    => partit%myDim_elem2D
eDim_elem2D     => partit%eDim_elem2D
eXDim_elem2D    => partit%eXDim_elem2D
myDim_edge2D    => partit%myDim_edge2D
eDim_edge2D     => partit%eDim_edge2D
pe_status       => partit%pe_status
elem_full_flag  => partit%elem_full_flag
MPIERR          => partit%MPIERR
npes            => partit%npes
mype            => partit%mype
maxPEnum        => partit%maxPEnum
part            => partit%part

lb=lbound(partit%s_mpitype_elem3D, 2)
ub=ubound(partit%s_mpitype_elem3D, 2)

myList_nod2D (1:myDim_nod2D +eDim_nod2D)                => partit%myList_nod2D(:)
myList_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)  => partit%myList_elem2D(:)
myList_edge2D(1:myDim_edge2D+eDim_edge2D)               => partit%myList_edge2D(:)

if (allocated(partit%remPtr_nod2D)) then
   remPtr_nod2D  (1:npes)                => partit%remPtr_nod2D(:)
   remList_nod2D (1:remPtr_nod2D(npes))  => partit%remList_nod2D(:)
end if

if (allocated(partit%remPtr_elem2D)) then
remPtr_elem2D (1:npes)                => partit%remPtr_elem2D(:)
remList_elem2D(1:remPtr_elem2D(npes)) => partit%remList_elem2D(:)
end if

s_mpitype_elem2D(1:com_elem2D%sPEnum, 1:4) => partit%s_mpitype_elem2D(:,:)
r_mpitype_elem2D(1:com_elem2D%rPEnum, 1:4) => partit%r_mpitype_elem2D(:,:)

s_mpitype_elem2D_full_i(1:com_elem2D_full%sPEnum) => partit%s_mpitype_elem2D_full_i(:)
r_mpitype_elem2D_full_i(1:com_elem2D_full%rPEnum) => partit%r_mpitype_elem2D_full_i(:)

s_mpitype_elem2D_full(1:com_elem2D_full%sPEnum, 1:4) => partit%s_mpitype_elem2D_full(:,:)
r_mpitype_elem2D_full(1:com_elem2D_full%rPEnum, 1:4) => partit%r_mpitype_elem2D_full(:,:)

s_mpitype_elem3D(1:com_elem2D%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D(:,:,:)
r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D(:,:,:)

s_mpitype_elem3D_full(1:com_elem2D_full%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D_full(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4)           => partit%r_mpitype_elem3D(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

s_mpitype_nod2D(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D(:)
r_mpitype_nod2D(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D(:)

s_mpitype_nod2D_i(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D_i(:)
r_mpitype_nod2D_i(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D_i(:)

s_mpitype_nod3D(1:com_nod2D%sPEnum, lb:ub, 1:3) => partit%s_mpitype_nod3D(:,:,:)
r_mpitype_nod3D(1:com_nod2D%rPEnum, lb:ub, 1:3) => partit%r_mpitype_nod3D(:,:,:)
# 263 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_ass.h" 1
nod2D              => mesh%nod2D
elem2D             => mesh%elem2D
edge2D             => mesh%edge2D
edge2D_in          => mesh%edge2D_in
ocean_area         => mesh%ocean_area
nl                 => mesh%nl
nn_size            => mesh%nn_size
ocean_areawithcav  => mesh%ocean_areawithcav
coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)                  => mesh%coord_nod2D(:,:)
geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)              => mesh%geo_coord_nod2D(:,:)
elem2D_nodes(1:3, 1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%elem2D_nodes(:,:)
edges(1:2,1:myDim_edge2D+eDim_edge2D)                      => mesh%edges(:,:)
edge_tri(1:2,1:myDim_edge2D+eDim_edge2D)                   => mesh%edge_tri(:,:)
elem_edges(1:3,1:myDim_elem2D)                             => mesh%elem_edges(:,:)
elem_area(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)         => mesh%elem_area(:)
edge_dxdy(1:2,1:myDim_edge2D+eDim_edge2D)                  => mesh%edge_dxdy(:,:)
edge_cross_dxdy(1:4,1:myDim_edge2D+eDim_edge2D)            => mesh%edge_cross_dxdy(:,:)
elem_cos(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)          => mesh%elem_cos(:)
metric_factor(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%metric_factor(:)
elem_neighbors(1:3,1:myDim_elem2D)                         => mesh%elem_neighbors(:,:)
nod_in_elem2D      => mesh%nod_in_elem2D   ! (maxval(rmax),myDim_nod2D+eDim_nod2D)
x_corners          => mesh%x_corners   ! (myDim_nod2D, maxval(rmax))
y_corners          => mesh%y_corners   ! (myDim_nod2D, maxval(rmax))
nod_in_elem2D_num(1:myDim_nod2D+eDim_nod2D)                => mesh%nod_in_elem2D_num(:)
depth(1:myDim_nod2D+eDim_nod2D)                            => mesh%depth(:)
gradient_vec(1:6,1:myDim_elem2D)                           => mesh%gradient_vec(:,:)
gradient_sca(1:6,1:myDim_elem2D)                           => mesh%gradient_sca(:,:)
bc_index_nod2D(1:myDim_nod2D+eDim_nod2D)                   => mesh%bc_index_nod2D(:)
zbar(1:mesh%nl)                                            => mesh%zbar(:)
Z(1:mesh%nl-1)                                             => mesh%Z(:)
elem_depth         => mesh%elem_depth      ! never used, not even allocated
nlevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%nlevels(:)
nlevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%nlevels_nod2D(:)
nlevels_nod2D_min(1:myDim_nod2D+eDim_nod2D)                => mesh%nlevels_nod2D_min(:)
area(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)                   => mesh%area(:,:)
areasvol(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%areasvol(:,:)
area_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%area_inv(:,:)
areasvol_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)           => mesh%areasvol_inv(:,:)
mesh_resolution(1:myDim_nod2d+eDim_nod2D)                  => mesh%mesh_resolution(:)
ssh_stiff                                                  => mesh%ssh_stiff
lump2d_north(1:myDim_nod2d)                                => mesh%lump2d_north(:)
lump2d_south(1:myDim_nod2d)                                => mesh%lump2d_south(:)
cavity_flag_n(1:myDim_nod2D+eDim_nod2D)                    => mesh%cavity_flag_n(:)
cavity_flag_e(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%cavity_flag_e(:)
!!$cavity_lev_nod2D(1:myDim_nod2D+eDim_nod2D)                 => mesh%cavity_lev_nod2D
!!$cavity_lev_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%cavity_lev_elem2D
cavity_depth(1:myDim_nod2D+eDim_nod2D)                     => mesh%cavity_depth(:)
ulevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%ulevels(:)
ulevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%ulevels_nod2D(:)
ulevels_nod2D_max(1:myDim_nod2D+eDim_nod2D)                => mesh%ulevels_nod2D_max(:)
nn_num(1:myDim_nod2D)                                      => mesh%nn_num(:)
nn_pos(1:mesh%nn_size, 1:myDim_nod2D)                      => mesh%nn_pos(:,:)
hnode(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)               => mesh%hnode(:,:)
hnode_new(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)           => mesh%hnode_new(:,:)
zbar_3d_n(1:mesh%nl, 1:myDim_nod2D+eDim_nod2D)             => mesh%zbar_3d_n(:,:)
Z_3d_n(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)              => mesh%Z_3d_n(:,:)
helem(1:mesh%nl-1, 1:myDim_elem2D)                         => mesh%helem(:,:)
bottom_elem_thickness(1:myDim_elem2D)                      => mesh%bottom_elem_thickness(:)
bottom_node_thickness(1:myDim_nod2D+eDim_nod2D)            => mesh%bottom_node_thickness(:)
dhe(1:myDim_elem2D)                                        => mesh%dhe(:)
hbar(1:myDim_nod2D+eDim_nod2D)                             => mesh%hbar(:)
hbar_old(1:myDim_nod2D+eDim_nod2D)                         => mesh%hbar_old(:)
!zbar_n(1:mesh%nl)                                          => mesh%zbar_n
!Z_n(1:mesh%nl-1)                                           => mesh%Z_n
zbar_n_bot(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_bot(:)
zbar_e_bot(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_bot(:)
zbar_n_srf(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_srf(:)
zbar_e_srf(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_srf(:)
# 264 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

    l_init_zero=.true.
    if (present(o_init_zero)) then
       l_init_zero=o_init_zero
    end if
    if (l_init_zero) then
!$OMP PARALLEL DO
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
        do n=1, myDim_nod2D
          do nz=1,mesh%nl
            flux(nz, n)=0.0_WP
          end do
        end do
        !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO
    end if
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tvert, n, nz, nzmax, nzmin)
!$OMP DO

    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1, myDim_nod2D
       !_______________________________________________________________________
       nzmax=nlevels_nod2D(n)
       nzmin=ulevels_nod2D(n)

       !_______________________________________________________________________
       ! vert. flux at surface layer
       nz=nzmin
       flux(nz,n)=-W(nz,n)*ttf(nz,n)*area(nz,n)-flux(nz,n)

       !_______________________________________________________________________
       ! vert. flux at bottom layer --> zero bottom flux
       nz=nzmax
       flux(nz,n)= 0.0_WP-flux(nz,n)

       !_______________________________________________________________________
       ! Be carefull have to do vertical tracer advection here on old vertical grid
       ! also horizontal advection is done on old mesh (see helem contains old
       ! mesh information)
       !_______________________________________________________________________
       ! vert. flux at remaining levels
       !$ACC LOOP VECTOR
       do nz=nzmin+1,nzmax-1
          flux(nz,n)=-0.5*(                                                        &
                      ttf(nz  ,n)*(W(nz,n)+abs(W(nz,n)))+ &
                      ttf(nz-1,n)*(W(nz,n)-abs(W(nz,n))))*area(nz,n)-flux(nz,n)
       end do
       !$ACC END LOOP
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO
!$OMP END PARALLEL
end subroutine adv_tra_ver_upw1
!
!
!===============================================================================
subroutine adv_tra_ver_qr4c(w, ttf, partit, mesh, num_ord, flux, o_init_zero)
    use MOD_MESH
    use o_ARRAYS
    use o_PARAM
    USE MOD_PARTIT
    USE MOD_PARSUP
    implicit none
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(in)    :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
    real(kind=WP), intent(in)    :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)    :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout) :: flux(mesh%nl,  partit%myDim_nod2D)
    logical, optional            :: o_init_zero
    logical                      :: l_init_zero
    real(kind=WP)                :: tvert(mesh%nl)
    integer                      :: node, vertical_level, nzmax, nzmin
    real(kind=WP)                :: Tmean, Tmean1, Tmean2
    real(kind=WP)                :: qc, qu, qd


# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_def.h" 1

  integer,          pointer     :: MPI_COMM_FESOM ! FESOM communicator (for ocean only runs if often a copy of MPI_COMM_WORLD)
  type(com_struct), pointer     :: com_nod2D
  type(com_struct), pointer     :: com_elem2D
  type(com_struct), pointer     :: com_elem2D_full
  integer                       :: ub, lb ! to work with r(s)_mpitype_elem3D(nod3D)

  integer, dimension(:),     pointer :: s_mpitype_edge2D, r_mpitype_edge2D
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D, r_mpitype_elem2D
  integer, dimension(:),     pointer :: s_mpitype_elem2D_full_i, r_mpitype_elem2D_full_i
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D_full,   r_mpitype_elem2D_full
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D,        r_mpitype_elem3D
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D_full,   r_mpitype_elem3D_full

  integer, dimension(:),     pointer :: s_mpitype_nod2D, r_mpitype_nod2D
  integer, dimension(:),     pointer :: s_mpitype_nod2D_i, r_mpitype_nod2D_i
  integer, dimension(:,:,:), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

  integer, pointer   :: MPIERR
  integer, pointer   :: npes
  integer, pointer   :: mype
  integer, pointer   :: maxPEnum

  integer, dimension(:), pointer     :: part

  ! Mesh partition
  integer, pointer                :: myDim_nod2D, eDim_nod2D
  integer, dimension(:), pointer  :: myList_nod2D
  integer, pointer                :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, dimension(:), pointer  :: myList_elem2D
  integer, pointer                :: myDim_edge2D, eDim_edge2D
  integer, dimension(:), pointer  :: myList_edge2D

  integer, pointer                :: pe_status

  integer, dimension(:), pointer  ::  remPtr_nod2D(:),  remList_nod2D(:)
  integer, dimension(:), pointer  ::  remPtr_elem2D(:), remList_elem2D(:)

  logical, pointer                :: elem_full_flag
# 341 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_def.h" 1
integer         , pointer :: nod2D
integer         , pointer :: elem2D
integer         , pointer :: edge2D
integer         , pointer :: edge2D_in
real(kind=WP)   , pointer :: ocean_area
real(kind=WP)   , pointer :: ocean_areawithcav
integer         , pointer :: nl
integer         , pointer :: nn_size
real(kind=WP), dimension(:,:), pointer :: coord_nod2D, geo_coord_nod2D
integer, dimension(:,:)      , pointer :: elem2D_nodes
integer, dimension(:,:)      , pointer :: edges
integer, dimension(:,:)      , pointer :: edge_tri
integer, dimension(:,:)      , pointer :: elem_edges
real(kind=WP), dimension(:)  , pointer :: elem_area
real(kind=WP), dimension(:,:), pointer :: edge_dxdy, edge_cross_dxdy
real(kind=WP), dimension(:)  , pointer :: elem_cos, metric_factor
integer,       dimension(:,:), pointer :: elem_neighbors
integer,       dimension(:,:), pointer :: nod_in_elem2D
real(kind=WP), dimension(:,:), pointer :: x_corners, y_corners
integer,       dimension(:)  , pointer :: nod_in_elem2D_num
real(kind=WP), dimension(:)  , pointer :: depth
real(kind=WP), dimension(:,:), pointer :: gradient_vec
real(kind=WP), dimension(:,:), pointer :: gradient_sca
integer,       dimension(:)  , pointer :: bc_index_nod2D
real(kind=WP), dimension(:)  , pointer :: zbar, Z, elem_depth
integer,       dimension(:)  , pointer :: nlevels, nlevels_nod2D, nlevels_nod2D_min
real(kind=WP), dimension(:,:), pointer :: area, area_inv, areasvol, areasvol_inv
real(kind=WP), dimension(:)  , pointer :: mesh_resolution
real(kind=WP), dimension(:)  , pointer :: lump2d_north, lump2d_south
type(sparse_matrix)          , pointer :: ssh_stiff
integer,       dimension(:)  , pointer :: cavity_flag_n, cavity_flag_e
real(kind=WP), dimension(:)  , pointer :: cavity_depth
integer,       dimension(:)  , pointer :: ulevels, ulevels_nod2D, ulevels_nod2D_max
integer,       dimension(:)  , pointer :: nn_num
integer,       dimension(:,:), pointer :: nn_pos

real(kind=WP), dimension(:,:), pointer :: hnode
real(kind=WP), dimension(:,:), pointer :: hnode_new
real(kind=WP), dimension(:,:), pointer :: zbar_3d_n
real(kind=WP), dimension(:,:), pointer :: Z_3d_n
real(kind=WP), dimension(:,:), pointer :: helem
real(kind=WP), dimension(:)  , pointer :: bottom_elem_thickness
real(kind=WP), dimension(:)  , pointer :: bottom_node_thickness
real(kind=WP), dimension(:)  , pointer :: dhe
real(kind=WP), dimension(:)  , pointer :: hbar
real(kind=WP), dimension(:)  , pointer :: hbar_old
!real(kind=WP), dimension(:)  , pointer :: zbar_n
!real(kind=WP), dimension(:)  , pointer :: Z_n
real(kind=WP), dimension(:)  , pointer :: zbar_n_bot
real(kind=WP), dimension(:)  , pointer :: zbar_e_bot
real(kind=WP), dimension(:)  , pointer :: zbar_n_srf
real(kind=WP), dimension(:)  , pointer :: zbar_e_srf
# 342 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_ass.h" 1
MPI_COMM_FESOM  => partit%MPI_COMM_FESOM
com_nod2D       => partit%com_nod2D
com_elem2D      => partit%com_elem2D
com_elem2D_full => partit%com_elem2D_full
myDim_nod2D     => partit%myDim_nod2D
eDim_nod2D      => partit%eDim_nod2D
myDim_elem2D    => partit%myDim_elem2D
eDim_elem2D     => partit%eDim_elem2D
eXDim_elem2D    => partit%eXDim_elem2D
myDim_edge2D    => partit%myDim_edge2D
eDim_edge2D     => partit%eDim_edge2D
pe_status       => partit%pe_status
elem_full_flag  => partit%elem_full_flag
MPIERR          => partit%MPIERR
npes            => partit%npes
mype            => partit%mype
maxPEnum        => partit%maxPEnum
part            => partit%part

lb=lbound(partit%s_mpitype_elem3D, 2)
ub=ubound(partit%s_mpitype_elem3D, 2)

myList_nod2D (1:myDim_nod2D +eDim_nod2D)                => partit%myList_nod2D(:)
myList_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)  => partit%myList_elem2D(:)
myList_edge2D(1:myDim_edge2D+eDim_edge2D)               => partit%myList_edge2D(:)

if (allocated(partit%remPtr_nod2D)) then
   remPtr_nod2D  (1:npes)                => partit%remPtr_nod2D(:)
   remList_nod2D (1:remPtr_nod2D(npes))  => partit%remList_nod2D(:)
end if

if (allocated(partit%remPtr_elem2D)) then
remPtr_elem2D (1:npes)                => partit%remPtr_elem2D(:)
remList_elem2D(1:remPtr_elem2D(npes)) => partit%remList_elem2D(:)
end if

s_mpitype_elem2D(1:com_elem2D%sPEnum, 1:4) => partit%s_mpitype_elem2D(:,:)
r_mpitype_elem2D(1:com_elem2D%rPEnum, 1:4) => partit%r_mpitype_elem2D(:,:)

s_mpitype_elem2D_full_i(1:com_elem2D_full%sPEnum) => partit%s_mpitype_elem2D_full_i(:)
r_mpitype_elem2D_full_i(1:com_elem2D_full%rPEnum) => partit%r_mpitype_elem2D_full_i(:)

s_mpitype_elem2D_full(1:com_elem2D_full%sPEnum, 1:4) => partit%s_mpitype_elem2D_full(:,:)
r_mpitype_elem2D_full(1:com_elem2D_full%rPEnum, 1:4) => partit%r_mpitype_elem2D_full(:,:)

s_mpitype_elem3D(1:com_elem2D%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D(:,:,:)
r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D(:,:,:)

s_mpitype_elem3D_full(1:com_elem2D_full%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D_full(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4)           => partit%r_mpitype_elem3D(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

s_mpitype_nod2D(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D(:)
r_mpitype_nod2D(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D(:)

s_mpitype_nod2D_i(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D_i(:)
r_mpitype_nod2D_i(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D_i(:)

s_mpitype_nod3D(1:com_nod2D%sPEnum, lb:ub, 1:3) => partit%s_mpitype_nod3D(:,:,:)
r_mpitype_nod3D(1:com_nod2D%rPEnum, lb:ub, 1:3) => partit%r_mpitype_nod3D(:,:,:)
# 343 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_ass.h" 1
nod2D              => mesh%nod2D
elem2D             => mesh%elem2D
edge2D             => mesh%edge2D
edge2D_in          => mesh%edge2D_in
ocean_area         => mesh%ocean_area
nl                 => mesh%nl
nn_size            => mesh%nn_size
ocean_areawithcav  => mesh%ocean_areawithcav
coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)                  => mesh%coord_nod2D(:,:)
geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)              => mesh%geo_coord_nod2D(:,:)
elem2D_nodes(1:3, 1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%elem2D_nodes(:,:)
edges(1:2,1:myDim_edge2D+eDim_edge2D)                      => mesh%edges(:,:)
edge_tri(1:2,1:myDim_edge2D+eDim_edge2D)                   => mesh%edge_tri(:,:)
elem_edges(1:3,1:myDim_elem2D)                             => mesh%elem_edges(:,:)
elem_area(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)         => mesh%elem_area(:)
edge_dxdy(1:2,1:myDim_edge2D+eDim_edge2D)                  => mesh%edge_dxdy(:,:)
edge_cross_dxdy(1:4,1:myDim_edge2D+eDim_edge2D)            => mesh%edge_cross_dxdy(:,:)
elem_cos(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)          => mesh%elem_cos(:)
metric_factor(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%metric_factor(:)
elem_neighbors(1:3,1:myDim_elem2D)                         => mesh%elem_neighbors(:,:)
nod_in_elem2D      => mesh%nod_in_elem2D   ! (maxval(rmax),myDim_nod2D+eDim_nod2D)
x_corners          => mesh%x_corners   ! (myDim_nod2D, maxval(rmax))
y_corners          => mesh%y_corners   ! (myDim_nod2D, maxval(rmax))
nod_in_elem2D_num(1:myDim_nod2D+eDim_nod2D)                => mesh%nod_in_elem2D_num(:)
depth(1:myDim_nod2D+eDim_nod2D)                            => mesh%depth(:)
gradient_vec(1:6,1:myDim_elem2D)                           => mesh%gradient_vec(:,:)
gradient_sca(1:6,1:myDim_elem2D)                           => mesh%gradient_sca(:,:)
bc_index_nod2D(1:myDim_nod2D+eDim_nod2D)                   => mesh%bc_index_nod2D(:)
zbar(1:mesh%nl)                                            => mesh%zbar(:)
Z(1:mesh%nl-1)                                             => mesh%Z(:)
elem_depth         => mesh%elem_depth      ! never used, not even allocated
nlevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%nlevels(:)
nlevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%nlevels_nod2D(:)
nlevels_nod2D_min(1:myDim_nod2D+eDim_nod2D)                => mesh%nlevels_nod2D_min(:)
area(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)                   => mesh%area(:,:)
areasvol(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%areasvol(:,:)
area_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%area_inv(:,:)
areasvol_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)           => mesh%areasvol_inv(:,:)
mesh_resolution(1:myDim_nod2d+eDim_nod2D)                  => mesh%mesh_resolution(:)
ssh_stiff                                                  => mesh%ssh_stiff
lump2d_north(1:myDim_nod2d)                                => mesh%lump2d_north(:)
lump2d_south(1:myDim_nod2d)                                => mesh%lump2d_south(:)
cavity_flag_n(1:myDim_nod2D+eDim_nod2D)                    => mesh%cavity_flag_n(:)
cavity_flag_e(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%cavity_flag_e(:)
!!$cavity_lev_nod2D(1:myDim_nod2D+eDim_nod2D)                 => mesh%cavity_lev_nod2D
!!$cavity_lev_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%cavity_lev_elem2D
cavity_depth(1:myDim_nod2D+eDim_nod2D)                     => mesh%cavity_depth(:)
ulevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%ulevels(:)
ulevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%ulevels_nod2D(:)
ulevels_nod2D_max(1:myDim_nod2D+eDim_nod2D)                => mesh%ulevels_nod2D_max(:)
nn_num(1:myDim_nod2D)                                      => mesh%nn_num(:)
nn_pos(1:mesh%nn_size, 1:myDim_nod2D)                      => mesh%nn_pos(:,:)
hnode(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)               => mesh%hnode(:,:)
hnode_new(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)           => mesh%hnode_new(:,:)
zbar_3d_n(1:mesh%nl, 1:myDim_nod2D+eDim_nod2D)             => mesh%zbar_3d_n(:,:)
Z_3d_n(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)              => mesh%Z_3d_n(:,:)
helem(1:mesh%nl-1, 1:myDim_elem2D)                         => mesh%helem(:,:)
bottom_elem_thickness(1:myDim_elem2D)                      => mesh%bottom_elem_thickness(:)
bottom_node_thickness(1:myDim_nod2D+eDim_nod2D)            => mesh%bottom_node_thickness(:)
dhe(1:myDim_elem2D)                                        => mesh%dhe(:)
hbar(1:myDim_nod2D+eDim_nod2D)                             => mesh%hbar(:)
hbar_old(1:myDim_nod2D+eDim_nod2D)                         => mesh%hbar_old(:)
!zbar_n(1:mesh%nl)                                          => mesh%zbar_n
!Z_n(1:mesh%nl-1)                                           => mesh%Z_n
zbar_n_bot(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_bot(:)
zbar_e_bot(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_bot(:)
zbar_n_srf(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_srf(:)
zbar_e_srf(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_srf(:)
# 344 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/edsl/interface.h" 1







# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/edsl/cpu_backend.h" 1









# 8 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/edsl/interface.h" 2

# 345 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

    l_init_zero=.true.
    if (present(o_init_zero)) then
       l_init_zero=o_init_zero
    end if
    if (l_init_zero) then
!$OMP PARALLEL DO
       !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
       do node = 1, myDim_nod2D
          do vertical_level = 1, mesh%nl
             flux(vertical_level, node) = 0.0_WP
          end do
       end do
       !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO
    end if
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tvert, node, vertical_level, nzmax, nzmin, Tmean, Tmean1, Tmean2, qc, qu,qd)
!$OMP DO
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do node =  1,  myDim_nod2D
        !_______________________________________________________________________
        nzmax=nlevels_nod2D(node)
        nzmin=ulevels_nod2D(node)

        !_______________________________________________________________________
        ! vert. flux at surface layer
        vertical_level =  nzmin
        flux(vertical_level, node) = - ttf (vertical_level, node) &
                                     * W   (vertical_level, node) &
                                     * area(vertical_level, node) &
                                     - flux(vertical_level, node)

        !_______________________________________________________________________
        ! vert. flux 2nd layer --> centered differences
        vertical_level =  nzmin + 1
        flux(vertical_level, node) = - 0.5_WP                            &
                                     * ( ttf(vertical_level - 1, node)   &
                                       + ttf(vertical_level,     node) ) &
                                     * W    (vertical_level,     node)   &
                                     * area (vertical_level,     node)   &
                                     - flux (vertical_level,     node)

        !_______________________________________________________________________
        ! vert. flux at bottom - 1 layer --> centered differences
        vertical_level =  nzmax - 1
        flux(vertical_level, node) = - 0.5_WP                            &
                                     * ( ttf(vertical_level - 1, node)   &
                                       + ttf(vertical_level,     node) ) &
                                     * W    (vertical_level,     node)   &
                                     * area (vertical_level,     node)   &
                                     - flux (vertical_level,     node)

        !_______________________________________________________________________
        ! vert. flux at bottom layer --> zero bottom flux
        vertical_level =  nzmax
        flux(vertical_level, node) = 0.0_WP - flux(vertical_level, node)


        !_______________________________________________________________________
        ! Be carefull have to do vertical tracer advection here on old vertical grid
        ! also horizontal advection is done on old mesh (see helem contains old
        ! mesh information)
        !_______________________________________________________________________
        ! vert. flux at remaining levels
        do vertical_level =  nzmin + 2,  nzmax - 2
            !centered (4th order)
            qc = ( ttf   (vertical_level - 1, node) - ttf   (vertical_level,     node) ) &
               / ( Z_3d_n(vertical_level - 1, node) - Z_3d_n(vertical_level,     node) )
            qu = ( ttf   (vertical_level,     node) - ttf   (vertical_level + 1, node) ) &
               / ( Z_3d_n(vertical_level,     node) - Z_3d_n(vertical_level + 1, node) )
            qd = ( ttf   (vertical_level - 2, node) - ttf   (vertical_level - 1, node) ) &
               / ( Z_3d_n(vertical_level - 2, node) - Z_3d_n(vertical_level - 1, node) )

            Tmean1 = ttf(vertical_level    , node) + (2 * qc + qu) &
                   * (zbar_3d_n(vertical_level,node) - Z_3d_n(vertical_level  , node)) &
                   / 3.0_WP
            Tmean2 = ttf(vertical_level - 1, node) + (2 * qc + qd) &
                   * (zbar_3d_n(vertical_level,node) - Z_3d_n(vertical_level - 1, node)) &
                   / 3.0_WP
            Tmean = ( W(vertical_level,node) + abs(W(vertical_level, node)) ) * Tmean1 &
                  + ( W(vertical_level,node) - abs(W(vertical_level, node)) ) * Tmean2

            flux(vertical_level, node) = (                                                          &
                                           -  0.5_WP * (1.0_WP - num_ord) * Tmean - num_ord         &
                                           * (0.5_WP * (Tmean1 + Tmean2)) * W(vertical_level, node) &
                                         )                                                          &
                                       * area(vertical_level, node)                                 &
                                       - flux(vertical_level, node)

        end do
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO
!$OMP END PARALLEL
end subroutine adv_tra_ver_qr4c
!
!
!===============================================================================
subroutine adv_tra_vert_ppm(dt, w, ttf, partit, mesh, flux, o_init_zero)
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto
    implicit none
    real(kind=WP), intent(in),  target :: dt
    type(t_partit),intent(in), target  :: partit
    type(t_mesh),  intent(in) , target :: mesh
    real(kind=WP), intent(in)          :: ttf (mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)          :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)       :: flux(mesh%nl,   partit%myDim_nod2D)
    logical, optional                  :: o_init_zero
    logical                            :: l_init_zero
    real(kind=WP)                      :: tvert(mesh%nl), tv(mesh%nl), aL, aR, aj, x
    real(kind=WP)                      :: dzjm1, dzj, dzjp1, dzjp2, deltaj, deltajp1
    integer                            :: n, nz, nzmax, nzmin
!   integer                            :: overshoot_counter, counter


# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_def.h" 1

  integer,          pointer     :: MPI_COMM_FESOM ! FESOM communicator (for ocean only runs if often a copy of MPI_COMM_WORLD)
  type(com_struct), pointer     :: com_nod2D
  type(com_struct), pointer     :: com_elem2D
  type(com_struct), pointer     :: com_elem2D_full
  integer                       :: ub, lb ! to work with r(s)_mpitype_elem3D(nod3D)

  integer, dimension(:),     pointer :: s_mpitype_edge2D, r_mpitype_edge2D
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D, r_mpitype_elem2D
  integer, dimension(:),     pointer :: s_mpitype_elem2D_full_i, r_mpitype_elem2D_full_i
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D_full,   r_mpitype_elem2D_full
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D,        r_mpitype_elem3D
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D_full,   r_mpitype_elem3D_full

  integer, dimension(:),     pointer :: s_mpitype_nod2D, r_mpitype_nod2D
  integer, dimension(:),     pointer :: s_mpitype_nod2D_i, r_mpitype_nod2D_i
  integer, dimension(:,:,:), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

  integer, pointer   :: MPIERR
  integer, pointer   :: npes
  integer, pointer   :: mype
  integer, pointer   :: maxPEnum

  integer, dimension(:), pointer     :: part

  ! Mesh partition
  integer, pointer                :: myDim_nod2D, eDim_nod2D
  integer, dimension(:), pointer  :: myList_nod2D
  integer, pointer                :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, dimension(:), pointer  :: myList_elem2D
  integer, pointer                :: myDim_edge2D, eDim_edge2D
  integer, dimension(:), pointer  :: myList_edge2D

  integer, pointer                :: pe_status

  integer, dimension(:), pointer  ::  remPtr_nod2D(:),  remList_nod2D(:)
  integer, dimension(:), pointer  ::  remPtr_elem2D(:), remList_elem2D(:)

  logical, pointer                :: elem_full_flag
# 467 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_def.h" 1
integer         , pointer :: nod2D
integer         , pointer :: elem2D
integer         , pointer :: edge2D
integer         , pointer :: edge2D_in
real(kind=WP)   , pointer :: ocean_area
real(kind=WP)   , pointer :: ocean_areawithcav
integer         , pointer :: nl
integer         , pointer :: nn_size
real(kind=WP), dimension(:,:), pointer :: coord_nod2D, geo_coord_nod2D
integer, dimension(:,:)      , pointer :: elem2D_nodes
integer, dimension(:,:)      , pointer :: edges
integer, dimension(:,:)      , pointer :: edge_tri
integer, dimension(:,:)      , pointer :: elem_edges
real(kind=WP), dimension(:)  , pointer :: elem_area
real(kind=WP), dimension(:,:), pointer :: edge_dxdy, edge_cross_dxdy
real(kind=WP), dimension(:)  , pointer :: elem_cos, metric_factor
integer,       dimension(:,:), pointer :: elem_neighbors
integer,       dimension(:,:), pointer :: nod_in_elem2D
real(kind=WP), dimension(:,:), pointer :: x_corners, y_corners
integer,       dimension(:)  , pointer :: nod_in_elem2D_num
real(kind=WP), dimension(:)  , pointer :: depth
real(kind=WP), dimension(:,:), pointer :: gradient_vec
real(kind=WP), dimension(:,:), pointer :: gradient_sca
integer,       dimension(:)  , pointer :: bc_index_nod2D
real(kind=WP), dimension(:)  , pointer :: zbar, Z, elem_depth
integer,       dimension(:)  , pointer :: nlevels, nlevels_nod2D, nlevels_nod2D_min
real(kind=WP), dimension(:,:), pointer :: area, area_inv, areasvol, areasvol_inv
real(kind=WP), dimension(:)  , pointer :: mesh_resolution
real(kind=WP), dimension(:)  , pointer :: lump2d_north, lump2d_south
type(sparse_matrix)          , pointer :: ssh_stiff
integer,       dimension(:)  , pointer :: cavity_flag_n, cavity_flag_e
real(kind=WP), dimension(:)  , pointer :: cavity_depth
integer,       dimension(:)  , pointer :: ulevels, ulevels_nod2D, ulevels_nod2D_max
integer,       dimension(:)  , pointer :: nn_num
integer,       dimension(:,:), pointer :: nn_pos

real(kind=WP), dimension(:,:), pointer :: hnode
real(kind=WP), dimension(:,:), pointer :: hnode_new
real(kind=WP), dimension(:,:), pointer :: zbar_3d_n
real(kind=WP), dimension(:,:), pointer :: Z_3d_n
real(kind=WP), dimension(:,:), pointer :: helem
real(kind=WP), dimension(:)  , pointer :: bottom_elem_thickness
real(kind=WP), dimension(:)  , pointer :: bottom_node_thickness
real(kind=WP), dimension(:)  , pointer :: dhe
real(kind=WP), dimension(:)  , pointer :: hbar
real(kind=WP), dimension(:)  , pointer :: hbar_old
!real(kind=WP), dimension(:)  , pointer :: zbar_n
!real(kind=WP), dimension(:)  , pointer :: Z_n
real(kind=WP), dimension(:)  , pointer :: zbar_n_bot
real(kind=WP), dimension(:)  , pointer :: zbar_e_bot
real(kind=WP), dimension(:)  , pointer :: zbar_n_srf
real(kind=WP), dimension(:)  , pointer :: zbar_e_srf
# 468 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_ass.h" 1
MPI_COMM_FESOM  => partit%MPI_COMM_FESOM
com_nod2D       => partit%com_nod2D
com_elem2D      => partit%com_elem2D
com_elem2D_full => partit%com_elem2D_full
myDim_nod2D     => partit%myDim_nod2D
eDim_nod2D      => partit%eDim_nod2D
myDim_elem2D    => partit%myDim_elem2D
eDim_elem2D     => partit%eDim_elem2D
eXDim_elem2D    => partit%eXDim_elem2D
myDim_edge2D    => partit%myDim_edge2D
eDim_edge2D     => partit%eDim_edge2D
pe_status       => partit%pe_status
elem_full_flag  => partit%elem_full_flag
MPIERR          => partit%MPIERR
npes            => partit%npes
mype            => partit%mype
maxPEnum        => partit%maxPEnum
part            => partit%part

lb=lbound(partit%s_mpitype_elem3D, 2)
ub=ubound(partit%s_mpitype_elem3D, 2)

myList_nod2D (1:myDim_nod2D +eDim_nod2D)                => partit%myList_nod2D(:)
myList_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)  => partit%myList_elem2D(:)
myList_edge2D(1:myDim_edge2D+eDim_edge2D)               => partit%myList_edge2D(:)

if (allocated(partit%remPtr_nod2D)) then
   remPtr_nod2D  (1:npes)                => partit%remPtr_nod2D(:)
   remList_nod2D (1:remPtr_nod2D(npes))  => partit%remList_nod2D(:)
end if

if (allocated(partit%remPtr_elem2D)) then
remPtr_elem2D (1:npes)                => partit%remPtr_elem2D(:)
remList_elem2D(1:remPtr_elem2D(npes)) => partit%remList_elem2D(:)
end if

s_mpitype_elem2D(1:com_elem2D%sPEnum, 1:4) => partit%s_mpitype_elem2D(:,:)
r_mpitype_elem2D(1:com_elem2D%rPEnum, 1:4) => partit%r_mpitype_elem2D(:,:)

s_mpitype_elem2D_full_i(1:com_elem2D_full%sPEnum) => partit%s_mpitype_elem2D_full_i(:)
r_mpitype_elem2D_full_i(1:com_elem2D_full%rPEnum) => partit%r_mpitype_elem2D_full_i(:)

s_mpitype_elem2D_full(1:com_elem2D_full%sPEnum, 1:4) => partit%s_mpitype_elem2D_full(:,:)
r_mpitype_elem2D_full(1:com_elem2D_full%rPEnum, 1:4) => partit%r_mpitype_elem2D_full(:,:)

s_mpitype_elem3D(1:com_elem2D%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D(:,:,:)
r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D(:,:,:)

s_mpitype_elem3D_full(1:com_elem2D_full%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D_full(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4)           => partit%r_mpitype_elem3D(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

s_mpitype_nod2D(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D(:)
r_mpitype_nod2D(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D(:)

s_mpitype_nod2D_i(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D_i(:)
r_mpitype_nod2D_i(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D_i(:)

s_mpitype_nod3D(1:com_nod2D%sPEnum, lb:ub, 1:3) => partit%s_mpitype_nod3D(:,:,:)
r_mpitype_nod3D(1:com_nod2D%rPEnum, lb:ub, 1:3) => partit%r_mpitype_nod3D(:,:,:)
# 469 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_ass.h" 1
nod2D              => mesh%nod2D
elem2D             => mesh%elem2D
edge2D             => mesh%edge2D
edge2D_in          => mesh%edge2D_in
ocean_area         => mesh%ocean_area
nl                 => mesh%nl
nn_size            => mesh%nn_size
ocean_areawithcav  => mesh%ocean_areawithcav
coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)                  => mesh%coord_nod2D(:,:)
geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)              => mesh%geo_coord_nod2D(:,:)
elem2D_nodes(1:3, 1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%elem2D_nodes(:,:)
edges(1:2,1:myDim_edge2D+eDim_edge2D)                      => mesh%edges(:,:)
edge_tri(1:2,1:myDim_edge2D+eDim_edge2D)                   => mesh%edge_tri(:,:)
elem_edges(1:3,1:myDim_elem2D)                             => mesh%elem_edges(:,:)
elem_area(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)         => mesh%elem_area(:)
edge_dxdy(1:2,1:myDim_edge2D+eDim_edge2D)                  => mesh%edge_dxdy(:,:)
edge_cross_dxdy(1:4,1:myDim_edge2D+eDim_edge2D)            => mesh%edge_cross_dxdy(:,:)
elem_cos(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)          => mesh%elem_cos(:)
metric_factor(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%metric_factor(:)
elem_neighbors(1:3,1:myDim_elem2D)                         => mesh%elem_neighbors(:,:)
nod_in_elem2D      => mesh%nod_in_elem2D   ! (maxval(rmax),myDim_nod2D+eDim_nod2D)
x_corners          => mesh%x_corners   ! (myDim_nod2D, maxval(rmax))
y_corners          => mesh%y_corners   ! (myDim_nod2D, maxval(rmax))
nod_in_elem2D_num(1:myDim_nod2D+eDim_nod2D)                => mesh%nod_in_elem2D_num(:)
depth(1:myDim_nod2D+eDim_nod2D)                            => mesh%depth(:)
gradient_vec(1:6,1:myDim_elem2D)                           => mesh%gradient_vec(:,:)
gradient_sca(1:6,1:myDim_elem2D)                           => mesh%gradient_sca(:,:)
bc_index_nod2D(1:myDim_nod2D+eDim_nod2D)                   => mesh%bc_index_nod2D(:)
zbar(1:mesh%nl)                                            => mesh%zbar(:)
Z(1:mesh%nl-1)                                             => mesh%Z(:)
elem_depth         => mesh%elem_depth      ! never used, not even allocated
nlevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%nlevels(:)
nlevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%nlevels_nod2D(:)
nlevels_nod2D_min(1:myDim_nod2D+eDim_nod2D)                => mesh%nlevels_nod2D_min(:)
area(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)                   => mesh%area(:,:)
areasvol(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%areasvol(:,:)
area_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%area_inv(:,:)
areasvol_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)           => mesh%areasvol_inv(:,:)
mesh_resolution(1:myDim_nod2d+eDim_nod2D)                  => mesh%mesh_resolution(:)
ssh_stiff                                                  => mesh%ssh_stiff
lump2d_north(1:myDim_nod2d)                                => mesh%lump2d_north(:)
lump2d_south(1:myDim_nod2d)                                => mesh%lump2d_south(:)
cavity_flag_n(1:myDim_nod2D+eDim_nod2D)                    => mesh%cavity_flag_n(:)
cavity_flag_e(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%cavity_flag_e(:)
!!$cavity_lev_nod2D(1:myDim_nod2D+eDim_nod2D)                 => mesh%cavity_lev_nod2D
!!$cavity_lev_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%cavity_lev_elem2D
cavity_depth(1:myDim_nod2D+eDim_nod2D)                     => mesh%cavity_depth(:)
ulevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%ulevels(:)
ulevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%ulevels_nod2D(:)
ulevels_nod2D_max(1:myDim_nod2D+eDim_nod2D)                => mesh%ulevels_nod2D_max(:)
nn_num(1:myDim_nod2D)                                      => mesh%nn_num(:)
nn_pos(1:mesh%nn_size, 1:myDim_nod2D)                      => mesh%nn_pos(:,:)
hnode(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)               => mesh%hnode(:,:)
hnode_new(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)           => mesh%hnode_new(:,:)
zbar_3d_n(1:mesh%nl, 1:myDim_nod2D+eDim_nod2D)             => mesh%zbar_3d_n(:,:)
Z_3d_n(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)              => mesh%Z_3d_n(:,:)
helem(1:mesh%nl-1, 1:myDim_elem2D)                         => mesh%helem(:,:)
bottom_elem_thickness(1:myDim_elem2D)                      => mesh%bottom_elem_thickness(:)
bottom_node_thickness(1:myDim_nod2D+eDim_nod2D)            => mesh%bottom_node_thickness(:)
dhe(1:myDim_elem2D)                                        => mesh%dhe(:)
hbar(1:myDim_nod2D+eDim_nod2D)                             => mesh%hbar(:)
hbar_old(1:myDim_nod2D+eDim_nod2D)                         => mesh%hbar_old(:)
!zbar_n(1:mesh%nl)                                          => mesh%zbar_n
!Z_n(1:mesh%nl-1)                                           => mesh%Z_n
zbar_n_bot(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_bot(:)
zbar_e_bot(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_bot(:)
zbar_n_srf(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_srf(:)
zbar_e_srf(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_srf(:)
# 470 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

    l_init_zero=.true.
    if (present(o_init_zero)) then
       l_init_zero=o_init_zero
    end if
    if (l_init_zero) then
!$OMP PARALLEL DO
       do n=1, myDim_nod2D
          flux(:, n)=0.0_WP
       end do
!$OMP END PARALLEL DO
    end if

    ! --------------------------------------------------------------------------
    ! Vertical advection
    ! --------------------------------------------------------------------------
    ! A piecewise parabolic scheme for uniformly-spaced layers.
    ! See Colella and Woodward, JCP, 1984, 174-201. It can be coded so as to to take
    ! non-uniformity into account, but this is more cumbersome. This is the version for AB
    ! time stepping
    ! --------------------------------------------------------------------------
!   overshoot_counter=0
!   counter          =0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tvert, tv, aL, aR, aj, x, dzjm1, dzj, dzjp1, dzjp2, deltaj, deltajp1, n, nz, nzmax, nzmin)
!$OMP DO
    do n=1, myDim_nod2D
        !_______________________________________________________________________
        !Interpolate to zbar...depth levels --> all quantities (tracer ...) are
        ! calculated on mid depth levels
        ! nzmax ... number of depth levels at node n
        nzmax=nlevels_nod2D(n)
        nzmin=ulevels_nod2D(n)

        ! tracer at surface level
        tv(nzmin)=ttf(nzmin,n)
        ! tracer at surface+1 level
!       tv(2)=-ttf(1,n)*min(sign(1.0, W(2,n)), 0._WP)+ttf(2,n)*max(sign(1.0, W(2,n)), 0._WP)
!       tv(3)=-ttf(2,n)*min(sign(1.0, W(3,n)), 0._WP)+ttf(3,n)*max(sign(1.0, W(3,n)), 0._WP)
        tv(nzmin+1)=0.5*(ttf(nzmin,  n)+ttf(nzmin+1,n))
        ! tacer at bottom-1 level
        tv(nzmax-1)=-ttf(nzmax-2,n)*min(sign(1.0_wp, W(nzmax-1,n)), 0._WP)+ttf(nzmax-1,n)*max(sign(1.0_wp, W(nzmax-1,n)), 0._WP)
!       tv(nzmax-1)=0.5_WP*(ttf(nzmax-2,n)+ttf(nzmax-1,n))
        ! tracer at bottom level
        tv(nzmax)=ttf(nzmax-1,n)

        !_______________________________________________________________________
        ! calc tracer for surface+2 until depth-2 layer
        ! see Colella and Woodward, JCP, 1984, 174-201 --> equation (1.9)
        ! loop over layers (segments)
        !!PS do nz=3, nzmax-3
        do nz=nzmin+1, nzmax-3
            !___________________________________________________________________
            ! for uniform spaced vertical grids --> piecewise parabolic method (ppm)
            ! equation (1.9)
            ! tv(nz)=(7.0_WP*(ttf(nz-1,n)+ttf(nz,n))-(ttf(nz-2,n)+ttf(nz+1,n)))/12.0_WP

            !___________________________________________________________________
            ! for non-uniformity spaced vertical grids --> piecewise parabolic
            ! method (ppm) see see Colella and Woodward, JCP, 1984, 174-201
            ! --> full equation (1.6), (1.7) and (1.8)
            dzjm1    = hnode_new(nz-1,n)
            dzj      = hnode_new(nz  ,n)
            dzjp1    = hnode_new(nz+1,n)
            dzjp2    = hnode_new(nz+2,n)
            ! Be carefull here vertical operation have to be done on NEW vertical mesh !!!

            !___________________________________________________________________
            ! equation (1.7)
            ! --> Here deltaj is the average slope in the jth zone of the parabola
            !     with zone averages a_(j-1) and a_j, a_(j+1)
            ! --> a_j^n
            deltaj   = dzj/(dzjm1+dzj+dzjp1)* &
                      ( &
                       (2._WP*dzjm1+dzj    )/(dzjp1+dzj)*(ttf(nz+1,n)-ttf(nz  ,n)) +  &
                       (dzj    +2._WP*dzjp1)/(dzjm1+dzj)*(ttf(nz  ,n)-ttf(nz-1,n)) &
                      )
            ! --> a_(j+1)^n
            deltajp1 = dzjp1/(dzj+dzjp1+dzjp2)* &
                      ( &
                       (2._WP*dzj+dzjp1  )/(dzjp2+dzjp1)*(ttf(nz+2,n)-ttf(nz+1,n)) +  &
                       (dzjp1+2._WP*dzjp2)/(dzj  +dzjp1)*(ttf(nz+1,n)-ttf(nz  ,n)) &
                      )
            !___________________________________________________________________
            ! condition (1.8)
            ! --> This modification leads to a somewhat steeper representation of
            !     discontinuities in the solution. It also guarantees that a_(j+0.5)
            !     lies in the range of values defined by a_j; and a_(j+1);
            if ( (ttf(nz+1,n)-ttf(nz  ,n))*(ttf(nz  ,n)-ttf(nz-1,n)) > 0._WP ) then
                deltaj = min(  abs(deltaj), &
                             2._WP*abs(ttf(nz+1,n)-ttf(nz  ,n)),&
                             2._WP*abs(ttf(nz  ,n)-ttf(nz-1,n)) &
                             )*sign(1.0_WP,deltaj)
            else
                deltaj = 0.0_WP
            endif
            if ( (ttf(nz+2,n)-ttf(nz+1,n))*(ttf(nz+1,n)-ttf(nz  ,n)) > 0._WP ) then
                deltajp1 = min(  abs(deltajp1),&
                               2._WP*abs(ttf(nz+2,n)-ttf(nz+1,n)),&
                               2._WP*abs(ttf(nz+1,n)-ttf(nz,n)) &
                               )*sign(1.0_WP,deltajp1)
            else
                deltajp1 = 0.0_WP
            endif
            !___________________________________________________________________
            ! equation (1.6)
            ! --> calcualte a_(j+0.5)
            ! nz+1 is the interface betweel layers (segments) nz and nz+1
            tv(nz+1)=    ttf(nz,n) &
                        + dzj/(dzj+dzjp1)*(ttf(nz+1,n)-ttf(nz,n)) &
                        + 1._WP/(dzjm1+dzj+dzjp1+dzjp2) * &
                        ( &
                            (2._WP*dzjp1*dzj)/(dzj+dzjp1)* &
                                ((dzjm1+dzj)/(2._WP*dzj+dzjp1) - (dzjp2+dzjp1)/(2._WP*dzjp1+dzj))*(ttf(nz+1,n)-ttf(nz,n)) &
                        - dzj*(dzjm1+dzj)/(2._WP*dzj+dzjp1)*deltajp1 &
                        + dzjp1*(dzjp1+dzjp2)/(dzj+2._WP*dzjp1)*deltaj &
                        )
                       !tv(nz+1)=max(min(ttf(nz, n), ttf(nz+1, n)), min(max(ttf(nz, n), ttf(nz+1, n)), tv(nz+1)))
        end do ! --> do nz=2,nzmax-3

        tvert(1:nzmax)=0._WP
        ! loop over layers (segments)
        do nz=nzmin, nzmax-1
            if ((W(nz,n)<=0._WP) .AND. (W(nz+1,n)>=0._WP)) CYCLE
            !counter=counter+1
            aL=tv(nz)
            aR=tv(nz+1)
            if ((aR-ttf(nz, n))*(ttf(nz, n)-aL)<=0._WP) then
                !   write(*,*) aL, ttf(nz, n), aR
                !   overshoot_counter=overshoot_counter+1
                aL =ttf(nz, n)
                aR =ttf(nz, n)
            end if
            if ((aR-aL)*(ttf(nz, n)-0.5_WP*(aL+aR))> (aR-aL)**2/6._WP) then
                aL =3._WP*ttf(nz, n)-2._WP*aR
            end if
            if ((aR-aL)*(ttf(nz, n)-0.5_WP*(aR+aL))<-(aR-aL)**2/6._WP) then
                aR =3._WP*ttf(nz, n)-2._WP*aL
            end if

            dzj   = hnode(nz,n)
            aj=6.0_WP*(ttf(nz, n)-0.5_WP*(aL+aR))

            if (W(nz,n)>0._WP) then
                x=min(W(nz,n)*dt/dzj, 1._WP)
                tvert(nz  )=(-aL-0.5_WP*x*(aR-aL+(1._WP-2._WP/3._WP*x)*aj))
                tvert(nz  )=tvert(nz) ! compute 2nd moment for DVD
                tvert(nz  )=tvert(nz)*area(nz,n)*W(nz,n)
            end if

            if (W(nz+1,n)<0._WP) then
                x=min(-W(nz+1,n)*dt/dzj, 1._WP)
                tvert(nz+1)=(-aR+0.5_WP*x*(aR-aL-(1._WP-2._WP/3._WP*x)*aj))
                tvert(nz+1)=tvert(nz+1) ! compute 2nd moment for DVD
                tvert(nz+1)=tvert(nz+1)*area(nz+1,n)*W(nz+1,n)
            end if
        end do

        !_______________________________________________________________________
        ! Surface flux
        tvert(nzmin)= -tv(nzmin)*W(nzmin,n)*area(nzmin,n)
        ! Zero bottom flux
        tvert(nzmax)=0.0_WP
        flux(nzmin:nzmax, n)=tvert(nzmin:nzmax)-flux(nzmin:nzmax, n)
    end do ! --> do n=1, myDim_nod2D
!       if (mype==0) write(*,*) 'PPM overshoot statistics:', real(overshoot_counter)/real(counter)
!$OMP END DO
!$OMP END PARALLEL
end subroutine adv_tra_vert_ppm
!
!
!===============================================================================
subroutine adv_tra_ver_cdiff(w, ttf, partit, mesh, flux, o_init_zero)
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto
    implicit none
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
    logical, optional                 :: o_init_zero
    logical                           :: l_init_zero
    integer                           :: n, nz, nzmax, nzmin
    real(kind=WP)                     :: tvert(mesh%nl), tv

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_def.h" 1

  integer,          pointer     :: MPI_COMM_FESOM ! FESOM communicator (for ocean only runs if often a copy of MPI_COMM_WORLD)
  type(com_struct), pointer     :: com_nod2D
  type(com_struct), pointer     :: com_elem2D
  type(com_struct), pointer     :: com_elem2D_full
  integer                       :: ub, lb ! to work with r(s)_mpitype_elem3D(nod3D)

  integer, dimension(:),     pointer :: s_mpitype_edge2D, r_mpitype_edge2D
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D, r_mpitype_elem2D
  integer, dimension(:),     pointer :: s_mpitype_elem2D_full_i, r_mpitype_elem2D_full_i
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D_full,   r_mpitype_elem2D_full
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D,        r_mpitype_elem3D
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D_full,   r_mpitype_elem3D_full

  integer, dimension(:),     pointer :: s_mpitype_nod2D, r_mpitype_nod2D
  integer, dimension(:),     pointer :: s_mpitype_nod2D_i, r_mpitype_nod2D_i
  integer, dimension(:,:,:), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

  integer, pointer   :: MPIERR
  integer, pointer   :: npes
  integer, pointer   :: mype
  integer, pointer   :: maxPEnum

  integer, dimension(:), pointer     :: part

  ! Mesh partition
  integer, pointer                :: myDim_nod2D, eDim_nod2D
  integer, dimension(:), pointer  :: myList_nod2D
  integer, pointer                :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, dimension(:), pointer  :: myList_elem2D
  integer, pointer                :: myDim_edge2D, eDim_edge2D
  integer, dimension(:), pointer  :: myList_edge2D

  integer, pointer                :: pe_status

  integer, dimension(:), pointer  ::  remPtr_nod2D(:),  remList_nod2D(:)
  integer, dimension(:), pointer  ::  remPtr_elem2D(:), remList_elem2D(:)

  logical, pointer                :: elem_full_flag
# 658 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_def.h" 1
integer         , pointer :: nod2D
integer         , pointer :: elem2D
integer         , pointer :: edge2D
integer         , pointer :: edge2D_in
real(kind=WP)   , pointer :: ocean_area
real(kind=WP)   , pointer :: ocean_areawithcav
integer         , pointer :: nl
integer         , pointer :: nn_size
real(kind=WP), dimension(:,:), pointer :: coord_nod2D, geo_coord_nod2D
integer, dimension(:,:)      , pointer :: elem2D_nodes
integer, dimension(:,:)      , pointer :: edges
integer, dimension(:,:)      , pointer :: edge_tri
integer, dimension(:,:)      , pointer :: elem_edges
real(kind=WP), dimension(:)  , pointer :: elem_area
real(kind=WP), dimension(:,:), pointer :: edge_dxdy, edge_cross_dxdy
real(kind=WP), dimension(:)  , pointer :: elem_cos, metric_factor
integer,       dimension(:,:), pointer :: elem_neighbors
integer,       dimension(:,:), pointer :: nod_in_elem2D
real(kind=WP), dimension(:,:), pointer :: x_corners, y_corners
integer,       dimension(:)  , pointer :: nod_in_elem2D_num
real(kind=WP), dimension(:)  , pointer :: depth
real(kind=WP), dimension(:,:), pointer :: gradient_vec
real(kind=WP), dimension(:,:), pointer :: gradient_sca
integer,       dimension(:)  , pointer :: bc_index_nod2D
real(kind=WP), dimension(:)  , pointer :: zbar, Z, elem_depth
integer,       dimension(:)  , pointer :: nlevels, nlevels_nod2D, nlevels_nod2D_min
real(kind=WP), dimension(:,:), pointer :: area, area_inv, areasvol, areasvol_inv
real(kind=WP), dimension(:)  , pointer :: mesh_resolution
real(kind=WP), dimension(:)  , pointer :: lump2d_north, lump2d_south
type(sparse_matrix)          , pointer :: ssh_stiff
integer,       dimension(:)  , pointer :: cavity_flag_n, cavity_flag_e
real(kind=WP), dimension(:)  , pointer :: cavity_depth
integer,       dimension(:)  , pointer :: ulevels, ulevels_nod2D, ulevels_nod2D_max
integer,       dimension(:)  , pointer :: nn_num
integer,       dimension(:,:), pointer :: nn_pos

real(kind=WP), dimension(:,:), pointer :: hnode
real(kind=WP), dimension(:,:), pointer :: hnode_new
real(kind=WP), dimension(:,:), pointer :: zbar_3d_n
real(kind=WP), dimension(:,:), pointer :: Z_3d_n
real(kind=WP), dimension(:,:), pointer :: helem
real(kind=WP), dimension(:)  , pointer :: bottom_elem_thickness
real(kind=WP), dimension(:)  , pointer :: bottom_node_thickness
real(kind=WP), dimension(:)  , pointer :: dhe
real(kind=WP), dimension(:)  , pointer :: hbar
real(kind=WP), dimension(:)  , pointer :: hbar_old
!real(kind=WP), dimension(:)  , pointer :: zbar_n
!real(kind=WP), dimension(:)  , pointer :: Z_n
real(kind=WP), dimension(:)  , pointer :: zbar_n_bot
real(kind=WP), dimension(:)  , pointer :: zbar_e_bot
real(kind=WP), dimension(:)  , pointer :: zbar_n_srf
real(kind=WP), dimension(:)  , pointer :: zbar_e_srf
# 659 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_part_ass.h" 1
MPI_COMM_FESOM  => partit%MPI_COMM_FESOM
com_nod2D       => partit%com_nod2D
com_elem2D      => partit%com_elem2D
com_elem2D_full => partit%com_elem2D_full
myDim_nod2D     => partit%myDim_nod2D
eDim_nod2D      => partit%eDim_nod2D
myDim_elem2D    => partit%myDim_elem2D
eDim_elem2D     => partit%eDim_elem2D
eXDim_elem2D    => partit%eXDim_elem2D
myDim_edge2D    => partit%myDim_edge2D
eDim_edge2D     => partit%eDim_edge2D
pe_status       => partit%pe_status
elem_full_flag  => partit%elem_full_flag
MPIERR          => partit%MPIERR
npes            => partit%npes
mype            => partit%mype
maxPEnum        => partit%maxPEnum
part            => partit%part

lb=lbound(partit%s_mpitype_elem3D, 2)
ub=ubound(partit%s_mpitype_elem3D, 2)

myList_nod2D (1:myDim_nod2D +eDim_nod2D)                => partit%myList_nod2D(:)
myList_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)  => partit%myList_elem2D(:)
myList_edge2D(1:myDim_edge2D+eDim_edge2D)               => partit%myList_edge2D(:)

if (allocated(partit%remPtr_nod2D)) then
   remPtr_nod2D  (1:npes)                => partit%remPtr_nod2D(:)
   remList_nod2D (1:remPtr_nod2D(npes))  => partit%remList_nod2D(:)
end if

if (allocated(partit%remPtr_elem2D)) then
remPtr_elem2D (1:npes)                => partit%remPtr_elem2D(:)
remList_elem2D(1:remPtr_elem2D(npes)) => partit%remList_elem2D(:)
end if

s_mpitype_elem2D(1:com_elem2D%sPEnum, 1:4) => partit%s_mpitype_elem2D(:,:)
r_mpitype_elem2D(1:com_elem2D%rPEnum, 1:4) => partit%r_mpitype_elem2D(:,:)

s_mpitype_elem2D_full_i(1:com_elem2D_full%sPEnum) => partit%s_mpitype_elem2D_full_i(:)
r_mpitype_elem2D_full_i(1:com_elem2D_full%rPEnum) => partit%r_mpitype_elem2D_full_i(:)

s_mpitype_elem2D_full(1:com_elem2D_full%sPEnum, 1:4) => partit%s_mpitype_elem2D_full(:,:)
r_mpitype_elem2D_full(1:com_elem2D_full%rPEnum, 1:4) => partit%r_mpitype_elem2D_full(:,:)

s_mpitype_elem3D(1:com_elem2D%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D(:,:,:)
r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D(:,:,:)

s_mpitype_elem3D_full(1:com_elem2D_full%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D_full(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4)           => partit%r_mpitype_elem3D(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

s_mpitype_nod2D(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D(:)
r_mpitype_nod2D(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D(:)

s_mpitype_nod2D_i(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D_i(:)
r_mpitype_nod2D_i(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D_i(:)

s_mpitype_nod3D(1:com_nod2D%sPEnum, lb:ub, 1:3) => partit%s_mpitype_nod3D(:,:,:)
r_mpitype_nod3D(1:com_nod2D%rPEnum, lb:ub, 1:3) => partit%r_mpitype_nod3D(:,:,:)
# 660 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

# 1 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/associate_mesh_ass.h" 1
nod2D              => mesh%nod2D
elem2D             => mesh%elem2D
edge2D             => mesh%edge2D
edge2D_in          => mesh%edge2D_in
ocean_area         => mesh%ocean_area
nl                 => mesh%nl
nn_size            => mesh%nn_size
ocean_areawithcav  => mesh%ocean_areawithcav
coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)                  => mesh%coord_nod2D(:,:)
geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)              => mesh%geo_coord_nod2D(:,:)
elem2D_nodes(1:3, 1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%elem2D_nodes(:,:)
edges(1:2,1:myDim_edge2D+eDim_edge2D)                      => mesh%edges(:,:)
edge_tri(1:2,1:myDim_edge2D+eDim_edge2D)                   => mesh%edge_tri(:,:)
elem_edges(1:3,1:myDim_elem2D)                             => mesh%elem_edges(:,:)
elem_area(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)         => mesh%elem_area(:)
edge_dxdy(1:2,1:myDim_edge2D+eDim_edge2D)                  => mesh%edge_dxdy(:,:)
edge_cross_dxdy(1:4,1:myDim_edge2D+eDim_edge2D)            => mesh%edge_cross_dxdy(:,:)
elem_cos(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)          => mesh%elem_cos(:)
metric_factor(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%metric_factor(:)
elem_neighbors(1:3,1:myDim_elem2D)                         => mesh%elem_neighbors(:,:)
nod_in_elem2D      => mesh%nod_in_elem2D   ! (maxval(rmax),myDim_nod2D+eDim_nod2D)
x_corners          => mesh%x_corners   ! (myDim_nod2D, maxval(rmax))
y_corners          => mesh%y_corners   ! (myDim_nod2D, maxval(rmax))
nod_in_elem2D_num(1:myDim_nod2D+eDim_nod2D)                => mesh%nod_in_elem2D_num(:)
depth(1:myDim_nod2D+eDim_nod2D)                            => mesh%depth(:)
gradient_vec(1:6,1:myDim_elem2D)                           => mesh%gradient_vec(:,:)
gradient_sca(1:6,1:myDim_elem2D)                           => mesh%gradient_sca(:,:)
bc_index_nod2D(1:myDim_nod2D+eDim_nod2D)                   => mesh%bc_index_nod2D(:)
zbar(1:mesh%nl)                                            => mesh%zbar(:)
Z(1:mesh%nl-1)                                             => mesh%Z(:)
elem_depth         => mesh%elem_depth      ! never used, not even allocated
nlevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%nlevels(:)
nlevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%nlevels_nod2D(:)
nlevels_nod2D_min(1:myDim_nod2D+eDim_nod2D)                => mesh%nlevels_nod2D_min(:)
area(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)                   => mesh%area(:,:)
areasvol(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%areasvol(:,:)
area_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%area_inv(:,:)
areasvol_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)           => mesh%areasvol_inv(:,:)
mesh_resolution(1:myDim_nod2d+eDim_nod2D)                  => mesh%mesh_resolution(:)
ssh_stiff                                                  => mesh%ssh_stiff
lump2d_north(1:myDim_nod2d)                                => mesh%lump2d_north(:)
lump2d_south(1:myDim_nod2d)                                => mesh%lump2d_south(:)
cavity_flag_n(1:myDim_nod2D+eDim_nod2D)                    => mesh%cavity_flag_n(:)
cavity_flag_e(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%cavity_flag_e(:)
!!$cavity_lev_nod2D(1:myDim_nod2D+eDim_nod2D)                 => mesh%cavity_lev_nod2D
!!$cavity_lev_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%cavity_lev_elem2D
cavity_depth(1:myDim_nod2D+eDim_nod2D)                     => mesh%cavity_depth(:)
ulevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%ulevels(:)
ulevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%ulevels_nod2D(:)
ulevels_nod2D_max(1:myDim_nod2D+eDim_nod2D)                => mesh%ulevels_nod2D_max(:)
nn_num(1:myDim_nod2D)                                      => mesh%nn_num(:)
nn_pos(1:mesh%nn_size, 1:myDim_nod2D)                      => mesh%nn_pos(:,:)
hnode(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)               => mesh%hnode(:,:)
hnode_new(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)           => mesh%hnode_new(:,:)
zbar_3d_n(1:mesh%nl, 1:myDim_nod2D+eDim_nod2D)             => mesh%zbar_3d_n(:,:)
Z_3d_n(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)              => mesh%Z_3d_n(:,:)
helem(1:mesh%nl-1, 1:myDim_elem2D)                         => mesh%helem(:,:)
bottom_elem_thickness(1:myDim_elem2D)                      => mesh%bottom_elem_thickness(:)
bottom_node_thickness(1:myDim_nod2D+eDim_nod2D)            => mesh%bottom_node_thickness(:)
dhe(1:myDim_elem2D)                                        => mesh%dhe(:)
hbar(1:myDim_nod2D+eDim_nod2D)                             => mesh%hbar(:)
hbar_old(1:myDim_nod2D+eDim_nod2D)                         => mesh%hbar_old(:)
!zbar_n(1:mesh%nl)                                          => mesh%zbar_n
!Z_n(1:mesh%nl-1)                                           => mesh%Z_n
zbar_n_bot(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_bot(:)
zbar_e_bot(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_bot(:)
zbar_n_srf(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_srf(:)
zbar_e_srf(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_srf(:)
# 661 "/work/ka1298/k202167/fesom_edsl/fesom2/dwarf/dwarf_tracer/src/oce_adv_tra_ver.F90" 2

    l_init_zero=.true.
    if (present(o_init_zero)) then
       l_init_zero=o_init_zero
    end if
    if (l_init_zero) then
!$OMP PARALLEL DO
       do n=1, myDim_nod2D
          flux(:, n)=0.0_WP
       end do
!$OMP END PARALLEL DO
    end if

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, nzmax, nzmin, tv, tvert)
!$OMP DO
    do n=1, myDim_nod2D
        !_______________________________________________________________________
        nzmax=nlevels_nod2D(n)-1
        nzmin=ulevels_nod2D(n)

        !_______________________________________________________________________
        ! Surface flux
        tvert(nzmin)= -W(nzmin,n)*ttf(nzmin,n)*area(nzmin,n)

        !_______________________________________________________________________
        ! Zero bottom flux
        tvert(nzmax+1)=0.0_WP

        !_______________________________________________________________________
        ! Other levels
        do nz=nzmin+1, nzmax
            tv=0.5_WP*(ttf(nz-1,n)+ttf(nz,n))
            tvert(nz)= -tv*W(nz,n)*area(nz,n)
        end do

        !_______________________________________________________________________
        flux(nzmin:nzmax, n)=tvert(nzmin:nzmax)-flux(nzmin:nzmax, n)
    end do ! --> do n=1, myDim_nod2D
!$OMP END DO
!$OMP END PARALLEL
end subroutine adv_tra_ver_cdiff
