#ifdef PUSHE_KERNEL
#include "adios_stub_macro.h"
#else
#include "adios_macro.h"
#endif

subroutine xgc_read
  use eq_module
  use sml_module, only : sml_2pi,sml_concentric,sml_read_eq,sml_bp_mult,sml_bt_mult,sml_mype, sml_comm, &
       sml_read_eq_m3dc1, sml_bt_sign, sml_bp_sign
#ifdef FUSION_IO
  use fusion_io
  use EZspline_obj
  use EZspline
  use itp_module
  use hdf5
#endif
  implicit none
  include 'mpif.h'
  real (kind=8):: deltheta,R0,B0,a,maxpsi,r,z,epsilon,a_min
  real (kind=8):: theta
  real (kind=8), external :: datanh
  integer :: i,j,ij,cn
#ifdef FUSION_IO
  integer(HID_T) :: file_id, dset_id      ! File identifier for m3d-c1 H5 file
  real (kind=8) :: data_out(2)
  integer(HSIZE_T), dimension(2) :: data_dims
#endif
  integer :: isrc, imag, ivecp, eq_g_idum, m3dc1_bp_sign, m3dc1_bt_sign
  integer :: ipsi_lcfs, ipsi_axis, ir_axis, iz_axis, cs, idum
  real (kind=8) :: rr,zz,psi,dr,dz,dpsi, period, rdum
  real (kind=8) :: bfield(3), vecpot(3), xx(3)
  real (kind=8) :: eq_axis_psi, eq_g_rleft, eq_g_rdim, eq_g_zdim, eq_g_zmid, eq_g_rcentr
  character (len=80) :: eq_g_header
  real (kind=8), allocatable :: psi_mid(:), I_mid(:)
  integer :: end_flag, ierr, h5err
  type(EZspline1_r8) :: spl_I
  !
  real (kind=8), external :: interpol_spl_1d


  if(sml_read_eq) then
    ! read equlibrium
    if(sml_mype==0) then
      open(9,file=eq_filename, status='old')
      read(9,300) eq_header
      read(9,200) eq_mr, eq_mz, eq_mpsi
    endif

    call mpi_bcast(eq_mr,   1, MPI_INTEGER, 0, sml_comm, ierr)
    call mpi_bcast(eq_mz,   1, MPI_INTEGER, 0, sml_comm, ierr)
    call mpi_bcast(eq_mpsi, 1, MPI_INTEGER, 0, sml_comm, ierr)

    allocate(eq_psi_grid(eq_mpsi), eq_rgrid(eq_mr),eq_zgrid(eq_mz))
    allocate(eq_I(eq_mpsi), eq_psi_rz(eq_mr,eq_mz))

    if(sml_mype==0) then
      read(9,100) eq_min_r, eq_max_r, eq_min_z, eq_max_z
      !print *,'min_r, max_r, min_z, max_z', eq_min_r,eq_max_r, eq_min_z, eq_max_z
      read(9,100) eq_axis_r, eq_axis_z,eq_axis_b  ! format changed 2006/02/24 - axis_z added
      !print *,'axis_r, axis_b', eq_axis_r, eq_axis_b
      read(9,100) eq_x_psi, eq_x_r, eq_x_z
      !print *, 'x_psi, x_r, x_z', eq_x_psi, eq_x_r, eq_x_z
      read(9,100) (eq_psi_grid(i), i=1, eq_mpsi)
      !print *, 'psi_grid', eq_psi_grid
      !do i=1, eq_mpsi
      !  write(112,*) i, eq_psi_grid(i)
      !enddo
      !close(112)

      read(9,100) (eq_I(i), i=1,eq_mpsi)
      !read(9,100) (eq_rgrid(i), i=1, eq_mr)
      !read(9,100) (eq_zgrid(i), i=1, eq_mz)
      read(9,100) ((eq_psi_rz(i,j) , i=1, eq_mr),j=1, eq_mz)
      read(9,200) end_flag
      print *, 'eq end_flag', end_flag
      if(end_flag/=-1) then
        !additional eq. input
        print *, 'wrong file format'
        stop
      endif
      close(9)
    endif

    ! Broadcast
    call mpi_bcast(eq_min_r, 1, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_max_r, 1, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_min_z, 1, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_max_z, 1, MPI_REAL8, 0, sml_comm, ierr)
    !
    call mpi_bcast(eq_axis_r, 1, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_axis_z, 1, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_axis_b, 1, MPI_REAL8, 0, sml_comm, ierr)
    !
    call mpi_bcast(eq_x_psi, 1, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_x_r, 1, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_x_z, 1, MPI_REAL8, 0, sml_comm, ierr)
    !
    call mpi_bcast(eq_psi_grid, eq_mpsi, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_I, eq_mpsi, MPI_REAL8, 0, sml_comm, ierr)
    call mpi_bcast(eq_psi_rz, eq_mr*eq_mz, MPI_REAL8, 0, sml_comm, ierr)

     do i=1, eq_mr
      eq_rgrid(i)=eq_min_r+(eq_max_r-eq_min_r)/real(eq_mr-1) * real(i-1)
     enddo
     do i=1, eq_mz
        eq_zgrid(i)=eq_min_z+(eq_max_z-eq_min_z)/real(eq_mz-1) * real(i-1)
     enddo

#ifdef FUSION_IO
  ! Read equilibrium from m3dc1 output:
  !
  ! The input should contain
  !     sml_read_eq_m3dc1 = .true.
  !     eq_m3dc1_filename       = <M3D-C1 file>
  !     eq_filename = <name of eqd file based on the M3D-C1 input written by this routine>
  !     eq_g_filename     = <g-file generated from M3D-C1 file>
  !
  ! The work flow is
  !     1) experiment EFIT g-file is used as input to M3D-C1 to generate equilibrium
  !     2) the original EFIT g-file and the M3D-C1 output files are used to generate new g-file with a tool such as m3dc12g
  !     3) the new M3D-C1 g-file is used to generate XGC mesh
  !     4) the new XGC mesh along with the new M3D-C1 g-file and the M3D-C1 output files are specified in the XGC input file
  !
  ! The reasons to specify the newly generated M3D-C1 g-file as well as the M3D-C1 output files are
  !     1) to read the limiter data, which is copied from the original experimental EFIT file
  !     2) so that the resolution which was used to generate XGC mesh is the same as the one used to generate XGC fields
  !
  else if (sml_read_eq_m3dc1) then

     if (sml_mype==0) then
        open(9,file=eq_g_filename, status='old')
        read(9,301) eq_g_header, eq_g_idum, eq_mr, eq_mz
        read(9,202) eq_g_rdim, eq_g_zdim, eq_g_rcentr, eq_g_rleft, eq_g_zmid
        close(9)

        ! resolution setting
        eq_mpsi = eq_mr

        ! domain limits
        eq_min_r = eq_g_rleft
        eq_max_r = eq_g_rleft + eq_g_rdim
        eq_min_z = eq_g_zmid  - eq_g_zdim/2d0
        eq_max_z = eq_min_z   + eq_g_zdim
     endif

     call mpi_bcast(eq_mr,   1, MPI_INTEGER, 0, sml_comm, ierr)
     call mpi_bcast(eq_mz,   1, MPI_INTEGER, 0, sml_comm, ierr)
     call mpi_bcast(eq_mpsi, 1, MPI_INTEGER, 0, sml_comm, ierr)

     allocate(eq_psi_grid(eq_mpsi), eq_rgrid(eq_mr),eq_zgrid(eq_mz))
     allocate(eq_I(eq_mpsi), eq_psi_rz(eq_mr,eq_mz))

     call mpi_bcast(eq_min_r, 1, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_max_r, 1, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_min_z, 1, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_max_z, 1, MPI_REAL8, 0, sml_comm, ierr)

     ! UNIFORM GRID
     do i=1, eq_mr
        eq_rgrid(i)=eq_min_r+(eq_max_r-eq_min_r)/real(eq_mr-1) * real(i-1)
     enddo
     do i=1, eq_mz
        eq_zgrid(i)=eq_min_z+(eq_max_z-eq_min_z)/real(eq_mz-1) * real(i-1)
     enddo

     if (sml_mype==0) then
        ! location of x-point, must use hdf5 for this as fusion-io doesn't have it yet
        call h5open_f(h5err)
        call h5fopen_f(trim(eq_m3dc1_filename), H5F_ACC_RDWR_F, file_id, h5err)
        call h5dopen_f(file_id,"/scalars/xnull",dset_id, h5err)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, h5err)
        eq_x_r = data_out(1)
        call h5dopen_f(file_id,"/scalars/znull",dset_id, h5err)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, h5err)
        eq_x_z = data_out(1)
        call h5dclose_f(dset_id, h5err)
        call h5fclose_f(file_id, h5err)
        call h5close_f(h5err)


        ! OPEN FILE
        call fio_open_source_f(FIO_M3DC1_SOURCE, trim(eq_m3dc1_filename), isrc, ierr)
        write(*,*) '-- M3D-C1File ', trim(eq_m3dc1_filename), ' opened'

        ! Print information about coordinate system
        call fio_get_int_parameter_f(isrc, FIO_GEOMETRY, cs, ierr)
        call fio_get_real_parameter_f(isrc, FIO_PERIOD, period, ierr)
        if(cs.eq.FIO_CARTESIAN) then
           print *, 'xgc_read: Using CARTESIAN coordinate system'
           print *, 'xgc_read: Toroidal period (in m) = ', period
           stop 'xgc_read: Cartesian coordinates are not supported for M3D-C1 input!'
        else if(cs.eq.FIO_CYLINDRICAL) then
           print *, 'xgc_read: Using CYLINDRICAL coordinate system'
           print *, 'xgc_read: Toroidal period (in rad) = ', period
        else
           print *, 'xgc_read: ERROR: Unrecognized coordinate system'
           stop
        end if

        ! GET VALUES ON MAGNETIC AXIS AND LAST FLUX SURFACE
        call fio_get_series_f(isrc,FIO_LCFS_PSI, ipsi_lcfs, ierr)
        psi=0D0
        call fio_eval_series_f(ipsi_lcfs, psi, eq_x_psi, ierr)
        write(*,*) '-- Value of last close flux surface read, eq_x_psi = ', eq_x_psi

        call fio_get_series_f(isrc,FIO_MAGAXIS_PSI, ipsi_axis, ierr)
        call fio_eval_series_f(ipsi_axis, psi, eq_axis_psi, ierr)
        write(*,*) '-- Value of psi on magnetic axis read, eq_axis_psi = ', eq_axis_psi

        call fio_get_series_f(isrc,FIO_MAGAXIS_R, ir_axis, ierr)
        call fio_eval_series_f(ir_axis, psi, eq_axis_r, ierr)
        write(*,*) '-- Value of R on magnetic axis read, eq_axis_r = ',eq_axis_r

        call fio_get_series_f(isrc,FIO_MAGAXIS_Z, iz_axis, ierr)
        call fio_eval_series_f(iz_axis, psi, eq_axis_z, ierr)
        write(*,*) '-- Value of Z on magnetic axis read, eq_axis_z = ',eq_axis_z

        ! SETUP FIELDS
        ! set option to read equilibrium fields only
        call fio_get_options_f(isrc, ierr)
        idum=0
        call fio_set_int_option_f(FIO_TIMESLICE, idum, ierr)
        rdum=1D0
        call fio_set_real_option_f(FIO_LINEAR_SCALE, rdum, ierr)
        call fio_set_int_option_f(FIO_PART, FIO_EQUILIBRIUM_ONLY, ierr)
        write(*,*) '-- Set up to read equilibrium quantities'

        ! magnetic field and total pressure are species-independent
        call fio_get_field_f(isrc, FIO_MAGNETIC_FIELD, imag, ierr)
        call fio_get_field_f(isrc, FIO_VECTOR_POTENTIAL, ivecp, ierr)
        write(*,*) '-- Read magnetic field and vector potential'

        ! GET FIELD QUANTITIES
        xx=(/eq_axis_r, 0d0, eq_axis_z/)
        call fio_eval_field_f(imag, xx, bfield, ierr)
        eq_axis_b = sqrt(sum(bfield**2))

        do j = 1,eq_mz
           do i = 1,eq_mr
              rr = eq_rgrid(i)
              zz = eq_zgrid(j)

              ! poloidal flux found from vector potential
              xx=(/rr, 0d0, zz/)
              call fio_eval_field_f(ivecp,xx,vecpot,ierr)
              eq_psi_rz(i,j) = vecpot(2) * rr
           enddo
        enddo

        ! CALCULATION OF eq_I
        ! - calculate I and psi at the outboard midplane, and then interpolate I to eq_psi_grid
        allocate(I_mid(eq_mr), psi_mid(eq_mr))
        dr = (eq_max_r - eq_axis_r)/real(eq_mr-1.0)
        do i = 1,eq_mr
           rr = eq_axis_r + dr*real(i-1)
           xx=(/rr, 0d0, eq_axis_z/)
           call fio_eval_field_f(imag, xx, bfield, ierr)
           I_mid(i) = bfield(2) * rr

           call fio_eval_field_f(ivecp,xx ,vecpot ,ierr)
           psi_mid(i) = vecpot(2) * rr
        enddo

        if (eq_x_psi<eq_axis_psi) then
           write(*,*) 'xgc_read: CAUTION: psi_sep < psi_axis ---> Inverting psi'
           m3dc1_bp_sign=-1
           eq_psi_rz = eq_axis_psi - eq_psi_rz
           psi_mid = eq_axis_psi - psi_mid
           eq_x_psi = eq_axis_psi - eq_x_psi
           eq_axis_psi = 0d0
        else
           m3dc1_bp_sign=1
           eq_psi_rz = eq_psi_rz - eq_axis_psi
           psi_mid = psi_mid -eq_axis_psi
           eq_x_psi = eq_x_psi - eq_axis_psi
           eq_axis_psi = 0D0
        endif
        do i=1, eq_mpsi
           eq_psi_grid(i)=eq_axis_psi+(eq_x_psi-eq_axis_psi)/real(eq_mpsi-1) * real(i-1)
        enddo

        if (I_mid(1) .lt. 0D0) then
           write(*,*) 'xgc_read: CAUTION: I(psi) < 0 ---> Inverting I(psi)'
           m3dc1_bt_sign=-1
           I_mid = -I_mid
        else
           m3dc1_bt_sign=1
        endif

        if (m3dc1_bp_sign .ne. sml_bp_sign .or.  &
            m3dc1_bt_sign .ne. sml_bt_sign) then
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*) 'xgc_read: WARNING: Magnetic field direction not consistent with M3D-C1 data!'
           write(*,*) 'xgc_read: WARNING: Make sure that this is really what you want!'
           write(*,'(a,i2,a,i2,a)') 'xgc_read: WARNING: (m3dc1_bp_sign,sml_bp_sign)=(',m3dc1_bp_sign,',',floor(sml_bp_sign),')'
           write(*,'(a,i2,a,i2,a)') 'xgc_read: WARNING: (m3dc1_bt_sign,sml_bt_sign)=(',m3dc1_bt_sign,',',floor(sml_bt_sign),')'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        endif

        call init_1d_interpolation(spl_I,psi_mid,I_mid,eq_mr,.false.)
        do i = 1,eq_mpsi
           psi = eq_psi_grid(i)
           eq_I(i)=interpol_spl_1d(psi,0,spl_I)
        enddo
        call finalize_1d_interpolation(spl_I)

        deallocate(I_mid, psi_mid)

        ! CLEAN UP
        call fio_close_field_f(imag, ierr)
        call fio_close_field_f(ivecp, ierr)
        call fio_close_source_f(isrc, ierr)

        ! Write eqd output file for later post-processing
        open(unit=4211, file=eq_filename,action='write',status='replace')
        write(4211,300) 'EQD file based on M3D-C1 output'
        write(4211,200) eq_mr, eq_mz, eq_mpsi
        write(4211,100) eq_min_r, eq_max_r, eq_min_z, eq_max_z
        write(4211,100) eq_axis_r, eq_axis_z,eq_axis_b
        write(4211,100) eq_x_psi, eq_x_r, eq_x_z
        write(4211,100) (eq_psi_grid(i), i=1, eq_mpsi)
        write(4211,100) (eq_I(i), i=1,eq_mpsi)
        write(4211,100) ((eq_psi_rz(i,j) , i=1, eq_mr),j=1, eq_mz)
        write(4211,200) -1
        write(4211,1101) 1
        write(4211,1102) 0D0, 0D0
        write(4211,200) -1
        write(4211,1101) 1
        write(4211,1102) 0D0, 0D0
        write(4211,200) -1
        close(4211)
     endif

     call mpi_bcast(eq_axis_r, 1, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_axis_z, 1, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_axis_b, 1, MPI_REAL8, 0, sml_comm, ierr)
     !
     call mpi_bcast(eq_x_psi, 1, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_x_r, 1, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_x_z, 1, MPI_REAL8, 0, sml_comm, ierr)
     !
     call mpi_bcast(eq_psi_grid, eq_mpsi, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_I, eq_mpsi, MPI_REAL8, 0, sml_comm, ierr)
     call mpi_bcast(eq_psi_rz, eq_mr*eq_mz, MPI_REAL8, 0, sml_comm, ierr)
     !
     if (sml_mype==0) then
        write(*,*) '--Finished reading in M3D-C1 file'
     endif
#endif

  ! concentric flux surface from analytic formula is removed.
  else
     print *, 'No equlibrium data'
     stop
  endif

  ! B field multiplication
  eq_x_psi=sml_bp_mult*eq_x_psi
  eq_psi_grid=sml_bp_mult*eq_psi_grid
  eq_psi_rz=sml_bp_mult*eq_psi_rz
  eq_I=sml_bt_mult*eq_I
  eq_axis_b=sml_bt_mult*eq_axis_b

  if(sml_mype==0) then
     call write_eq_used
  endif

100  format(4(e19.13,1x))
200  format(8I8)
201  format(3I8)
202  format(5(e16.9))
300  format(a80)
301  format(a48,3i4)
1100 format (a2)
1101 format (i8)
1102 format (e19.13,1x,e19.13)

end subroutine xgc_read

! write out the equilibrim file (which was read or generated)
! read-in eq file and write-out eq file should be identical
subroutine write_eq_used
  use eq_module
  use sml_module
  implicit none
  integer :: i,j
  integer*8 :: buf_id, buf_size, total_size
  integer :: err

  buf_size = 3 * 4 + 10 * 8 + 2 * eq_mpsi * 8 + eq_mr * eq_mz * 8
  ADIOS_OPEN(buf_id,'eq_used','xgc.equil.bp','w',sml_comm_null,err)
  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
  ADIOS_WRITE(buf_id,eq_mr,err)
  ADIOS_WRITE(buf_id,eq_mz,err)
  ADIOS_WRITE(buf_id,eq_mpsi,err)
  ADIOS_WRITE(buf_id,eq_min_r,err)
  ADIOS_WRITE(buf_id,eq_max_r,err)
  ADIOS_WRITE(buf_id,eq_min_z,err)
  ADIOS_WRITE(buf_id,eq_max_z,err)
  ADIOS_WRITE(buf_id,eq_axis_r,err)
  ADIOS_WRITE(buf_id,eq_axis_z,err)
  ADIOS_WRITE(buf_id,eq_axis_b,err)
  ADIOS_WRITE(buf_id,eq_x_psi,err)
  ADIOS_WRITE(buf_id,eq_x_r,err)
  ADIOS_WRITE(buf_id,eq_x_z,err)
  ADIOS_WRITE(buf_id,eq_psi_grid,err)
  ADIOS_WRITE(buf_id,eq_I,err)
  ADIOS_WRITE(buf_id,eq_psi_rz,err)
  ADIOS_CLOSE(buf_id,err)

  !if (sml_mype==0) then
  !  open(11, file='fort.used.eq', status='replace')
  !  write(11,300) eq_header
  !  write(11,200) eq_mr, eq_mz, eq_mpsi
  !  write(11,100) eq_min_r, eq_max_r, eq_min_z, eq_max_z
  !  write(11,100) eq_axis_r, eq_axis_z, eq_axis_b
  !  write(11,100) eq_x_psi, eq_x_r, eq_x_z
  !  write(11,100) ( eq_psi_grid(i), i=1, eq_mpsi )
  !  write(11,100) ( eq_I(i), i=1, eq_mpsi)
  !  !write(9,100) (eq_rgrid(i)*sml_norm_r, i=1, eq_mr)
  !  !write(9,100) (eq_zgrid(i)*sml_norm_r, i=1, eq_mz
  !  write(11,100) ((eq_psi_rz(i,j), i=1, eq_mr),j=1, eq_mz)
  !  write(11,200) -1

  !  !!!!!!!!! Limiter and separatrix information should be added

  !  close(11)
  !endif

100 format(4(e19.13,1x))
200 format(8I8)
300 format(a80)

end subroutine



