subroutine xgc_read
  use eq_module
  use sml_module, only : sml_2pi,sml_concentric,sml_read_eq,sml_bp_mult,sml_bt_mult,sml_mype
  implicit none
  real (kind=8):: deltheta,R0,B0,a,maxpsi,r,z,epsilon,a_min
  real (kind=8):: theta
  real (kind=8), external :: datanh
  integer :: i,j,ij,cn
  integer :: end_flag


  if(sml_read_eq) then
     ! read equlibrium
     open(9,file=eq_filename, status='old')
     read(9,300) eq_header
     read(9,200) eq_mr, eq_mz, eq_mpsi

     allocate(eq_psi_grid(eq_mpsi), eq_rgrid(eq_mr),eq_zgrid(eq_mz))
     allocate(eq_I(eq_mpsi), eq_psi_rz(eq_mr,eq_mz))

     read(9,100) eq_min_r, eq_max_r, eq_min_z, eq_max_z
!     print *,'min_r, max_r, min_z, max_z', eq_min_r,eq_max_r, eq_min_z, eq_max_z
     read(9,100) eq_axis_r, eq_axis_z,eq_axis_b  ! format changed 2006/02/24 - axis_z added
!     print *,'axis_r, axis_b', eq_axis_r, eq_axis_b
     read(9,100) eq_x_psi, eq_x_r, eq_x_z
!     print *, 'x_psi, x_r, x_z', eq_x_psi, eq_x_r, eq_x_z
     read(9,100) (eq_psi_grid(i), i=1, eq_mpsi)
!     print *, 'psi_grid', eq_psi_grid
!     do i=1, eq_mpsi
!       write(112,*) i, eq_psi_grid(i)
!    enddo
!    close(112)

     read(9,100) (eq_I(i), i=1,eq_mpsi)
     !     read(9,100) (eq_rgrid(i), i=1, eq_mr)
     !     read(9,100) (eq_zgrid(i), i=1, eq_mz)
     read(9,100) ((eq_psi_rz(i,j) , i=1, eq_mr),j=1, eq_mz)
     read(9,200) end_flag
     if(sml_mype==0) print *, 'eq end_flag', end_flag
     if(end_flag/=-1) then
        !additional eq. input
        print *, 'wrong file format'
        stop
     endif

     ! read limiter -----------------------------------------------------------------------
!!$     read(9,1101) lim_mdata
!!$     allocate(lim_org_r(lim_mdata),lim_org_z(lim_mdata))
!!$
!!$     do i=1, lim_mdata
!!$        read(9,1102) lim_org_r(i),lim_org_z(i)
!!$        !   print *, lim_org_r(i), lim_org_z(i)
!!$     enddo
!!$
!!$     read(9,200) end_flag
!!$     if(sml_mype==0) print *, 'eq end_flag of lim_data', end_flag
!!$     if(end_flag/=-1) then
!!$        !additional eq. input
!!$        print *, 'wrong file format'
!!$        stop
!!$     endif
!!$
!!$     ! read separatrix
!!$     read(9,200) neu_sep_mtheta_file
!!$     allocate(neu_sep_r_file(neu_sep_mtheta_file+1),neu_sep_z_file(neu_sep_mtheta_file+1))
!!$     do i=1, neu_sep_mtheta_file
!!$        read(9,1102) neu_sep_r_file(i),neu_sep_z_file(i)
!!$     enddo
!!$     neu_sep_r_file(neu_sep_mtheta_file+1)=neu_sep_r_file(1)
!!$     neu_sep_z_file(neu_sep_mtheta_file+1)=neu_sep_z_file(1)
!!$
!!$
!!$     read(9,200) end_flag
!!$     if(sml_mype==0) print *, 'eq end_flag of sep_data', end_flag
!!$     if(end_flag/=-1) then
!!$        !additional eq. input
!!$        print *, 'wrong file format'
!!$        stop
!!$     endif

     close(9)

     do i=1, eq_mr
        eq_rgrid(i)=eq_min_r+(eq_max_r-eq_min_r)/real(eq_mr-1) * real(i-1)
     enddo
     do i=1, eq_mz
        eq_zgrid(i)=eq_min_z+(eq_max_z-eq_min_z)/real(eq_mz-1) * real(i-1)
     enddo

100  format(4(e19.13,1x))
200  format(8I8)
201  format(3I8)
300  format(a80)
1100 format (a2)
1101 format (i8)
1102 format (e19.13,1x,e19.13)
  ! concentric flux surface
  else if(sml_concentric) then
     eq_header="concentric example 1"
     R0=1.70D0   !meter
     B0=1.91D0/R0 !Tesla
     a=0.358D0
     !a_min=a*0.1D0
     !q0=0.854
     !q1=0
     !q2=2.184/a^2  !

     !q2=(3.-q0)/a^2 --ATC

     ! read equilibrium file
     ! allocate eq_module
     eq_mr=150
     eq_mz=150
     eq_mpsi=40
     eq_sep=0
!     eq_x_psi=1D40 x-psi defined later - last mesh flux
     eq_x_z=-1D40
     eq_min_r=R0-a*2.
     eq_min_z=-a*2.
     eq_max_r=R0+a*2.
     eq_max_z=a*2.


     eq_axis_r=R0
     eq_axis_z=0D0
     eq_axis_b=b0

     allocate(eq_psi_grid(eq_mpsi), eq_rgrid(eq_mr),eq_zgrid(eq_mz))
     allocate(eq_I(eq_mpsi), eq_psi_rz(eq_mr,eq_mz))

     do i=1, eq_mr
        eq_rgrid(i)=eq_min_r+(eq_max_r-eq_min_r)/real(eq_mr-1) * real(i-1)
     enddo
     do i=1, eq_mz
        eq_zgrid(i)=eq_min_z+(eq_max_z-eq_min_z)/real(eq_mz-1) * real(i-1)
     enddo
!     print *, eq_rgrid
!     print *, eq_zgrid
     do i=1, eq_mr
        do j=1, eq_mz
           r=eq_rgrid(i)
           z=eq_zgrid(j)
           epsilon=sqrt((r-r0)**2 + z**2)/r0
!           print *, epsilon
!           print *,'aa', 0.975846*dSqrt(1.-epsilon**2)
!           print *,'bb', datanh(0.975846*dSqrt(1.-epsilon**2)),datanh(0.88195d0)
           eq_psi_rz(i,j)=B0*R0**2 * ( &
!                0.0614001-0.0239453*datanh(0.988217*sqrt(1.-epsilon**2)) ) --ATC
                0.126108-0.0572657*datanh(0.975846*sqrt(1.-epsilon**2))  )

        enddo
     enddo
     epsilon=a/r0
     eq_x_psi=B0*R0**2 * ( &
          0.126108-0.0572657*datanh(0.975846*sqrt(1.-epsilon**2)) )

     maxpsi=eq_psi_rz(eq_mr,eq_mz)
!     print *, maxpsi
     do i=1, eq_mpsi
        eq_psi_grid(i)=maxpsi*real(i-1)/real(eq_mpsi-1)
     enddo

!     print *, eq_psi_grid(1:2)
     eq_I=R0*B0

     !for neu separatrix information
     epsilon=a/r0
!     do i=1, neu_sep_mtheta
!        theta=sml_2pi*0.01D0*real(i-1)
!        neu_sep_r(i)=1+ epsilon * cos(theta)
!        neu_sep_z(i)=sin(theta)
!     enddo

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

end subroutine xgc_read