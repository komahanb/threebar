! =============================================================================
!
!                            Main driver program
!
! =============================================================================


! Three Bar Truss Optimization
! Author : Komahan Boopathy, University of Dayton, OH 24/9/2013
! Thanks to IPOPT


program Threebar
!
      implicit none
!
!     include the Ipopt return codes
!
      include 'IpReturnCodes.inc'
!
!     Size of the problem (number of variables and equality constraints)
!
      integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
      parameter  (N = 3, M = 8, NELE_JAC = 24, NELE_HESS = 6)
      parameter  (IDX_STY = 1 )
!
!     Space for multipliers and constraints
!
      double precision LAM(M)
      double precision G(M)
!
!     Vector of variables
!
      double precision X(N)
!
!     Vector of lower and upper bounds
!
      double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
      double precision G_L(M), G_U(M)
!
!     Private data for evaluation routines
!     This could be used to pass double precision and integer arrays untouched
!     to the evaluation subroutines EVAL_*
!
      double precision DAT(20)
      integer IDAT(20)
!
!     Place for storing the Ipopt Problem Handle
!
      integer*8 IPROBLEM
      integer*8 IPCREATE
!
      integer IERR
      integer IPSOLVE, IPADDSTROPTION
      integer IPADDNUMOPTION, IPADDINTOPTION
      integer IPOPENOUTPUTFILE
!
      double precision F,pi
      integer i

      double precision  infbound
      parameter        (infbound = 1.d+20)

!
!     The following are the Fortran routines for computing the model
!     functions and their derivatives - their code can be found further
!     down in this file.
!
      external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS
!!
!!     The next is an optional callback method.  It is called once per
!!     iteration.
!!
      external ITER_CB
!
!     Set initial point and bounds:
!

      pi=4.0*atan(1.0) ! constant for later use (visible globally)
      
      !Problem data and other constants
      dat(1)=10.0 !length
      dat(2)=1.0e7 !E
      dat(3)=0.1 !gamma
      dat(4)=135.0*pi/180.0
      dat(5)=20000.0

!      I am thinking of changing "theta" between [0,180] degree to model uncertainty in operating environment (like a see-saw/ride)

      !!! Max constraint values
      !Tensile
      dat(6)=5000.0    ! psi
      dat(7)=20000.0   ! psi
      dat(8)=5000.0    ! psi  
      !Compressive
      dat(9)=5000.0    ! psi
      dat(10)=20000.0  ! psi
      dat(11)=5000.0   ! psi
      !Displacement
      dat(12)=0.005 ! in
      dat(13)=0.005 ! in

!!$      tensile_sigma1_max=dat(6)      
!!$      tensile_sigma2_max=dat(7)
!!$      tensile_sigma3_max=dat(8)
!!$
!!$      comp_sigma1_max=dat(9)      
!!$      comp_sigma2_max=dat(10)
!!$      comp_sigma3_max=dat(11)
!!$
!!$      max_u_disp=dat(12)
!!$      max_v_disp=dat(13)
!!$      

      do i=1,N
         X(i)   = 1.0     !500 mm
         X_L(i) = 0.010  !0 mm
         X_U(i) = infbound !max 1 metre =1000 mm
      end do
!
!     Set bounds for the constraints
!
      do i=1,M
         G_L(i)=-infbound
         G_U(i)=0.d0
      end do
!
!     First create a handle for the Ipopt problem (and read the options
!     file)
!

      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating an Ipopt Problem handle.'
         stop
      endif
!
!     Open an output file
!
      IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 50)
      if (IERR.ne.0 ) then
         write(*,*) 'Error opening the Ipopt output file.'
         goto 9000
      endif

!
!!
!!     Set a callback function to give you control once per iteration.
!!     You can use it if you want to generate some output, or to stop
!!     the optimization early.
!!
      call IPSETCALLBACK(IPROBLEM, ITER_CB)

!
!     As a simple example, we pass the constants in the constraints to
!     the EVAL_C routine via the "private" DAT array.
!

!      DAT(2) = 0.d0
!
!     Call optimization routine
!
      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)
!
!     Output:
!
      if( IERR.eq.IP_SOLVE_SUCCEEDED ) then
         write(*,*)
         write(*,*) 'The solution was found.'
         write(*,*)
         write(*,*) 'The final value of the objective function is ',F
         write(*,*)
         write(*,*) 'The optimal values of X are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'X  (',i,') = ',X(i)
         enddo
!!$         write(*,*)
!!$         write(*,*) 'The multipliers for the lower bounds are:'
!!$         write(*,*)
!!$         do i = 1, N
!!$            write(*,*) 'Z_L(',i,') = ',Z_L(i)
!!$         enddo
!!$         write(*,*)
!!$         write(*,*) 'The multipliers for the upper bounds are:'
!!$         write(*,*)
!!$         do i = 1, N
!!$            write(*,*) 'Z_U(',i,') = ',Z_U(i)
!!$         enddo
         write(*,*)
         write(*,*) 'The multipliers for the equality constraints are:'
         write(*,*)
         do i = 1, M
            write(*,*) 'LAM(',i,') = ',LAM(i)
         enddo

      else
         write(*,*)
         write(*,*) 'An error occoured.'
         write(*,*) 'The error code is ',IERR
         write(*,*)
      endif
!
9000  continue
!
!     Clean up
!
      call IPFREE(IPROBLEM)
      stop
!
9990  continue
      write(*,*) 'Error setting an option'
      goto 9000
    end program Threebar
!
! =============================================================================
!
!                    Computation of objective function
!
! =============================================================================
!
      subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X,I
      double precision F, X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      real *8::gamma, L
!!$
!!$      pi=4.0*atan(1.0) ! constant for later use (visible globally)
!!$      
!!$      dat(1)=10.0 !length
!!$      dat(2)=1.0e7 !E
!!$      dat(3)=0.1 !gamma
!!$      dat(4)=135*pi/180.0
!!$      dat(5)=20000.0

      gamma=dat(3)
      L=dat(1)
       
      F= x(1)*gamma*L*sqrt(2.0) + x(2)*gamma*L +  x(3)*gamma*L*sqrt(2.0)
      
      IERR = 0

      return
    end subroutine EV_F
!
! =============================================================================
!
!                Computation of gradient of objective function
!
! =============================================================================
!
      subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X,i
      double precision GRAD(N), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      real*8::L,gamma


      gamma=dat(3)
      L=dat(1)
      
      grad(1) = L*sqrt(2.0)*gamma
      grad(2) = L*gamma
      grad(3) = L*sqrt(2.0)*gamma
  

      IERR = 0

      return
    end subroutine EV_GRAD_F
!
! =============================================================================
!
!                     Computation of equality constraints
!
! =============================================================================
!
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X, M,i
      double precision G(M), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR

      real*8::u,v,sigma(3),theta,pu,pv,L,E,P

      real*8::tensile_sigma1_max,tensile_sigma2_max,tensile_sigma3_max

      real*8::comp_sigma1_max,comp_sigma2_max,comp_sigma3_max

      real*8::max_u_disp,max_v_disp

      

      !pi=4.0*atan(1.0) ! constant for later use (visible globally)
      
      P=dat(5)
      theta=dat(4)
      L=dat(1)
      E=dat(2)

      tensile_sigma1_max=dat(6)      
      tensile_sigma2_max=dat(7)
      tensile_sigma3_max=dat(8)

      comp_sigma1_max=dat(9)      
      comp_sigma2_max=dat(10)
      comp_sigma3_max=dat(11)

      max_u_disp=dat(12)
      max_v_disp=dat(13)
      

      pu=P*cos(theta)
      pv=P*sin(theta)

      u=(L/E)*(x(1)*pu + 2*sqrt(2.0)*x(2)*pu + x(3)*pu + x(3)*pv - x(1)*pv)/(x(1)*x(2) + sqrt(2.0)*x(1)*x(3) + x(2)*x(3))

      v=(L/E)*(-x(1)*pu + x(3)*pu + x(1)*pv + x(3)*pv)/(x(1)*x(2) + sqrt(2.0)*x(1)*x(3) + x(2)*x(3))

!!$if (loadcase.eq.2) then

      sigma(1)=-1.0*(sqrt(2.0)*x(2)*pu + x(3)*pu  +x(3)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))

!!$else

!!$   sigma(1)=(sqrt(2.0)*x(2)*pu + x(3)*pu  +x(3)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))

!!$ end if

      sigma(2)= (-(x(1)-x(3))*pu+(x(1)+x(3))*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))

      sigma(3)=(-sqrt(2.0)*x(2)*pu -x(1)*pu  +x(1)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))


      ! stress constraints (normalized)

      G(1) = (sigma(1) - tensile_sigma1_max)/tensile_sigma1_max    !tensile 1
      G(4) = (-1.0*sigma(1)/comp_sigma1_max) -1.0     !compressive 1

      G(2) = (sigma(2) - tensile_sigma2_max)/tensile_sigma2_max   !tensile 2
      G(5) = (-1.0*sigma(2)/comp_sigma2_max) -1.0  !compressive 2

      G(3) = (sigma(3) - tensile_sigma3_max)/tensile_sigma3_max    ! tensile 3
      G(6) = (-1.0*sigma(3) / comp_sigma3_max) -1.0 !compressive 3

      ! displacement constraints (normalized)

      G(7) = (-1.0*u -max_u_disp)/max_u_disp
      G(8) = (v -max_v_disp)/max_v_disp
      print*,''
      write(*,'(4x,a)') '>>Normalized Constraint Values:'
      do i=1,8
        write(*,'(E13.2)'),g(i)
      end do
      print*,''
      IERR = 0

      return
    end subroutine EV_G
!
! =============================================================================
!
!                Computation of Jacobian of equality constraints
!
! =============================================================================
!
    subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,IDAT, DAT, IERR)
      implicit none
      integer TASK, N, NEW_X, M, NZ
      double precision X(N), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*),cgrad(m,n)
      integer IDAT(*)
      integer IERR

      real*8::u,v,sigma(3),theta,pu,pv,L,E,P,pi


      real*8::tensile_sigma1_max,tensile_sigma2_max,tensile_sigma3_max

      real*8::comp_sigma1_max,comp_sigma2_max,comp_sigma3_max

      real*8::max_u_disp,max_v_disp


      pi=4.0*atan(1.0) ! constant for later use (visible globally)

      P=dat(5)
      theta=dat(4)
      L=dat(1)
      E=dat(2)

      tensile_sigma1_max=dat(6)      
      tensile_sigma2_max=dat(7)
      tensile_sigma3_max=dat(8)

      comp_sigma1_max=dat(9)      
      comp_sigma2_max=dat(10)
      comp_sigma3_max=dat(11)

      max_u_disp=dat(12)
      max_v_disp=dat(13)

      pu=P*cos(theta)
      pv=P*sin(theta)

      if( TASK.eq.0 ) then 
         !
         !     structure of Jacobian:
         !     


         ACON(1) = 1
         AVAR(1) = 1

         ACON(2) = 1
         AVAR(2) = 2

         ACON(3) = 1
         AVAR(3) = 3



         ACON(4) = 2
         AVAR(4) = 1

         ACON(5) = 2
         AVAR(5) = 2

         ACON(6) = 2
         AVAR(6) = 3



         ACON(7) = 3
         AVAR(7) = 1

         ACON(8) = 3
         AVAR(8) = 2

         ACON(9) = 3
         AVAR(9) = 3



         ACON(10) = 4
         AVAR(10) = 1

         ACON(11) = 4
         AVAR(11) = 2

         ACON(12) = 4
         AVAR(12) = 3



         ACON(13) = 5
         AVAR(13) = 1

         ACON(14) = 5
         AVAR(14) = 2

         ACON(15) = 5
         AVAR(15) = 3



         ACON(16) = 6
         AVAR(16) = 1

         ACON(17) = 6
         AVAR(17) = 2

         ACON(18) = 6
         AVAR(18) = 3



         ACON(19) = 7
         AVAR(19) = 1

         ACON(20) = 7
         AVAR(20) = 2

         ACON(21) = 7
         AVAR(21) = 3

         ACON(22) = 8
         AVAR(22) = 1

         ACON(23) = 8
         AVAR(23) = 2

         ACON(24) = 8
         AVAR(24) = 3

      else

         !---- GRADIENT OF CONSTRAINTS


         cgrad(:,:)=0.0

         ! Tensile Stress 1

         cgrad(1,1)=((x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(1,2)=((x(1) + x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*pu)/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

         cgrad(1,3)=((x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu + pv)/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))    

         ! Tensile Stress on 2          

         cgrad(2,1)=((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(3)))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu - pv)/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

         cgrad(2,2)=((x(1) + x(3))*(pu*(x(1) - x(3)) - pv*(x(1) + x(3))))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(2,3)=(pu + pv)/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) + ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(1)))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         ! Tensile Stress on 3          

         cgrad(3,1)=((x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu - pv)/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

         cgrad(3,2)=((x(1) + x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*pu)/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

         cgrad(3,3)=((x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         ! Compressive Stress on 1 

         cgrad(4,1)=-((x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(4,2)=(2.0**(1.0/2.0)*pu)/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(1) + x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(4,3)=(pu + pv)/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         ! Compressive Stress on 2 

         cgrad(5,1)=(pu - pv)/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(3)))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(5,2)=-((x(1) + x(3))*(pu*(x(1) - x(3)) - pv*(x(1) + x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(5,3)=- (pu + pv)/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(1)))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         ! Compressive Stress on 3 
         cgrad(6,1)=(pu - pv)/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(6,2)=(2.0**(1.0/2.0)*pu)/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(1) + x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(6,3)=-((x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         !X-Displacement
         cgrad(7,1) = ((1.0/max_u_disp)*L*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - ((1.0/max_u_disp)*L*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

         cgrad(7,2) = ((1.0/max_u_disp)*L*(x(1) + x(3))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (400.0*2.0**(1.0/2.0)*L*pu)/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

         cgrad(7,3) = ((1.0/max_u_disp)*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - ((1.0/max_u_disp)*L*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))


         ! Y-displacement

         cgrad(8,1)= -((1.0/max_v_disp)*L*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((1.0/max_v_disp)*L*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(8,2)= -((1.0/max_v_disp)*L*(x(1) + x(3))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

         cgrad(8,3)= ((1.0/max_v_disp)*L*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((1.0/max_v_disp)*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)




         A(1)=cgrad(1,1)
         A(2)=cgrad(1,2)
         A(3)=cgrad(1,3)

         A(4)=cgrad(2,1)
         A(5)=cgrad(2,2)
         A(6)=cgrad(2,3)

         A(7)=cgrad(3,1)
         A(8)=cgrad(3,2)
         A(9)=cgrad(3,3)

         A(10)=cgrad(4,1)
         A(11)=cgrad(4,2)
         A(12)=cgrad(4,3)

         A(13)=cgrad(5,1)
         A(14)=cgrad(5,2)
         A(15)=cgrad(5,3)

         A(16)=cgrad(6,1)
         A(17)=cgrad(6,2)
         A(18)=cgrad(6,3)

         A(19)=cgrad(7,1)
         A(20)=cgrad(7,2)
         A(21)=cgrad(7,3)

         A(22)=cgrad(8,1)
         A(23)=cgrad(8,2)
         A(24)=cgrad(8,3) 

      end if

      IERR = 0
      return
    end subroutine EV_JAC_G

!
! =============================================================================
!
!                Computation of Hessian of Lagrangian
!
! =============================================================================
!
      subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
      implicit none
      integer TASK, N, NEW_X, M, NEW_LAM, NNZH, i, ir,j
      double precision X(N), OBJFACT, LAM(M), HESS(NNZH),OBJHESS(NNZH),CONHESS(M,NNZH)

      integer IRNH(NNZH), ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      double precision :: hesstmp

      if( TASK.eq.0 ) then
!
!     structure of sparse Hessian (lower triangle):
!

         IRNH(1) = 1
         ICNH(1) = 1

         IRNH(2) = 2
         ICNH(2) = 2

         IRNH(3) = 3
         ICNH(3) = 3

         IRNH(4) = 2
         ICNH(4) = 1

         IRNH(5) = 3
         ICNH(5) = 2

         IRNH(6) = 3
         ICNH(6) = 1
        
      else
!!$!
!!$!     calculate Hessian:
!!$!
!!$
!!$         !Objective function
!!$
!!$         objhess(1)=0.0
!!$         objhess(2)=0.0
!!$         objhess(3)=1.0
!!$         
!!$         
!!$         ! first constraint
!!$         
!!$         conhess(1,1)= (12.0*BM*fs)/((x(1)**3)*(x(2)**2)*sigma_allow)
!!$         conhess(1,2)= (36.0*BM*fs)/(x(1)*(x(2)**4)*sigma_allow)
!!$         conhess(1,3)= (12.0*BM*Fs)/((x(1)**2)*(x(2)**3)*sigma_allow)
!!$         
!!$         ! Second constraint
!!$         
!!$         conhess(2,1)=(3.0*V*fs)/((x(1)**3)*x(2)*tau_allow)
!!$         conhess(2,2)=(3.0*V*fs)/(x(1)*(x(2)**3)*tau_allow)
!!$         conhess(2,3)=(3.0*V*fs)/(2.0*(x(1)**2)*(x(2)**2)*tau_allow)
!!$         
!!$         ! Third Constaint
!!$
!!$         conhess(3,1)=  x(2)/(x(1)**3)
!!$         conhess(3,2)=  0.0
!!$         conhess(3,3)= -1.0/(2.0*x(1)**2)
!!$
!!$         
!!$         ! Assemble
!!$         
!!$         HESS(:)=0.0
!!$         do i=1,NNZH
!!$            hesstmp=0.0
!!$            do j=1,m
!!$               hesstmp=hesstmp+lam(j)*conhess(j,i)
!!$            end do
!!$            hess(i)=hesstmp+objhess(i)
!!$         end do
         
      IERR = 0
      
   endif
   return
 end subroutine EV_HESS
 !
! =============================================================================
 !
 !                   Callback method called once per iteration
 !
 ! =============================================================================
!
      subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,DAT, ISTOP)
      implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      double precision DAT(*)
      integer IDAT(*)
      integer ISTOP

      if (ITER_COUNT .eq.0) then
         write(*,*) 
         write(*,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
      end if

      write(*,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU
      if (ITER_COUNT .gt. 1 .and. DNORM.le.1D-10) ISTOP = 1

      return
    end subroutine ITER_CB
