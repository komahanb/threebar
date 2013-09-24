Module dimthree
  implicit none
  
  integer::maxdim,maxcon,maxdat
  parameter (maxdim=10,maxcon=20,maxdat=50)
  !  integer:: loadcase
  
  real*8 :: x(3) ! Design variable vector
  !  real*8 :: dat(maxdat) ! Constants and others for the constraint and objective function
  
  real*8 :: pi
  
  real*8 :: P,theta,L(3),E,gamma,pu,pv
  
  ! Objective and constaint function values
  real*8::objF,G(maxcon)
  
  ! Their gradients
  real*8::grad(maxdim),cgrad(maxcon,maxdim)
  
End Module dimthree
