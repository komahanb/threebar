program threebar		
  use dimthree
  integer,parameter:: N = 3  ! Number of Design Varibles
  integer,parameter:: M = 16 ! Number of constaints  
  
  call loaddata

  call objfn
  call consfn

  call gradobj
  call gradconst
  

!  print*,objF,G(1:8)
!  print*,grad(1:3)
  print*,cgrad(1:8,1:3)

!  call setup
!  call solve
!  call results

!Author : Komahan Boopathy, University of Dayton, OH, 11/22/2013
end program threebar
