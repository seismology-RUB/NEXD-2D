! -*- f90 -*-
  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8 !it is a know issue, that double precision currently doesn't work


  integer, parameter :: CUSTOM_REAL = SIZE_REAL
  integer, parameter :: ORDER = 4
  integer, parameter :: NGLL= ORDER +1 !5
  integer, parameter :: Np=(ORDER+1)*(ORDER+2)/2
  integer, parameter :: NpF=ORDER+1
  real(kind=CUSTOM_REAL), parameter :: PI = 3.141592653589793 
  real(kind=CUSTOM_REAL), parameter :: EPS = 1.0e-5 ! 1e-5 
  real(kind=CUSTOM_REAL), dimension(15), parameter :: balpha=(/0.0000, 0.0000, 1.4152, 0.1001,&
       & 0.2751, 0.9800, 1.0999, 1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258/)

! attenuation
  integer, parameter :: nMB = 3 ! Number of SLS

! maximum number of neighbors per element
  integer, parameter :: max_neighbor = 5 !5
  integer, parameter :: nsize=10 ! elements share one node 10
