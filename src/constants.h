! -*- f90 -*-
  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8
  integer, parameter :: max_length_string = 400

  integer, parameter :: CUSTOM_REAL = SIZE_REAL
  integer, parameter :: ORDER = 4
  integer, parameter :: NGLL= ORDER +1 !5
  integer, parameter :: Np=(ORDER+1)*(ORDER+2)/2
  integer, parameter :: NpF=ORDER+1
  integer, parameter :: nsurface = 3
  real(kind=CUSTOM_REAL), parameter :: PI = 3.141592653589793
  real(kind=CUSTOM_REAL), parameter :: EPS = 1.0e-5 ! 1e-5
  real(kind=CUSTOM_REAL), dimension(15), parameter :: balpha=(/0.0000, 0.0000, 1.4152, 0.1001,&
       & 0.2751, 0.9800, 1.0999, 1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258/)

! usefulls:
  real(kind=custom_real), parameter :: zero = 0.0
  real(kind=custom_real), parameter :: one = 1.0
  real(kind=custom_real), parameter :: two = 2.0

! attenuation
  integer, parameter :: nMB = 3 ! Number of SLS

! maximum number of neighbors per element
  integer, parameter :: max_neighbor = 10 !5
  integer, parameter :: nsize=100 ! elements share one node 10

! Low storage Runge-Kutta coefficients
  real(kind=custom_real), dimension(5), parameter :: rk4a = (/0.0, &
        -567301805773.0/1357537059087.0, &
        -2404267990393.0/2016746695238.0, &
        -3550918686646.0/2091501179385.0, &
        -1275806237668.0/842570457699.0 /)
  real(kind=custom_real), dimension(5), parameter :: rk4b = (/1432997174477.0/9575080441755.0, &
         5161836677717.0/13612068292357.0, &
         1720146321549.0/2090206949498.0,  &
         3134564353537.0/4481467310338.0,  &
         2277821191437.0/14882151754819.0/)
  real(kind=custom_real), dimension(5), parameter :: rk4c = (/0.0,  &
         1432997174477.0/9575080441755.0, &
         2526269341429.0/6820363962896.0, &
         2006345519317.0/3224310063776.0, &
         2802321613138.0/2924317926251.0/)
