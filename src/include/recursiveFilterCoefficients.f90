!-----------------------------------------------------------------------
!   Copyright 2013 Wolfgang Friederich (Ruhr-Universität Bochum, GER)
!
!   This file is part of NEXD 2D.
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful, but
!   WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with NEXD 2D. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------
!----------------------------------------------------------------------------
!  calculation of filter coefficients for recursive filtering
!----------------------------------------------------------------------------
module recursiveFilterCoefficients
    use constantsMod
    implicit none
    type recursive_filter_coeff
        character (len=3) :: typ
        double precision :: f0,f1,f2,g1,g2
    end type
!
contains
    !---------------------------------------------------------------------------
    !  coefficients for first order low pass filter
    !
    subroutine lowPassO1RecursiveFilterCoefficients(this,tc,dt)
        type (recursive_filter_coeff) :: this
        real :: tc,dt
        double precision :: tcn,eps
        !
        this%typ = 'lp1'
        tcn = tc/dt
        eps =2.*pi/tcn
        this%f2 = 0.d0
        this%g2 = 0.d0
        this%g1 = (2.d0-eps)/(2.d0+eps)
        this%f0 = eps/(2.d0+eps)
        this%f1 = this%f0
    end subroutine lowPassO1RecursiveFilterCoefficients
    !---------------------------------------------------------------------------
    !  coefficients for first order high pass filter
    !
    subroutine highPassO1RecursiveFilterCoefficients(this,tc,dt)
        type (recursive_filter_coeff) :: this
        real :: tc,dt
        double precision :: tcn,eps
        !
        this%typ = 'hp1'
        tcn = tc/dt
        eps =2.d0*pi/tcn
        this%f2 = 0.d0
        this%g2 = 0.d0
        this%g1 = (2.d0-eps)/(2.d0+eps)
        this%f0 = 2.d0/(2.d0+eps)
        this%f1 = -this%f0
    end subroutine highPassO1RecursiveFilterCoefficients
    !---------------------------------------------------------------------------
    !  coefficients for second order low pass filter
    !
    subroutine lowPassO2RecursiveFilterCoefficients(this,tc,h,dt)
        type (recursive_filter_coeff) :: this
        real :: tc,dt,h
        double precision :: tcn,eps,epsq,a,b,c
        !
        this%typ = 'lp2'
        tcn = tc/dt
        eps =2.d0*pi/tcn
        epsq = eps*eps
        a = 1.d0-eps*h+epsq/4.d0
        b = -2.d0+epsq/2.d0
        c = 1.d0+eps*h+epsq/4.d0
        this%g1 = -b/c
        this%g2 = -a/c
        this%f0 = epsq/4.d0/c
        this%f1 = 2.d0*this%f0
        this%f2 = this%f0
    end subroutine lowPassO2RecursiveFilterCoefficients
    !----------------------------------------------------------------------------
    !  coefficients for second order high pass filter
    !
    subroutine highPassO2RecursiveFilterCoefficients(this,tc,h,dt)
        type (recursive_filter_coeff) :: this
        real :: tc,dt,h
        double precision :: tcn,eps,epsq,a,b,c
        !
        this%typ = 'hp2'
        tcn = tc/dt
        eps =2.d0*pi/tcn
        epsq = eps*eps
        a = 1.d0-eps*h+epsq/4.d0
        b = -2.d0+epsq/2.d0
        c = 1.d0+eps*h+epsq/4.d0
        this%g1 = -b/c
        this%g2 = -a/c
        this%f0 = 1.d0/c
        this%f1 = -2.d0*this%f0
        this%f2 = this%f0
    end subroutine highPassO2RecursiveFilterCoefficients
    !----------------------------------------------------------------------------
    !  coefficients for second order band pass filter
    !
    subroutine bandPassO2RecursiveFilterCoefficients(this,tc,h,dt)
        type (recursive_filter_coeff) :: this
        real :: tc,dt,h
        double precision :: tcn,eps,epsq,a,b,c
        !
        this%typ = 'bp2'
        tcn = tc/dt
        eps =2.d0*pi/tcn
        epsq = eps*eps
        a = 1.d0-eps*h+epsq/4.d0
        b = -2.d0+epsq/2.d0
        c = 1.d0+eps*h+epsq/4.d0
        this%g1 = -b/c
        this%g2 = -a/c
        this%f0 = eps/2.d0/c
        this%f1 = 0.d0
        this%f2 = -this%f0
    end subroutine bandPassO2RecursiveFilterCoefficients
!
end module recursiveFilterCoefficients
