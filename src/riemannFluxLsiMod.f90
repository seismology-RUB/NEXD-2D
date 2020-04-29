!-----------------------------------------------------------------------
!   Copyright 2014-2020 Thomas Möller (Ruhr-Universität Bochum, GER)
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
module riemannFluxLsiMod
    use constantsMod
    use slipInterfaceMod
    use parameterMod

    implicit none

    contains

    subroutine riemannExpansionCoefficients(zin_vp, zout_vp, zin_vs, zout_vs, dV, dSigma, gamma)
    !This subroutine calculates the expansioncoefficients for the Riemanflux using the parameters
    !of both the element and its neighbour for the calculation.
        !input
        real(kind=custom_real), intent(in) :: zin_vp, zout_vp, zin_vs, zout_vs  !Impedances for the element and neighbour
        real(kind=custom_real), dimension(:,:), intent(in) :: dV                !Solution in the neighbouring element for the face that makes up the Slip interface. Dimension is (5,5)
        real(kind=custom_real), dimension(:,:), intent(in) :: dSigma            !Solution inside of the element for the face that makes up the Slip interface. Dimension is (5,5)
        !output
        real(kind=custom_real), dimension(:,:), intent(out) :: gamma
        !local
        integer :: i

        do i = 1, NGLL
            !P-SV
            gamma(i,1) = calcNonLSIGamma(dV(i,2), dSigma(i,3), zin_vs, zout_vs)
            gamma(i,2) = calcNonLSIGamma(dV(i,1), dSigma(i,1), zin_vp, zout_vp)
        enddo
    end subroutine

    subroutine riemannExpansionCoefficientsAddLSI(lsipar, zin_vp, zout_vp, zin_vs, zout_vs, S_n, S_t, gamma)
        !Adds the Slip-Interface influence to the expansion-coefficient.
        !input
        type(lsi_parameter) :: lsipar
        real(kind=custom_real), intent(in) :: zin_vp, zin_vs            !Impedance of the interior element
        real(kind=custom_real), intent(in) :: zout_vp, zout_vs          !Impedance of the neighbouring element
        real(kind=custom_real), dimension(:), intent(in) :: S_n         !Jump function for a normal jump
        real(kind=custom_real), dimension(:), intent(in) :: S_t         !Jump function for a tangential jump
        !in/out
        real(kind=custom_real), dimension(:,:), intent(inout) :: gamma  !Expansioncoefficients for P-SV
        !local
        integer :: i

        do i = 1, NGLL
            !P-SV
            if (lsipar%normal) then
                gamma(i,2) =  gamma(i,2) - calcLSIGamma(zin_vp, zout_vp, S_n(i))
            end if
            if (lsipar%tangential) then
                gamma(i,1) = gamma(i,1) - calcLSIGamma(zin_vs, zout_vs, S_t(i))
            end if
        enddo
    end subroutine

    function calcNonLSIGamma(dV, dSigma, zin, zout) result(gamma)
        !This function calculates the Expansion coefficients for the right going LSI Riemann Fluxes,
        !without the Slip-Interface influence!
        !Input
        real(kind=custom_real), intent(in) :: dv    !Velocity difference over the Interface
        real(kind=custom_real), intent(in) :: dSigma !Stress component
        real(kind=custom_real), intent(in) :: zin   !Impedance of the interior Element
        real(kind=custom_real), intent(in) :: zout  !Impedance of the neighbouring element
        !output
        real(kind=custom_real) :: gamma             !expansion coefficient for the Riemann-Flux

        gamma = 1/(zin + zout) * (dSigma - zout * dV)
    end function calcNonLSIGamma

    function calcLSIGamma(zin, zout, S) result(gamma)
        !This function calculates the Slip-Interface influence for the expansion coefficients.
        !Input
        real(kind=custom_real), intent(in) :: zin   !Impedance of the interior Element
        real(kind=custom_real), intent(in) :: zout  !Impedance of the neighbouring element
        real(kind=custom_real), intent(in) :: S     !Jump-Function
        !output
        real(kind=custom_real) :: gamma             !expansion coefficient for the Riemann-Flux

        gamma = 2 * zout/(zin + zout) * S
    end function calcLSIGamma

end module riemannFluxLsiMod
