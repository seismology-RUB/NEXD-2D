!-----------------------------------------------------------------------
!   Copyright 2014-2018 Thomas Möller (Ruhr-Universität Bochum, GER)
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
!> \brief Handle slip interface properties characterized by eta and xi
!
module slipInterfaceMod
    use errorMessage
    use constantsMod
    use mpiMod

    implicit none

    type lsi_parameter
        !A Series of Switches to turn certain parts of the lsi on and off
        integer :: lsin
        logical :: lsi          !Calculate LSI influence
        logical :: normal       !Jump for P-Waves
        logical :: tangential   !Jump for SV-Waves
        logical :: cubitFile    !LSI file created by cubit/trelis is filled with data
    end type

    !Type that contains the specifications of the interfaces. Information from the "Interfaces" file is stored here.
    type :: interface_spec
        character(len=10) :: type                                   !Interface types: Choices are elastic
        real(kind=custom_real) :: e_mod                             !Youngs-Modulus
        real(kind=custom_real) :: h                                 !Width
        real(kind=custom_real) :: vp                                !P-wave velocity
        real(kind=custom_real) :: vs                                !S-Wave velocity
        real(kind=custom_real) :: rho                               !Density
        real(kind=custom_real) :: emod                              !Young's Modulus
        real(kind=custom_real) :: C_n                               !Relation between the Elastic modulus of the crack and the reference Modulus. (Normal to the crack)
        real(kind=custom_real) :: C_t                               !Relation between the Elastic modulus of the crack and the reference Modulus. (tangentical to the crack)
        real(kind=custom_real) :: eta_vp                            !Compliance (P-Wave)
        real(kind=custom_real) :: eta_vs                            !Compliance (S-Wave)
    end type

    !Type that contains the information neede to calculate the activity of each interface
    type :: lsiVar !initialized for each interface
        integer, dimension(2) :: elements                                     !Elements that connect over the lsip interface
        integer :: prop                                                       !Property index of the slip interface
        real(kind=custom_real), dimension(nsurface, NpF) :: S_t               !Storage variable for the Slipinterfaceactivity
        real(kind=custom_real), dimension(nsurface, NpF) :: rhsS_t            !right hand side of the evolution equation for S
        real(kind=custom_real), dimension(nsurface, NpF) :: resS_t            !Runge-Kutta residual (S)
        real(kind=custom_real), dimension(nsurface, NpF) :: S_n               !Storage variable for the Slipinterfaceactivity
        real(kind=custom_real), dimension(nsurface, NpF) :: rhsS_n            !Right-hand-side of the evolution equation for S
        real(kind=custom_real), dimension(nsurface, NpF) :: resS_n            !Runge-Kutta residual (S)
    end type lsiVar

    type elementToLSI
        logical :: isLSI                        !Boolean that is true if the element has at least on slip interface surface
        logical, dimension(3) :: lsisurface     !Boolean for each surface: true if the surface is an slip interface
        integer, dimension(3) :: lsiprop        !Index for the surfaces. 0 for no interface and > 0 for an interface
    end type elementToLSI

    type lsi_location
        ! This type contains values that define the type of the interface (property), and it's start and end point in the model.
        ! From these values the path is caluated that defines the interface.
        integer :: number
        integer :: property
        real(kind=custom_real) :: start_x
        real(kind=custom_real) :: start_z
        real(kind=custom_real) :: end_x
        real(kind=custom_real) :: end_z
    end type

    contains

    subroutine setup_lsi(nelem, etolsi, lsi, lsi_spec, lsipar, filename, myrank, errmsg)
        !Input
        integer :: myrank
        integer :: nelem
        !input/output
        type(lsi_parameter) :: lsipar
        type(error_message) :: errmsg
        type(lsiVar), dimension(:), allocatable :: lsi
        type(interface_spec), dimension(:), allocatable :: lsi_spec
        type(elementToLSI), dimension(:), allocatable :: etolsi
        character(len=*) :: filename
        !local
        character(len=15) :: myname = 'setup_lsi'

        call addTrace(errmsg,myname)

        !Read the specifications of all possible interface types
        call readSlipInterfaceSpecs(lsi_spec, 1, "data/interfaces" ,errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif

        !read element to slip interface mapping
        call readLSIDatabase(1, filename, nelem, lsi, lsipar, etolsi, myrank, errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    end subroutine

    subroutine readLSIlocation(lsi_loc, coord, set_pml, pml_delta, lu, filename, errmsg)
        !This subroutine is used to read in the file containing the information on start- and end-point of the interface.
        !input
        integer :: lu
        character(len=*) :: filename
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL) :: pml_delta
        real(kind=CUSTOM_REAL), dimension(:,:) :: coord
        logical :: set_pml
        !output
        type(lsi_location), dimension(:), allocatable :: lsi_loc
        !local
        character(len=27) :: myname = 'readLSILocation'
        real(kind = CUSTOM_REAL) :: pml_minX, pml_minZ, pml_maxX, pml_maxZ
        integer :: nlsi
        integer :: i
        integer :: maxprop                                          ! This is a dummy variable to read the maximum number
                                                                    ! of different slip interfaces are specified in the parfile
                                                                    ! This number is specfied in the property index, hence the name -> maximum property index

        call addTrace(errmsg,myname)

        if (set_pml) then
        !pml boundary
            pml_maxX = maxval(coord(1,:)) - pml_delta
            pml_minX = minval(coord(1,:)) + pml_delta
            pml_maxZ = maxval(coord(2,:)) - pml_delta
            pml_minZ = minval(coord(2,:)) + pml_delta
        end if

        maxprop = readFromLSIParfile("data/interfaces", myname, lu, errmsg)
        close(lu)

        nlsi = readFromLSIParfile(filename, myname, lu, errmsg)

        if (maxprop <= 0 ) then
            call add(errmsg,2,"The maximum number of slip interface types must be 1 or above.",myname, "data/interfaces")
        end if
        if (nlsi <= 0 ) then
            call add(errmsg,2,"The number of fractures must be 1 or above.",myname, filename)
        end if

        allocate (lsi_loc(nlsi))
        do i = 1, nlsi
            read(lu,*) lsi_loc(i)%number, lsi_loc(i)%property, lsi_loc(i)%start_x, lsi_loc(i)%start_z, lsi_loc(i)%end_x, lsi_loc(i)%end_z
            if (lsi_loc(i)%property > maxprop) then
                call add(errmsg,2,"The specified property index for this fracture is too large.",myname, filename)
            end if
            if (lsi_loc(i)%property <= 0) then
                call add(errmsg,2,"The specified property index for this fracture must be 1 or above.",myname, filename)
            end if
            call modelConflict(lsi_loc(i)%start_x, "start_x", minval(coord(1,:)), minval(coord(2,:)), maxval(coord(1,:)), maxval(coord(2,:)), myname, filename, errmsg)
            call modelConflict(lsi_loc(i)%start_z, "start_z", minval(coord(1,:)), minval(coord(2,:)), maxval(coord(1,:)), maxval(coord(2,:)), myname, filename, errmsg)
            call modelConflict(lsi_loc(i)%end_x, "end_x", minval(coord(1,:)), minval(coord(2,:)), maxval(coord(1,:)), maxval(coord(2,:)), myname, filename, errmsg)
            call modelConflict(lsi_loc(i)%end_z, "end_z", minval(coord(1,:)), minval(coord(2,:)), maxval(coord(1,:)), maxval(coord(2,:)), myname, filename, errmsg)
        enddo
        close(lu)
    end subroutine

!---------------------------------------------------------------------
!> \brief Read slip interface specifications from a file
!
    subroutine readSlipInterfaceSpecs(lsi_spec, lu, filename, errmsg)
        !input
        integer :: lu
        character(len=*) :: filename
        type (error_message) :: errmsg
        !output
        type(interface_spec), dimension(:), allocatable :: lsi_spec
        !local
        character(len=27) :: myname = 'readSlipInterfaceSpec'
        integer :: nslip
        integer :: i

        call addTrace(errmsg,myname)

        nslip = readFromLSIParfile(filename, myname, lu, errmsg)

        allocate (lsi_spec(nslip))

        do i = 1, nslip
            read(1,*) lsi_spec(i)%type,lsi_spec(i)%h,lsi_spec(i)%emod,lsi_spec(i)%C_n,lsi_spec(i)%C_t
            select case(trim(lsi_spec(i)%type))
                case ("elastic")
                    lsi_spec(i)%eta_vp = (2*lsi_spec(i)%C_n*lsi_spec(i)%h)/lsi_spec(i)%emod
                    lsi_spec(i)%eta_vs = (2*lsi_spec(i)%C_t*lsi_spec(i)%h)/lsi_spec(i)%emod
                case default
                    call add(errmsg,2,"No accepted interface-type entered.",myname, "data/interfaces")
                    call print(errmsg); return
            end select
        enddo
        close(lu)
    end subroutine readSlipInterfaceSpecs

    subroutine readLSIDatabase(lu, filename, nelem, lsi, lsipar, etolsi, myrank, errmsg)
        !Subroutine to read the lsi database files created in mesher!
        !Input
        type(error_message) :: errmsg
        character(len=*), intent(in) :: filename
        integer, intent(in) :: lu
        integer, intent(in):: nelem
        integer, intent(in) :: myrank
        !output
        type(lsiVar), dimension(:), allocatable :: lsi
        type(lsi_parameter) :: lsipar
        type(elementToLSI), dimension(:), allocatable :: etolsi
        !local
        character(len=15) :: myname = "readLSIDatabase"
        integer :: e, ios

        call addTrace(errmsg, myname)

        open(unit=lu,file=trim(filename), iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open file',myname, filename)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif

        read(lu,*) lsipar%lsin

        write (*, '(a51, i5, a2, i5, a17)') &
        "|                Number of slip interfaces in rank ", myrank, ": ", lsipar%lsin, "                |"

        allocate(lsi(lsipar%lsin))
        allocate(etolsi(nelem))
        do e = 1, nelem
            read(lu,*) etolsi(e)%islsi, etolsi(e)%lsisurface(1), etolsi(e)%lsisurface(2), etolsi(e)%lsisurface(3),&
                       etolsi(e)%lsiprop(1), etolsi(e)%lsiprop(2), etolsi(e)%lsiprop(3)
        enddo

    end subroutine readLSIDatabase

    function readFromLSIParfile(filename, myname, lu, errmsg) result(value)
        !This function is designed to read the gobal parameter from the slip interface parameter files "interfaces" and "fracs"
        !input
        type(error_message) :: errmsg
        character(len = *) :: filename
        character(len = *) :: myname
        integer :: lu
        !output
        integer :: value
        !local
        character(len=max_length_string) :: text
        integer :: ios

        open(unit=lu,file=trim(filename), iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open file!', myname, filename)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif
        text = 'none'
        do while(text /= 'BEGIN')
            read(lu,'(a)') text
        enddo
        read(lu,*) value
    end function
end module slipInterfaceMod
