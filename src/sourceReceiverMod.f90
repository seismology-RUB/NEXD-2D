!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2018 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2018 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2018 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
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
module sourceReceiverMod
    use constantsMod
    use parameterMod
    use nodesMod
    use triTrafoMod
    use vandermondeMod
    use matrixMod
    use geometricFactorsMod
    use errorMessage

    implicit none

    type :: srcVar
        integer :: nsrc                                                           !Total Number of Sources
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcxz => null()        !x-z-values of the source
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcrs => null()        !r-s-values of the source
        integer, dimension(:), pointer :: srcelem => null()                       !Element in which the source is located
        integer, dimension(:), pointer :: srci => null()                          !Source number
        integer, dimension(:), pointer :: srctype => null()                       !Type of source: 1 = Moment tensor, 0 = single force
        integer, dimension(:), pointer :: srcstf => null()                        !Type of source-time-function: 1 = gauss, 2 = ricker, 3 = sin^3, 4 = external
        character(len=80), dimension(:), pointer :: extwavelet => null()          !File with external source wavelet
        real(kind=CUSTOM_REAL), dimension(:), pointer :: srcf0 => null()          !Center frequency of sft
        real(kind=CUSTOM_REAL), dimension(:), pointer :: srcfactor => null()      !Factor of the sft
        real(kind=CUSTOM_REAL), dimension(:), pointer :: srcangle_force => null() !Angle of the force action in the media
        real(kind=custom_real), dimension(:), allocatable :: delay                !Factor that enacts a time delay on the source
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcM => null()         !Momenttensor of the source
    end type srcVar

    type :: recVar
        integer :: nrec                                                     !Total Number of receivers
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: recxz => null()  !x-z-values of the receiver
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: recrs => null()  !r-s-values of the receiver
        integer, dimension(:), pointer :: recelem => null()                 !Element in which the receiver is located
        integer, dimension(:), pointer :: reci => null()                    !Receiver number
        integer, dimension(:), pointer :: recnr => null()
    end type recVar

    contains

    subroutine initSource(par,this, coord, errmsg)
        !subroutine to read the source parameters from the corrosponding file.
        type(error_message) :: errmsg
        type (srcVar) :: this
        type (parameterVar) :: par
        real(kind=custom_real), dimension(:,:) :: coord
        !local variables
        character(len=80) :: filename
        character(len=10) :: myname = "initSource"

        integer :: ier, isrc, pos
        logical :: file_exists

        call addTrace(errmsg, myname)

        filename=trim('data/source')
        open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
        if (ier /= 0) then
            call add(errmsg,2,'Could not open file!',myname,filename)
            call print(errmsg)
            stop
        endif

        !Cycle to read the parameters form the source file. The order of appearence of the parameters is not important
        call readIntPar(this%nsrc, "nsrc", filename, 0, errmsg)
        allocate(this%srctype(this%nsrc))
        allocate(this%srcxz(2,this%nsrc))
        allocate(this%srcrs(2,this%nsrc))
        allocate(this%srcstf(this%nsrc))
        allocate(this%extwavelet(this%nsrc))
        allocate(this%srcf0(this%nsrc))
        allocate(this%srcfactor(this%nsrc))
        allocate(this%srcangle_force(this%nsrc))
        allocate(this%srcM(3,this%nsrc))
        allocate(this%srci(this%nsrc))
        allocate(this%srcelem(this%nsrc))
        allocate(this%delay(this%nsrc))

        if (par%log) then
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, i10, a30)")   "|                    Number of sources: ", this%nsrc, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        end if

        do isrc = 1, this%nsrc
            pos = setFilePosition("source", filename, isrc, errmsg)
            call readIntPar(this%srctype(isrc), "sourcetype", filename, pos, errmsg)
            if (this%srctype(isrc) == 0 .or. this%srctype(isrc) == 1) then
                continue
            else
                call add(errmsg, 2, "Parameter to select the source type is out of range. Select either 0 or 1.", myname, filename)
            end if
            call readIntPar(this%srcstf(isrc), "stf", filename, pos, errmsg)
            if( this%srcstf(isrc) < 1 .or. this%srcstf(isrc) > 4) then
                !this message needs to be adjusted if new wavelets are added to the Program
                call add(errmsg, 2, "Parameter to select the source time function is out of range. Select either 1, 2, 3 or 4.", myname, filename)
            end if
            if (this%srcstf(isrc) == 4) then
                call readStringPar(this%extwavelet(isrc), "extwavelet", filename, pos, errmsg)
                inquire(file=trim(this%extwavelet(isrc)), exist=file_exists)
                if (.not. file_exists) then
                    call add(errmsg, 2, "File does not exist!", myname, filename)
                end if
            end if
            call readFloatPar(this%srcxz(1, isrc), "xsource", filename, pos, errmsg)
            call readFloatPar(this%srcxz(2, isrc), "zsource", filename, pos, errmsg)
            call modelConflict(this%srcxz(1, isrc), "xsource", minval(coord(1,:)), minval(coord(2,:)), maxval(coord(1,:)), maxval(coord(2,:)), myname, filename, errmsg)
            call modelConflict(this%srcxz(2, isrc), "zsource", minval(coord(1,:)), minval(coord(2,:)), maxval(coord(1,:)), maxval(coord(2,:)), myname, filename, errmsg)
            call readFloatPar(this%srcf0(isrc), "f0", filename, pos, errmsg)
            call readFloatPar(this%srcfactor(isrc), "factor", filename, pos, errmsg)
            call readFloatPar(this%srcangle_force(isrc), "angle_force", filename, pos, errmsg)
            call readFloatPar(this%srcM(1, isrc), "Mxx", filename, pos, errmsg)
            call readFloatPar(this%srcM(2, isrc), "Mzz", filename, pos, errmsg)
            call readFloatPar(this%srcM(3, isrc), "Mxz", filename, pos, errmsg)
            call readFloatPar(this%delay(isrc), "delay", filename, pos, errmsg)
        enddo
        close(19)



        if (par%log) then
            do isrc = 1, this%nsrc
                write (*,"(a40, f10.2, a30)")&
               "|                X-value of the source: ", this%srcxz(1, isrc), "                             |"
                write (*,"(a40, f10.2, a30)")&
               "|                Z-value of the source: ", this%srcxz(2, isrc), "                             |"
                write (*,"(a40, es10.3, a30)")&
               "|                           time delay: ", this%delay(isrc), "                             |"
                write (*,"(a40, i10, a30)")&
               "|                       Type of source: ", this%srctype(isrc), "                             |"
                write (*,"(a40, i10, a30)")&
               "|           source-time-function (stf): ", this%srcstf(isrc), "                             |"
                if (this%srcstf(isrc) == 4) then
                    write (*,"(a40, a38, a2)")&
                    "|      file with external sourcewavlet: ", this%extwavelet(isrc), " |"
                end if
                write (*,"(a40, f10.1, a30)")&
               "|               Centerfrequency of stf: ", this%srcf0(isrc), "                             |"
                write (*,"(a40, f10.1, a30)")&
               "|                    Factor of the stf: ", this%srcfactor(isrc), "                             |"
                write (*,"(a40, f10.2, a30)")&
               "|            Angle of the force action: ", this%srcangle_force(isrc), "                             |"
                write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Mxx: ", this%srcM(1, isrc), "                             |"
                write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Mzz: ", this%srcM(2, isrc), "                             |"
                write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Mxz: ", this%srcM(3, isrc), "                             |"
                write (*, "(a80)") "|------------------------------------------------------------------------------|"
            enddo
        endif
    end subroutine initSource

    subroutine initReceiver(par, this, coord, errmsg)
        type(error_message) :: errmsg
        type (parameterVar) :: par
        type(recVar) :: this
        integer :: irec, ier
        real(kind=custom_real), dimension(:,:) :: coord
        character(len=256) :: filename, dummy
        character(len=12) :: myname = "initReceiver"

        call addTrace(errmsg, myname)

        filename=trim('data/stations')
        open(unit=19,file=trim(filename), status ='old', iostat = ier)
        if (ier /= 0) then
            call add(errmsg,2,'Could not open file!',myname, filename)
            call print(errmsg)
            stop
        endif

        call readIntPar(this%nrec, "nrec", filename, 0, errmsg)
        if (.level.errmsg ==2 ) then; call print(errmsg);stop; endif
        read(19,*) dummy !This line exists to skip the first line

        allocate(this%recxz(2,this%nrec))
        allocate(this%recrs(2,this%nrec))
        allocate(this%recelem(this%nrec))
        allocate(this%reci(this%nrec))
        allocate(this%recnr(this%nrec))
        do irec=1,this%nrec
            read(19,*) this%recnr(irec),this%recxz(1,irec),this%recxz(2,irec)
            call modelConflict(this%recxz(1, irec), "xrec", minval(coord(1,:)), minval(coord(2,:)), maxval(coord(1,:)), maxval(coord(2,:)), myname, filename, errmsg)
            call modelConflict(this%recxz(2, irec), "zrec", minval(coord(1,:)), minval(coord(2,:)), maxval(coord(1,:)), maxval(coord(2,:)), myname, filename, errmsg)
        end do
        close(19)
        this%recrs=0
        this%recelem=0
        this%reci=0

        if(par%log) then
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a40, i10, a30)")   "|                  Number of receivers: ", this%nrec, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            do irec=1,this%nrec
                write (*,"(a40, i10, a30)") &
               "|                             reciever: ", this%recnr(irec), "                             |"
                write (*,"(a40, f10.2, a30)")&
               "|              X-value of the reciever: ", this%recxz(1, irec), "                             |"
                write (*,"(a40, f10.2, a30)")&
               "|              Z-value of the reciever: ", this%recxz(2, irec), "                             |"
                write (*, "(a80)") "|------------------------------------------------------------------------------|"

            enddo
        endif
    end subroutine initReceiver

    subroutine findSource(par,this,vx,vz,nelem,ibool,coord,elem,dr,ds,pml,errmsg)
        type(parameterVar) :: par
        type(srcVar) :: this
        type(error_message) :: errmsg
        integer, dimension(:) :: pml
        integer, dimension(:,:), pointer, intent(in) :: elem
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: coord
        integer, pointer, dimension(:,:), intent(in) :: ibool
        real(kind=CUSTOM_REAL), dimension(:), pointer, intent(in) :: vx,vz
        integer :: nelem
        real(kind=CUSTOM_REAL), dimension(Np) :: x,z,r,s
        logical, dimension(nelem,this%nsrc) :: checkelem
        integer :: isrc, i, ie, j, iglob, besave
        integer, dimension(NP) :: iv
        real(kind=CUSTOM_REAL) :: epsilon
        real(kind=CUSTOM_REAL) :: src_r, src_s
        real(kind=CUSTOM_REAL), dimension(1) :: src_rr, src_ss
        real(kind=CUSTOM_REAL), dimension(2) :: xztemp
        integer :: num_iter=4
        real(kind=CUSTOM_REAL), dimension(2) :: xz
        real(kind=CUSTOM_REAL), dimension(Np) :: rx,rz,sx,sz, jacobian
        real(kind=CUSTOM_REAL) :: dx,dz
        real(kind=CUSTOM_REAL), dimension(Np) :: srctemp, srcInt
        real(kind=CUSTOM_REAL), dimension(Np,NP) :: vdm, VdmTinv
        real(kind=CUSTOM_REAL), dimension(:,:) :: dr,ds
        real(kind=CUSTOM_REAL) :: drdx,drdz,dsdx,dsdz
        real(kind=CUSTOM_REAL)  :: ddr,dds

        character(len=10) :: myname = "findSource"
        character(len=256) :: errstr


        call addTrace(errmsg, myname)

        checkelem=.true.
        ! get local element
        call nodes(balpha(NGLL),x,z)
        call xyToRs(x,z,r,s)

        call vdm2d(vdm,r,s)
        vdmTinv=vdm
        vdmTinv=transpose(vdmTinv)
        call invert(vdmTinv, errmsg)

        src_r = 0.
        src_s = 0.

        ! do initial guess for source points
        isrc=1

        besave=1
        do while (isrc <= this%nsrc)
            epsilon=1e7
            do ie=1,nelem
                do i=1,Np
                    iglob=ibool(i,ie)
                    if ((sqrt((vx(iglob)-this%srcxz(1,isrc))**2+(vz(iglob)-this%srcxz(2,isrc))**2) < epsilon).and.checkelem(ie,isrc)) then
                        epsilon= sqrt((vx(iglob)-this%srcxz(1,isrc))**2+(vz(iglob)-this%srcxz(2,isrc))**2)
                        this%srci=iglob !initial guess
                        src_r=r(i)
                        src_s=s(i)
                        this%srcelem(isrc)=ie
                    end if
                end do
            end do

            xztemp(1) = this%srcxz(1,isrc)
            xztemp(2) = this%srcxz(2,isrc)
            ie=this%srcelem(isrc)
            iv = ibool(:,ie)
            call geometricFactors2d(rx,sx,rz,sz,jacobian,vx(iv),vz(iv),dr,ds, errmsg)

            do i=1, num_iter
                xz(1) = 0.5 * ( -(src_r+src_s) * coord(1,elem(1,ie)) + (1.0+src_r) * coord(1,elem(2,ie)) +&
                        (1+src_s) * coord(1,elem(3,ie)))
                xz(2) = 0.5 * ( -(src_r+src_s) * coord(2,elem(1,ie)) + (1.0+src_r) * coord(2,elem(2,ie)) +&
                        (1+src_s) * coord(2,elem(3,ie)))

                src_rr(1)=src_r
                src_ss(1)=src_s

                call vdm2D(srcTemp(:),src_rr,src_ss)
                srcInt(:)=matmul(vdmTinv,srcTemp(:))

                drdx=0
                drdz=0
                dsdx=0
                dsdz=0
                do j=1,NP
                    drdx=drdx+srcInt(j)*rx(j)
                    drdz=drdz+srcInt(j)*rz(j)
                    dsdx=dsdx+srcInt(j)*sx(j)
                    dsdz=dsdz+srcInt(j)*sz(j)
                end do

                dx = (xztemp(1)-xz(1))
                dz = (xztemp(2)-xz(2))

                ddr= drdx*dx+drdz*dz
                dds= dsdx*dx+dsdz*dz

                src_r = src_r + ddr
                src_s = src_s + dds

                if (i==num_iter) then
                    if (src_r < -1.01) then
                        isrc=isrc
                        checkelem(ie,isrc) = .false.
                        besave=besave+1
                    else if (src_r > 1.01) then
                        isrc=isrc
                        checkelem(ie,isrc) = .false.
                        besave=besave+1
                    else if (src_s < -1.01) then
                        isrc=isrc
                        checkelem(ie,isrc) = .false.
                        besave=besave+1
                    else if (src_s > 1.01) then
                        isrc=isrc
                        checkelem(ie,isrc) = .false.
                        besave=besave+1
                    else
                        this%srcrs(1,isrc) = src_r
                        this%srcrs(2,isrc) = src_s
                        if (pml(ie) == 1 .and. par%set_pml) then
                            write(errstr,*) "Source ", isrc, " is located inside the PML layer. Check the location of this source."
                            call add(errmsg, 2, errstr, myname, "data/source")
                        end if
                        if (par%log) then
                            write(*,'(a40, i4, a36)') "|                          Find source: ",isrc, "                                   |"
                            write(*,'(a41, f8.2, a2, f8.2, a21)') &
                                "|                   global coordinates: (", xz(1), ", ",  xz(2), " )                  |"
                            write(*,'(a40, i4, a36)') "|               Final guess for source: ",isrc, "                                   |"
                            write(*,'(a40, i6, a34)') "|                       global element: ", this%srcelem(isrc), "                                 |"
                            write(*,'(a41, f8.2, a2, f8.2, a21)') &
                                "|                    local coordinates: (", this%srcrs(1,isrc),", ", this%srcrs(2,isrc), " )                  |"
                        end if
                        isrc=isrc+1
                        besave=1
                    end if
                    if (besave > 10) then !10 versuche sonst abbruch
                        write(errstr,*) "ERROR, could not find source", isrc,"at global coordinates", xz(1),xz(2)
                        call add(errmsg, 2, errstr, myname, "data/source")
                        call print(errmsg)
                        stop
                    end if
                end if
            end do
        end do !isrc
    end subroutine findSource

    subroutine findReceiver(par,this,vx,vz,nelem,ibool,coord,elem,dr,ds,pml,errmsg)
        type(parameterVar) :: par
        type(recVar), intent(inout) :: this
        type(error_message) :: errmsg
        integer, dimension(:) :: pml
        integer, dimension(:,:), pointer, intent(in) :: elem
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: coord
        integer, pointer, dimension(:,:), intent(in) :: ibool
        real(kind=CUSTOM_REAL), dimension(:), pointer, intent(in) :: vx,vz
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: dr,ds
        integer :: nelem
        real(kind=CUSTOM_REAL), dimension(Np) :: x,z,r,s

        logical, dimension(nelem,this%nrec) :: checkelem
        integer :: irec, i, ie, j, iglob, besave
        integer, dimension(NP) :: iv
        real(kind=CUSTOM_REAL) :: epsilon
        real(kind=CUSTOM_REAL) :: rec_r, rec_s
        real(kind=CUSTOM_REAL), dimension(1) :: rec_rr, rec_ss
        real(kind=CUSTOM_REAL), dimension(2) :: xztemp
        integer :: num_iter=4
        real(kind=CUSTOM_REAL), dimension(2) :: xz
        real(kind=CUSTOM_REAL), dimension(Np) :: rx,rz,sx,sz, jacobian
        real(kind=CUSTOM_REAL) :: dx,dz
        real(kind=CUSTOM_REAL), dimension(Np) :: rectemp, recInt
        real(kind=CUSTOM_REAL), dimension(Np,NP) :: vdm, VdmTinv

        real(kind=CUSTOM_REAL) :: drdx,drdz,dsdx,dsdz
        real(kind=CUSTOM_REAL)  :: ddr,dds

        character(len=12) :: myname = "findReceiver"
        character(len=256) :: errstr

        call addTrace(errmsg, myname)

        checkelem=.true.
        ! get local element
        call nodes(balpha(NGLL),x,z)
        call xyToRs(x,z,r,s)
        call vdm2d(vdm,r,s)
        vdmTinv=vdm
        vdmTinv=transpose(vdmTinv)
        call invert(vdmTinv, errmsg)

        rec_r = 0.
        rec_s = 0.

        ! do initial guess for receiver points
        irec=1

        besave=1
        do while (irec <= this%nrec)
            epsilon=1e7
            do ie=1,nelem
                do i=1,Np
                    iglob=ibool(i,ie)
                    if ((sqrt((vx(iglob)-this%recxz(1,irec))**2+(vz(iglob)-this%recxz(2,irec))**2) < epsilon).and.checkelem(ie,irec)) then
                        epsilon= sqrt((vx(iglob)-this%recxz(1,irec))**2+(vz(iglob)-this%recxz(2,irec))**2)
                        this%reci=iglob !initial guess
                        rec_r=r(i)
                        rec_s=s(i)
                        this%recelem(irec)=ie
                    end if
                end do
            end do
            xztemp(1) = this%recxz(1,irec)
            xztemp(2) = this%recxz(2,irec)
            ie=this%recelem(irec)
            iv = ibool(:,ie)
            call geometricFactors2d(rx,sx,rz,sz,jacobian,vx(iv),vz(iv),dr,ds, errmsg)

            do i=1, num_iter
                xz(1) = 0.5 * ( -(rec_r+rec_s) * coord(1,elem(1,ie)) + (1.0+rec_r) * coord(1,elem(2,ie)) +&
                        (1+rec_s) * coord(1,elem(3,ie)) )
                xz(2) = 0.5 * ( -(rec_r+rec_s) * coord(2,elem(1,ie)) + (1.0+rec_r) * coord(2,elem(2,ie)) +&
                        (1+rec_s) * coord(2,elem(3,ie)) )

                rec_rr(1)=rec_r
                rec_ss(1)=rec_s

                call vdm2D(recTemp(:),rec_rr,rec_ss)
                recInt(:)=matmul(vdmTinv,recTemp(:))

                drdx=0
                drdz=0
                dsdx=0
                dsdz=0
                do j=1,NP
                    drdx=drdx+recInt(j)*rx(j)
                    drdz=drdz+recInt(j)*rz(j)
                    dsdx=dsdx+recInt(j)*sx(j)
                    dsdz=dsdz+recInt(j)*sz(j)
                end do

                dx = (xztemp(1)-xz(1))
                dz = (xztemp(2)-xz(2))

                ddr= drdx*dx+drdz*dz
                dds= dsdx*dx+dsdz*dz

                rec_r = rec_r + ddr
                rec_s = rec_s + dds

                if (i==num_iter) then
                    if (rec_r < -1.01) then
                        irec=irec
                        checkelem(ie,irec) = .false.
                        besave=besave+1
                    else if (rec_r > 1.01) then
                        irec=irec
                        checkelem(ie,irec) = .false.
                        besave=besave+1
                    else if (rec_s < -1.01) then
                        irec=irec
                        checkelem(ie,irec) = .false.
                        besave=besave+1
                    else if (rec_s > 1.01) then
                        irec=irec
                        checkelem(ie,irec) = .false.
                        besave=besave+1
                    else
                        this%recrs(1,irec) = rec_r
                        this%recrs(2,irec) = rec_s
                        if (pml(ie) == 1 .and. par%set_pml) then
                            write(errstr,*) "Station ", irec, " is located inside the PML layer. Check the location of this station."
                            call add(errmsg, 1, errstr, myname)
                        end if
                        if (par%log) then
                            write(*,'(a40, i4, a36)') "|                        Find receiver: ",irec, "                                   |"
                            write(*,'(a41, f8.2, a2, f8.2, a21)') &
                                "|                   global coordinates: (", xz(1), ", ",  xz(2), " )                  |"
                            write(*,'(a40, i4, a36)') "|             Final guess for receiver: ",irec, "                                   |"
                            write(*,'(a40, i6, a34)') "|                       global element: ", this%recelem(irec), "                                 |"
                            write(*,'(a41, f8.2, a2, f8.2, a21)') &
                                "|                    local coordinates: (", this%recrs(1,irec),", ", this%recrs(2,irec), " )                  |"
                        end if
                        irec=irec+1
                        besave=1
                    end if
                    if (besave > 10) then !10 versuche sonst abbruch
                        write(errstr,*) "ERROR, could not find receiver", irec,"at global coordinates", xz(1),xz(2)
                        call add(errmsg, 2, errstr, myname, "data/stations")
                        call print(errmsg)
                        stop
                    end if
                end if
            end do
        end do !irec
    end subroutine findReceiver

!
    !"--------------------------------------------------------------------------------------"
    ! change src element
    subroutine changeSrcElement(this,i,element)
        type(srcVar) :: this
        integer :: element,i
        this%srcelem(i) = element
    end subroutine changeSrcElement

    !"--------------------------------------------------------------------------------------"
    ! change rec element
    subroutine changeRecElement(this,i,element)
        type(recVar) :: this
        integer :: element,i
        this%recelem(i) = element
    end subroutine changeRecElement

    !"--------------------------------------------------------------------------------------"
    ! prepare src recalculate

    subroutine prepareRecalcSrc(this,nsrc,srcxz,srctype,srcstf,extwavelet,srcf0,srcfactor,srcangle_force,srcM, delay)
        type(srcVar) :: this
        integer :: nsrc
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcxz, srcM
        real(kind=CUSTOM_REAL), dimension(:), pointer :: srcf0,srcfactor
        real(kind=CUSTOM_REAL), dimension(:), pointer :: srcangle_force
        real(kind=custom_real), dimension(:), pointer :: delay
        integer, dimension(:), pointer :: srctype, srcstf
        character(len=80), dimension(:), pointer :: extwavelet
        real(kind=CUSTOM_REAL), dimension(2,nsrc) :: tempxz
        real(kind=CUSTOM_REAL), dimension(3,nsrc) :: tempm
        real(kind=CUSTOM_REAL), dimension(nsrc) :: tempf0,tempfactor
        real(kind=CUSTOM_REAL), dimension(nsrc) :: tempangle_force
        integer, dimension(nsrc) :: temptype, tempstf
        character(len=80), dimension(nsrc) :: tempextwavelet
        tempxz = srcxz
        temptype = srctype
        tempstf = srcstf
        tempextwavelet = extwavelet
        tempf0 = srcf0
        tempfactor = srcfactor
        tempangle_force = srcangle_force
        tempm = srcm

        ! create new pointer to avoid haning pointers
        allocate(this%srcxz(2,nsrc))
        allocate(this%srcrs(2,nsrc))
        allocate(this%srcelem(nsrc))
        allocate(this%srci(nsrc))
        allocate(this%srctype(nsrc))
        allocate(this%srcstf(nsrc))
        allocate(this%extwavelet(nsrc))
        allocate(this%srcf0(nsrc))
        allocate(this%srcfactor(nsrc))
        allocate(this%srcangle_force(nsrc))
        allocate(this%srcM(3,nsrc))
        this%nsrc = nsrc
        this%srcxz = tempxz
        this%srctype = temptype
        this%srcstf = tempstf
        this%extwavelet = tempextwavelet
        this%srcf0 = tempf0
        this%srcfactor = tempfactor
        this%srcangle_force = tempangle_force
        this%srcM = tempm
        this%delay = delay
    end subroutine prepareRecalcSrc

    !"--------------------------------------------------------------------------------------"
    ! prepare rec recalculate
    subroutine prepareRecalcRec(this,nrec,recxz,recnr)
        type(recVar) :: this
        integer :: nrec
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: recxz
        integer, dimension(:), pointer :: recnr
        real(kind=CUSTOM_REAL), dimension(2,nrec), target :: tempxz
        integer, dimension(nrec), target :: tempnr

        tempxz = recxz
        tempnr = recnr
        !produce new pointer to awoid haning pointers
        allocate(this%recxz(2,nrec))
        allocate(this%recrs(2,nrec))
        allocate(this%recelem(nrec))
        allocate(this%reci(nrec))
        allocate(this%recnr(nrec))
        this%nrec = nrec
        this%recxz = tempxz
        this%recnr = tempnr
    end subroutine prepareRecalcRec

    !"--------------------------------------------------------------------------------------"
    ! write srcVar
    subroutine writeSrcVar(this,filename)
        type(srcVar) :: this
        character(len=*) :: filename
        !write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a40, a16, a24)') "|            write soure database-file: ", trim(filename), "                       |"
        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        write(27) this%nsrc
        write(27) this%srcxz
        write(27) this%srcrs
        write(27) this%srcelem
        write(27) this%srci
        write(27) this%srctype
        write(27) this%srcstf
        write(27) this%extwavelet
        write(27) this%srcf0
        write(27) this%srcfactor
        write(27) this%srcangle_force
        write(27) this%srcM
        write(27) this%delay
        close(27)
    end subroutine writeSrcVar

    !"--------------------------------------------------------------------------------------"
    ! write recVar
    subroutine writeRecVar(this,filename)
        type(recVar) :: this
        character(len=*) :: filename
        write(*,'(a40, a16, a24)') "|         write receiver database-file: ", trim(filename), "                       |"
        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        write(27) this%nrec
        write(27) this%recxz
        write(27) this%recrs
        write(27) this%recelem
        write(27) this%reci
        write(27) this%recnr
        close(27)
    end subroutine writeRecVar

    !"--------------------------------------------------------------------------------------"
    ! read srcVar
    subroutine readSrcVar(this,filename)
        type(srcVar) :: this
        character(len=*) :: filename

        write(*,'(a40, a16, a24)') "|            read source database-file: ", trim(filename), "                       |"

        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        read(27) this%nsrc

        allocate(this%srcxz(2,this%nsrc))
        allocate(this%srcrs(2,this%nsrc))
        allocate(this%srcelem(this%nsrc))
        allocate(this%srci(this%nsrc))
        allocate(this%srctype(this%nsrc))
        allocate(this%srcstf(this%nsrc))
        allocate(this%extwavelet(this%nsrc))
        allocate(this%srcf0(this%nsrc))
        allocate(this%srcfactor(this%nsrc))
        allocate(this%srcangle_force(this%nsrc))
        allocate(this%srcM(3,this%nsrc))
        allocate(this%delay(this%nsrc))

        read(27) this%srcxz
        read(27) this%srcrs
        read(27) this%srcelem
        read(27) this%srci
        read(27) this%srctype
        read(27) this%srcstf
        read(27) this%extwavelet
        read(27) this%srcf0
        read(27) this%srcfactor
        read(27) this%srcangle_force
        read(27) this%srcM
        read(27) this%delay
        close(27)
    end subroutine readSrcVar

    !"--------------------------------------------------------------------------------------"
    ! read recVar
    subroutine readRecVar(this,filename)
        type(recVar) :: this
        character(len=*) :: filename
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a40, a16, a24)') "|          read receiver database-file: ", trim(filename), "                       |"

        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        read(27) this%nrec

        allocate(this%recxz(2,this%nrec))
        allocate(this%recrs(2,this%nrec))
        allocate(this%recelem(this%nrec))
        allocate(this%reci(this%nrec))
        allocate(this%recnr(this%nrec))
        read(27) this%recxz
        read(27) this%recrs
        read(27) this%recelem
        read(27) this%reci
        read(27) this%recnr
        close(27)
    end subroutine readRecVar

    !"--------------------------------------------------------------------------------------"
    ! deallocate srcVar
    subroutine deallocSrcVar(this)
        type(srcVar) :: this
        if (associated(this%srcxz)) deallocate(this%srcxz)
        if (associated(this%srcrs)) deallocate(this%srcrs)
        if (associated(this%srcelem)) deallocate(this%srcelem)
        if (associated(this%srci)) deallocate(this%srci)
        if (associated(this%srctype)) deallocate(this%srctype)
        if (associated(this%srcstf)) deallocate(this%srcstf)
        if (associated(this%extwavelet)) deallocate(this%extwavelet)
        if (associated(this%srcf0)) deallocate(this%srcf0)
        if (associated(this%srcfactor)) deallocate(this%srcfactor)
        if (associated(this%srcangle_force)) deallocate(this%srcangle_force)
        if (associated(this%srcm)) deallocate(this%srcm)
        if (allocated(this%delay)) deallocate(this%delay)
    end subroutine deallocSrcVar

    !"--------------------------------------------------------------------------------------"
    ! deallocate recVar
    subroutine deallocRecVar(this)
        type(recVar) :: this
        if (associated(this%recxz)) deallocate(this%recxz)
        if (associated(this%recrs)) deallocate(this%recrs)
        if (associated(this%recelem)) deallocate(this%recelem)
        if (associated(this%reci)) deallocate(this%reci)
        if (associated(this%recnr)) deallocate(this%recnr)
    end subroutine deallocRecVar
end module sourceReceiverMod
