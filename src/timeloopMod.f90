!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Marc S. Boxberg (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Thomas Möller (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2015-2017 Andre Lamert (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2016 Elena Busch (Rice University, United States)
!
!   This file is part of NEXD 2D.
!
!   NEXD 2D is free software: you can redistribute it and/or modify it 
!   under the terms of the GNU General Public License as published by the 
!   Free Software Foundation, either version 3 of the License, or (at your 
!   option) any later version.
!
!   NEXD 2D is distributed in the hope that it will be useful, but WITHOUT
!   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!   FITNESS FOR A PARTICULAR PURPOSE. 
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License v3.0
!   along with NEXD 2D. If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------
module timeloopMod
    ! module to calculate the timeloop
    use constantsMod
    use meshMod
    use parameterMod
    use waveMod
    use stfMod
    use sourceReceiverMod
    use matrixMod
    use plotMod
    use mpiMod
    use fileunitMod
    use pmlMod
    use logoMod
    use errorMessage
    use, intrinsic :: iso_fortran_env

    implicit none

    contains

    subroutine timeloop2d(par,mesh,src,rec,movie, myrank, errmsg)
        type(error_message) :: errmsg
        type(parameterVar) :: par
        type(meshVar) :: mesh
        type(srcVar) ::src
        type(recVar) ::rec
        type(movie_parameter) :: movie
        ! time vatiables
        integer, pointer :: nt,timeint
        real(kind=CUSTOM_REAL) :: dt
        real(kind=CUSTOM_REAL), pointer :: cfl
        real(kind=CUSTOM_REAL) :: f0,f0tmp
        real(kind=CUSTOM_REAL) :: plott0
        real(kind=CUSTOM_REAL), pointer, dimension(:) :: t0
        real(kind=CUSTOM_REAL) :: time
        real(kind=CUSTOM_REAL) :: fcrit,avg_energy1,avg_energy2,sta_lta
        integer :: timecrit
        logical :: pmlcrit = .true.
        ! mesh
        real(kind=CUSTOM_REAL), dimension(:), pointer :: vx,v
        ! absorbing
        ! free
        real, dimension(5,5) :: free1,free2
        ! indices
        integer :: i,j,r
        integer :: iglob
        integer, dimension(Np) :: iv
        integer :: ie, it
        ! RK
        integer :: irk, nrk, wrk
        ! rec variabels
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: recInt,recTemp
        real, dimension(:,:), allocatable :: plotv,plotw,plotux,plotuz,plotax,plotaz
        real(kind=CUSTOM_REAL), dimension(1) :: r_v,s_v
        real(kind=CUSTOM_REAL) :: ux_temp, uz_temp, vx_temp, vz_temp, ax_temp, az_temp, angle_temp
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: energy,all_energy, energy_kin, energy_pot
        integer, dimension(:), pointer :: recelemv
        ! source
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: srcInt,srcTemp
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotstf
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotDiffstf
        real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: srcArray
        real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: srcArrayM
        integer, dimension(:), pointer :: srcelemv
        real(kind=CUSTOM_REAL), dimension(2,2) :: M,Ms,Ma
        integer :: whichsource
        ! PML
        logical, dimension(:), allocatable :: pmlcheck
        integer :: pmlcnt
        real(kind=CUSTOM_REAL) :: dpmlx, dpmlz, rc, kmax, amax, afac
        real(kind=CUSTOM_REAL) :: vptemp
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ddx,ddz,alphax,alphaz,kx,kz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: fprime, gprime
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: fprimen, fprimem, gprimen, gprimem
        ! fields
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ux,uz,ax,az, uplot, vplot
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rQ,Q,Qn,Qm,e
        real(kind=CUSTOM_REAL), dimension(Np,5) :: ftemp,gtemp,qtemp
        real(kind=CUSTOM_REAL), dimension(Np) :: dFdr,dFds,dGdr,dGds
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Dr,Ds
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: invmass, vdmTinv
        real(kind=CUSTOM_REAL), dimension(Ngll*3,5) :: flux
        real(kind=CUSTOM_REAL) , dimension(:), allocatable :: stf

        logical :: file_exists
        character(len=80) ::filename
        ! usefulls
        real(kind=CUSTOM_REAL) :: onehalf = 1./2.
        real(kind=CUSTOM_REAL) :: onethree = 1./3.
        real(kind=CUSTOM_REAL) :: twothree = 2./3.
        real(kind=CUSTOM_REAL) :: onefor = 1./4.
        real(kind=CUSTOM_REAL) :: threefor = 3./4.
        ! timer
        integer :: t1,t2,max,rate
        real :: dt1,dt_ges
        ! MPI
        integer :: myrank
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_send
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_rec
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi_test
        integer :: c,k
        integer :: dest
        integer :: req, reqrec, tag
        integer ,dimension(:), allocatable:: req1
        real(kind=CUSTOM_REAL) :: maxv
        real(kind=CUSTOM_REAL) :: maxu
        integer :: elem_sum
        ! attenuation
        real(kind=CUSTOM_REAL), dimension(Np,3*nMB) :: a_ftemp,a_gtemp
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: theta,thetam,thetan
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: a_rQ
        real(kind=CUSTOM_REAL), dimension(Ngll*3,3) :: aflux
        ! adjoint
        logical :: inv_model
        integer :: iit, isave

        real(kind=CUSTOM_REAL), dimension(:), allocatable :: adj_sxx,adj_szz,adj_sxz,adj_ax,adj_az
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: adj_kernel_lambda, adj_kernel_mu, adj_kernel_rho
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: adj_kernel_vp, adj_kernel_vs, adj_kernel_rhop
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: adj_src_x,adj_src_z
        real(kind=CUSTOM_REAL) :: tmp

        ! error message
        character(len=10) :: myname = "timeloop2d"

        call addTrace(errmsg, myname)

        call sync_mpi()
        write(*,*) "Start MPI job" , myrank
        call sync_mpi()
        if (myrank == 0) write(*,*)
        call sync_mpi()

        call deallocMeshVar(mesh)

        write(filename,"('out/meshVar',i6.6)") myrank+1
        inquire(file=trim(filename), exist=file_exists)
        if (file_exists) then
            call readMeshVar(mesh,filename, errmsg)
        else
            call add(errmsg, 2, "Error in databases, file "//trim(filename)//" does not exist! Abort...", myname)
            call print(errmsg)
            call stop_mpi()
        end if

        call sync_mpi()
        if (myrank == 0) write(*,*)
        call sync_mpi()

        if (par%adjoint) then
            par%save_forward = .false.
        end if

        if (mesh%has_src) then
            write(filename,"('out/srcVar',i6.6)") myrank+1
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readSrcVar(src,filename)
            else
                call add(errmsg, 2, "Error in srcVar, files do not exist! Abort...", myname)
                call print(errmsg)
                call stop_mpi()
            end if
            write(*,*) "source in rank ",myrank, "and element",src%srcelem
            f0tmp = maxval(src%srcf0)
        else
            f0tmp = 0.0
        endif
        
        ! get the maximum f0
        call maxval_real_all(f0tmp,f0)
        plott0=1.2/f0
        
        ! load receivers
        if (mesh%has_rec) then
            write(filename,"('out/recVar',i6.6)") myrank+1
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readRecVar(rec,filename)
            else
                call add(errmsg, 2, "Error in recVar, files do not exist! Abort...", myname)
                call print(errmsg)
                call stop_mpi()
            end if
            write(*,*) 'm', myrank
            write(*,*) "receiver in rank ",myrank, "and element",rec%recelem
        end if

        call sync_mpi()
        if (myrank == 0) write(*,*)
        call sync_mpi()

        ! allocate fields
        allocate(ux(mesh%nglob),uz(mesh%nglob),ax(mesh%nglob),az(mesh%nglob),uplot(mesh%nglob),vplot(mesh%nglob))
        allocate(rQ(mesh%nglob,5),Q(mesh%nglob,5),Qn(mesh%nglob,5),Qm(mesh%nglob,5))
        allocate(e(mesh%nglob,3))
        allocate(a_rQ(mesh%nglob,3*nMB),theta(mesh%nglob,3,nMB),thetan(mesh%nglob,3,nMB),thetam(mesh%nglob,3,nMB))
        allocate(stf(par%nt))
        allocate(energy(par%nt),energy_kin(par%nt), energy_pot(par%nt))
        allocate(all_energy(par%nt))
        energy     = 0.0
        energy_kin = 0.0
        energy_pot = 0.0
        all_energy = 0.0
        if (.not.par%adjoint) then
            if (mesh%has_src) then
                allocate(srcInt(Np,src%nsrc),srcTemp(Np,src%nsrc))
                allocate(t0(src%nsrc))
                allocate(srcArray(Np,2,src%nsrc))
                allocate(srcArrayM(Np,3,src%nsrc))
                allocate(plotstf(par%nt,2,src%nsrc))
                allocate(plotDiffstf(par%nt,2,src%nsrc))
            end if
            if (mesh%has_rec) then
                allocate(recInt(Np,rec%nrec),recTemp(Np,rec%nrec))
                allocate(plotv(rec%nrec,par%nt),plotw(rec%nrec,par%nt),plotux(rec%nrec,par%nt),plotuz(rec%nrec,par%nt),plotax(rec%nrec,par%nt),plotaz(rec%nrec,par%nt))
            end if
        else !adjoint
            if (mesh%has_rec) then
                allocate(recInt(Np,rec%nrec),recTemp(Np,rec%nrec))
                allocate(adj_src_x(par%nt,2,rec%nrec),adj_src_z(par%nt,2,rec%nrec))
                allocate(srcArray(Np,2,rec%nrec))
            end if
        endif

        ! mpi
        allocate(q_send(NGLL*mesh%mpi_ne*5,mesh%mpi_nn))
        allocate(q_rec(NGLL*mesh%mpi_ne*5,mesh%mpi_nn))
        allocate(qi(NGLL,5,mesh%mpi_ne,mesh%mpi_nn))
        allocate(qi_test(NGLL,5,mesh%mpi_ne,mesh%mpi_nn))
        allocate(req1(mesh%mpi_nnmax))
        ! pml
        allocate(pmlcheck(mesh%nelem))
        pmlcheck=.false.

        if(par%set_pml) then
            allocate(fprime(mesh%nglob,5))
            allocate(gprime(mesh%nglob,5))
            allocate(fprimen(mesh%nglob,5),fprimem(mesh%nglob,5))
            allocate(gprimen(mesh%nglob,5),gprimem(mesh%nglob,5))
        end if
        if (.not.par%adjoint) then
            if (mesh%has_rec) then
                plotux = 0.0
                plotuz = 0.0
                plotv  = 0.0
                plotw  = 0.0
                plotax = 0.0
                plotaz = 0.0
            end if
        endif
        q      = 1e-24!eps
        e      = 1e-24!eps
        rQ     = 1e-24!eps
        Qn     = 1e-24!eps
        if (par%attenuation) then
           a_rQ   = 1e-24!eps
           theta  = 1e-24!eps
           thetan = 1e-24!eps
        endif

        ux    = 0.0
        uz    = 0.0
        ax    = 0.0
        az    = 0.0
        uplot = 0.0
        vplot = 0.0

        if (par%autodt) then
            dt = mesh%dtfactor*par%cfl
        else
            dt = par%dt
        endif

        t1     = 0.0
        t2     = 0.0
        dt_ges = 0.0

        fcrit    = 5/f0
        timecrit = ((1.2/f0)*3)/dt
        ! sta lta properties for avoiding pml instabilities
        if (par%log.and.myrank==0) write(*,*) "timecrit: ",timecrit
        if (par%log.and.myrank==0) write(*,*) "avg_window1 :",par%avg_window1
        if (par%log.and.myrank==0) write(*,*) "avg_window2 :",par%avg_window2
        if (par%log.and.myrank==0) write(*,*) "sta_lta_trigger:",par%sta_lta_trigger

        if (par%timeint == 1) then
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) " use euler"
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            nrk = 1 !euler
            wrk = 1
        else if (par%timeint == 2) then
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) " use tvd runge kutta second order"
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            nrk = 2 ! rk2
            wrk = 1
        else if (par%timeint == 3) then
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) " use tvd runge kutta third order"
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            nrk = 3 ! rk3
            wrk = 2
        end if

        if (par%inv_model) then
            mesh%rho(:) = 0.0
            mesh%vp(:) = 0.0
            mesh%vs(:) = 0.0
            mesh%lambda(:) = 0.0
            mesh%mu(:) = 0.0
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) " read external velocity model"
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            write(filename,"('save/model/model_itter',i6.6,'.dat')") myrank+1
            open(unit = 76, file = trim(filename),status = 'unknown')
            do ie=1,mesh%nelem ! element loop
                do j=1,Np
                    i=mesh%ibool(j,ie)
                    read(76,*) tmp, tmp, tmp, mesh%rho(ie), mesh%vp(ie), mesh%vs(ie)
                end do
                mesh%mu(ie) = mesh%rho(ie) * mesh%vs(ie)**2
                mesh%lambda(ie) = mesh%rho(ie) * mesh%vp(ie)**2 - 2 * mesh%rho(ie)*mesh%vs(ie)**2
            end do
            close(76)
        endif

        if(par%set_pml) then
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) " using NPML boundary conditions ", par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,par%pml_kmax,par%pml_afac
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            pmlcnt=0
            amax=par%pml_afac*pi*f0
            dpmlx=par%pml_delta
            dpmlz=par%pml_delta

            allocate(ddx(mesh%nglob), ddz(mesh%nglob))
            allocate(alphax(mesh%nglob), alphaz(mesh%nglob))
            allocate(kx(mesh%nglob),kz(mesh%nglob))
            ddx=0
            ddz=0
            kx=1
            kz=1

            ! setup npml damping profiles
            ! 1: left
            ! 2: bottom
            ! 3: right
            ! 4: top
            ! 5: bottom left
            ! 6: bottom right
            ! 7: top right
            ! 8: top left

            do i=1,mesh%nelem
                vptemp=mesh%vp(i)
                ! assign elements to pml
                if (mesh%pml(i)< 1) then
                    pmlcheck(i)=.false.
                else
                    pmlcheck(i)=.true.
                    pmlcnt=pmlcnt+1
                end if
                ! damping functions for pml
                iv=mesh%ibool(:,i)
                if ( mesh%pml(i) == 1) then !left
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(i),par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,1)
                else if ( mesh%pml(i) == 2 ) then !bottom
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(i),par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,2)
                else if ( mesh%pml(i) == 3 ) then !right
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(i),par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,3)
                else if ( mesh%pml(i) == 4 ) then !top
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(i),par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,4)
                else if ( mesh%pml(i) == 5) then !bottom left
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(i),par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,5)
                else if ( mesh%pml(i) == 6 ) then !bottom right
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(i),par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,6)
                else if ( mesh%pml(i) == 7 ) then !top right
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(i),par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,7)
                else if ( mesh%pml(i) == 8 ) then !top left
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(i),par%xmin,par%xmax,par%zmin,par%zmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,8)
                end if
            end do
        end if

        ! free surface conditions
        if (par%strongform) then
            free1=0.0
            free1(1,1)=-2
            free1(2,2)=0
            free1(3,3)=-2
            free1(4,4)=0
            free1(5,5)=0

            free2=0.0
            free2(1,1)=0
            free2(2,2)=0
            free2(3,3)=0
            free2(4,4)=-2
            free2(5,5)=-2
        else
            free1=0.0
            free1(1,1)=-1
            free1(2,2)=1
            free1(3,3)=-1
            free1(4,4)=1
            free1(5,5)=1

            free2=0.0
            free2(1,1)=0
            free2(2,2)=0
            free2(3,3)=0
            free2(4,4)=-1
            free2(5,5)=-1
        end if

        vdmTinv=mesh%vdm
        vdmTinv=transpose(vdmTinv)
        call invert(vdmTinv, errmsg)
        allocate(srcelemV(mesh%nelem))
        allocate(recelemV(mesh%nelem))
        srcelemv=0
        recelemv=0
        if (.not.par%adjoint) then
            if (mesh%has_src) then
                do i=1,src%nsrc
                    if (src%srctype(i) == 0) then ! singleforce
                        whichsource = 0
                        r_v(1) = src%srcrs(1,i)
                        s_v(1) = src%srcrs(2,i)
                        call vdm2D(srcTemp(:,i),r_v,s_v)
                        srcInt(:,i)=matmul(vdmTinv,srcTemp(:,i))

                        ! set t0 to a save value
                        t0(i)=1.2/src%srcf0(i)
                        srcelemv(src%srcelem(i)) = i
                    else if (src%srctype(i) == 1) then !momenttensor
                        whichsource = 1
                        r_v(1) = src%srcrs(1,i)
                        s_v(1) = src%srcrs(2,i)
                        call vdm2D(srcTemp(:,i),r_v,s_v)
                        srcInt(:,i)=matmul(vdmTinv,srcTemp(:,i))
                        ! set t0 to a save value
                        t0(i)=1.2/src%srcf0(i)
                        srcelemv(src%srcelem(i)) = i

                        M(1,1)=src%srcm(1,i) !Mxx
                        M(1,2)=src%srcm(3,i) !Mxz
                        M(2,1)=src%srcm(3,i) !Mzx
                        M(2,2)=src%srcm(2,i) !Mzz

                        do j=1,2
                            do k=1,2
                                ms(j,k)=0.5*(M(j,k)+M(k,j))
                                ma(j,k)=0.5*(M(j,k)-M(k,j))
                            end do
                        end do
                    end if
                end do
            end if
        end if
        
        ! set up rec interpolation
        if (mesh%has_rec) then
            do i=1,rec%nrec
                r_v(1) = rec%recrs(1,i)
                s_v(1) = rec%recrs(2,i)
                call vdm2D(recTemp(:,i),r_v,s_v)
                recInt(:,i)=matmul(vdmTinv,recTemp(:,i))
                recelemv(rec%recelem(i)) = i
            end do
        end if

        ! choose source time function
        if (.not.par%adjoint) then
            if (mesh%has_src) then
                do i=1,src%nsrc
                    do it=1,par%nt
                        time = (float(it)-1.)*dt
                        if (src%srcstf(i) == 1) then
                            plotstf(it,1,i) = time
                            plotDiffstf(it,1,i) = time
                            plotstf(it,2,i) = -stfGauss(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                            plotDiffstf(it,2,i) = -stfDiffGauss(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                        else if (src%srcstf(i) == 2) then !RICKER
                            plotstf(it,1,i) = time
                            plotDiffstf(it,1,i) = time
                            plotstf(it,2,i) = -stfRicker(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                            plotDiffstf(it,2,i) = -stfDiffRicker(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                        else if (src%srcstf(i) == 3) then !Sin^3
                            plotstf(it,1,i) = time
                            plotstf(it,2,i) = -stfSin3(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                        endif
                    end do
                    if (src%srcstf(i) == 4) then
                        filename=src%extwavelet(i)
                        call stfExternal(plotstf(:,2,i),plotstf(:,1,i),dt,par%nt,.true.,6,3,trim(filename),0, errmsg)
                        call stfExternal(plotDiffstf(:,2,i),plotDiffstf(:,1,i),dt,par%nt,.true.,6,3,trim(filename),1, errmsg)
                    end if

                    ! write stf
                    write(filename,"('out/stf',i6.6)") i
                    write(*,*) filename
                    open(unit=27,file=trim(filename),status='unknown')
                    do it=1,par%nt
                        write(27,*) plotstf(it,1,i),plotstf(it,2,i)
                    end do
                    close(27)
                    write(filename,"('out/stfdiff',i6.6)") i
                    write(*,*) filename
                    open(unit=27,file=trim(filename),status='unknown')
                    do it=1,par%nt
                        write(27,*) plotDiffstf(it,1,i),plotDiffstf(it,2,i)
                    end do
                    close(27)
                end do
            end if
        endif

        ! LL LL adjoint save
        if (par%save_forward) then
            isave = 1
            write(filename,"('/ax',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=86, file=trim(filename), form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='replace')

            write(filename,"('/az',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=87, file=trim(filename), form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='replace')

            write(filename,"('/sxx',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=88, file=trim(filename), form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='replace')

            write(filename,"('/szz',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=89, file=trim(filename), form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='replace')

            write(filename,"('/sxz',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=90, file=trim(filename), form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='replace')
        endif

        if (par%adjoint) then
            isave = par%nt/par%adj_nstep
            allocate(adj_sxx(mesh%nglob),adj_szz(mesh%nglob),adj_sxz(mesh%nglob) )
            allocate(adj_ax(mesh%nglob),adj_az(mesh%nglob))
            allocate(adj_kernel_lambda(mesh%nglob),adj_kernel_mu(mesh%nglob),adj_kernel_rho(mesh%nglob))
            allocate(adj_kernel_vp(mesh%nglob),adj_kernel_vs(mesh%nglob),adj_kernel_rhop(mesh%nglob))
            adj_kernel_lambda = 0.0
            adj_kernel_mu = 0.0
            adj_kernel_rho = 0.0
            adj_kernel_vp = 0.0
            adj_kernel_vs = 0.0

            write(filename,"('/ax',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=86, file=trim(filename), action='readwrite', form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='old')
            write(filename,"('/az',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=87, file=trim(filename), action='readwrite', form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='old')

            write(filename,"('/sxx',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=88, file=trim(filename), action='readwrite', form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='old')
            write(filename,"('/szz',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=89, file=trim(filename), action='readwrite', form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='old')
            write(filename,"('/sxz',i6.6,'.bin')") myrank+1
            filename="out"//trim(filename)
            open(unit=90, file=trim(filename), action='readwrite', form='unformatted', &
                access='direct',recl=CUSTOM_REAL*(mesh%nglob),status='old')

            ! adjoint sources
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) "read adjoint sources"
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (mesh%has_rec) then
                do r=1,rec%nrec
                    write(filename,"('out/seismo.x.',i7.7,'.sda.adj')") rec%recnr(r)
                    open(unit=45,file=trim(filename),action = 'read')
                    do i=1,par%nt
                        iit = (par%nt + 1) - i
                        read(45,*) adj_src_x(iit,1,r), adj_src_x(iit,2,r)
                    end do
                    close(45)
                    write(filename,"('out/seismo.z.',i7.7,'.sda.adj')") rec%recnr(r)
                    open(unit=45,file=trim(filename),action = 'read')
                    do i=1,par%nt
                        iit = (par%nt + 1) - i
                        read(45,*) adj_src_z(iit,1,r), adj_src_z(iit,2,r)
                    end do
                    close(45)
                end do
            endif
        endif

        invmass=matmul(mesh%vdm,transpose(mesh%vdm))
        call sum_int(mesh%nelem, elem_sum)
        if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log.and.myrank==0) write(*,*) "Start timeloop over" ,elem_sum," elements with dt: ", dt," and nt : ",par%nt, " and t0: ",plott0
        if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"

        ! ---------------------------------------------------------------------------------------------
        ! ------------------------------------timeloop-------------------------------------------------
        ! ---------------------------------------------------------------------------------------------

        if (par%log .and. myrank==0) then
            call system_clock(t1, rate, max) !get timer
        end if

        do it=1,par%nt ! timeloop
            stf(it)=(it-1)*dt
            if (.not.par%adjoint) then
                if (mesh%has_src) then
                    if (whichsource==0) then
                        do i=1,src%nsrc
                            srcArray(:,1,i)=srcInt(:,i)*(-sin((src%srcangle_force(i)*PI)/180.)*plotstf(it,2,i))/mesh%rho(src%srcelem(i))
                            srcArray(:,2,i)=srcInt(:,i)*(+cos((src%srcangle_force(i)*PI)/180.)*plotstf(it,2,i))/mesh%rho(src%srcelem(i))

                            iv=mesh%ibool(:,src%srcelem(i))
                            srcArray(:,1,i)=matmul(invmass,srcArray(:,1,i))/mesh%jacobian(iv)
                            srcArray(:,2,i)=matmul(invmass,srcArray(:,2,i))/mesh%jacobian(iv)
                        end do
                    else if (whichsource==1) then
                        do i=1,src%nsrc
                            srcArrayM(:,1,i)=srcInt(:,i)*(ms(1,1)*plotDiffstf(it,2,i))
                            srcArrayM(:,2,i)=srcInt(:,i)*(ms(2,2)*plotDiffstf(it,2,i))
                            srcArrayM(:,3,i)=srcInt(:,i)*(ms(1,2)*plotDiffstf(it,2,i))

                            iv=mesh%ibool(:,src%srcelem(i))
                            srcArrayM(:,1,i)=matmul(invmass,srcArrayM(:,1,i))/mesh%jacobian(iv)
                            srcArrayM(:,2,i)=matmul(invmass,srcArrayM(:,2,i))/mesh%jacobian(iv)
                            srcArrayM(:,3,i)=matmul(invmass,srcArrayM(:,3,i))/mesh%jacobian(iv)
                        end do
                    end if
                end if
                ! interpolate receivers
                if (mesh%has_rec) then
                    do r=1,rec%nrec
                        do j=1,Np
                            iglob=mesh%ibool(j,rec%recelem(r))
                            plotux(r,it) = plotux(r,it) + ux(iglob)*recint(j,r)
                            plotuz(r,it) = plotuz(r,it) + uz(iglob)*recint(j,r)
                            plotv(r,it)  = plotv(r,it)  + q(iglob,4)*recint(j,r)
                            plotw(r,it)  = plotw(r,it)  + q(iglob,5)*recint(j,r)
                            plotax(r,it) = plotax(r,it) + ax(iglob)*recint(j,r)
                            plotaz(r,it) = plotaz(r,it) + az(iglob)*recint(j,r)
                        end do
                        ! rotate receivers
                        angle_temp = (par%rec_angle*PI)/180.
                        ux_temp = cos(angle_temp) * plotux(r,it) + sin(angle_temp) * plotuz(r,it)
                        uz_temp = -sin(angle_temp) * plotux(r,it) + cos(angle_temp) * plotuz(r,it)
                        vx_temp = cos(angle_temp) * plotv(r,it) + sin(angle_temp) * plotw(r,it)
                        vz_temp = -sin(angle_temp) * plotv(r,it) + cos(angle_temp) * plotw(r,it)
                        ax_temp = cos(angle_temp) * plotax(r,it) + sin(angle_temp) * plotaz(r,it)
                        az_temp = -sin(angle_temp) * plotax(r,it) + cos(angle_temp) * plotaz(r,it)
                        plotux(r,it) = ux_temp
                        plotuz(r,it) = uz_temp
                        plotv(r,it) = vx_temp
                        plotw(r,it) = vz_temp
                        plotax(r,it) = ax_temp
                        plotaz(r,it) = az_temp
                    end do
                end if
            else
                if (mesh%has_rec) then
                    do i=1,rec%nrec
                        srcArray(:,1,i) = recInt(:,i)*adj_src_x(it,2,i)/mesh%rho(rec%recelem(i))
                        srcArray(:,2,i) = recInt(:,i)*adj_src_z(it,2,i)/mesh%rho(rec%recelem(i))
                        iv = mesh%ibool(:,rec%recelem(i))
                        srcArray(:,1,i) = matmul(invmass,srcArray(:,1,i))/mesh%jacobian(iv)
                        srcArray(:,2,i) = matmul(invmass,srcArray(:,2,i))/mesh%jacobian(iv)
                    end do
                end if
            endif

            tag=0
            do irk=1,nrk ! runge-kutta loop
                qm=q
                if (par%attenuation) then
                    thetam=theta
                end if
                if(par%set_pml) then
                    fprimem=fprime
                    gprimem=gprime
                end if
                q_send(:,:) = 0

                ! mpi comunication
                ! build send buffer
                do i=1,mesh%mpi_nn
                    do ie=1,mesh%mpi_ne ! loop over interface elements
                        do k=1,5
                            do j=1,NGLL
                                if ( mesh%mpi_connection(i,ie,1) >0) then
                                    q_send((ie-1)*5*NGLL + (k-1)*NGLL + j,i) = &
                                        qm(mesh%ibt(j,mesh%mpi_connection(i,ie,2),mesh%mpi_connection(i,ie,1)),k)
                                end if
                            end do
                        end do
                    end do
                end do ! all interfaces

                ! send and rec
                do i=1,mesh%mpi_nn
                    dest=mesh%mpi_neighbor(i)-1
                    call isendV_real(q_send(:,i),(mesh%mpi_ne*5*NGLL),dest,tag,req)
                    call irecV_real(q_rec(:,i),(mesh%mpi_ne*5*NGLL),dest,tag,reqrec)
                    call wait_req(req)
                    call wait_req(reqrec)
                end do

                ! unpack mpi buffers
                do i=1,mesh%mpi_nn
                    c = 1
                    do ie=1,mesh%mpi_ne
                        do k=1,5
                            do j=1,NGLL
                                if ( mesh%mpi_connection(i,ie,1) > 0) then
                                    qi(j,k,ie,i) = q_rec(c,i)
                                end if
                                c = c+1
                            end do
                        end do
                    end do
                end do

                do ie=1,mesh%nelem ! element loop
                    iv=mesh%ibool(:,ie)
                    qtemp=qm(iv,:)

                    ! get fluxes
                    call elasticFluxes(qtemp,mesh%lambda(ie),mesh%mu(ie),mesh%rho(ie),ftemp,gtemp)

                    if (pmlcheck(ie)) then
                        do i=1,5
                            ftemp(:,i) = (fprimem(iv,i) + ftemp(:,i))/kx(iv)
                            gtemp(:,i) = (gprimem(iv,i) + gtemp(:,i))/kz(iv)
                        enddo
                    end if

                    if (par%strongform) then
                        ! comp strong deriverate
                        do i=1,5
                            dFdr = matmul(mesh%Dr,ftemp(:,i))
                            dFds = matmul(mesh%Ds,ftemp(:,i))
                            dGdr = matmul(mesh%Dr,gtemp(:,i))
                            dGds = matmul(mesh%Ds,gtemp(:,i))
                            rQ(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds) + (mesh%rz(iv)*dGdr + mesh%sz(iv)*dGds)
                        end do
                    else
                        ! comp weak deriverate
                        do i=1,5
                            dFdr = matmul(mesh%Drw,ftemp(:,i))
                            dFds = matmul(mesh%Dsw,ftemp(:,i))
                            dGdr = matmul(mesh%Drw,gtemp(:,i))
                            dGds = matmul(mesh%Dsw,gtemp(:,i))
                            rQ(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds) + (mesh%rz(iv)*dGdr + mesh%sz(iv)*dGds)
                        end do
                    endif

                    if (par%attenuation) then
                        call anelasticFluxes(qtemp,mesh%wl(:,ie),a_ftemp,a_gtemp)
                        if (par%strongform) then
                            do i=1,3*nMB
                                dFdr = matmul(mesh%Dr,a_ftemp(:,i))
                                dFds = matmul(mesh%Ds,a_ftemp(:,i))
                                dGdr = matmul(mesh%Dr,a_gtemp(:,i))
                                dGds = matmul(mesh%Ds,a_gtemp(:,i))
                                a_rQ(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds) + (mesh%rz(iv)*dGdr + mesh%sz(iv)*dGds)
                            end do
                        else
                            do i=1,3*nMB
                                dFdr = matmul(mesh%Drw,a_ftemp(:,i))
                                dFds = matmul(mesh%Dsw,a_ftemp(:,i))
                                dGdr = matmul(mesh%Drw,a_gtemp(:,i))
                                dGds = matmul(mesh%Dsw,a_gtemp(:,i))
                                a_rQ(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds) + (mesh%rz(iv)*dGdr + mesh%sz(iv)*dGds)
                            end do
                        end if
                    end if

                    ! compute fluxes on the surface
                    if (par%strongform) then
                        call computeExactRiemannSF(flux,qm,qi,mesh%neighbor(:,ie),mesh%vp(ie),mesh%vs(ie),mesh%rho(ie),&
                            mesh%lambda(ie),mesh%mu(ie),mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie) &
                            ,mesh%ibt(:,:,ie),mesh%ibn(:,:,ie),mesh%nx(:,ie),mesh%nz(:,ie),free1,.false.,par%set_pml)
                    else
                        call computeExactRiemannWF(flux,qm,qi,mesh%neighbor(:,ie),mesh%vp(ie),mesh%vs(ie),mesh%rho(ie),&
                            mesh%lambda(ie),mesh%mu(ie),mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie) &
                            ,mesh%ibt(:,:,ie),mesh%ibn(:,:,ie),mesh%nx(:,ie),mesh%nz(:,ie),free1,.false.,par%set_pml)
                    endif

                    if(par%attenuation) then
                        if (par%strongform) then
                            call computeExactRiemannSFAnelastic(aflux,qm,qi,mesh%neighbor(:,ie),mesh%vpu(ie),mesh%vsu(ie),mesh%rho(ie),&
                                mesh%lambdau(ie),mesh%muu(ie),mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie) &
                                ,mesh%ibt(:,:,ie),mesh%ibn(:,:,ie),mesh%nx(:,ie),mesh%nz(:,ie),free1,.false.,par%set_pml)
                        else
                            call computeExactRiemannWFAnelastic(aflux,qm,qi,mesh%neighbor(:,ie),mesh%vpu(ie),mesh%vsu(ie),mesh%rho(ie),&
                                mesh%lambdau(ie),mesh%muu(ie),mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie) &
                                ,mesh%ibt(:,:,ie),mesh%ibn(:,:,ie),mesh%nx(:,ie),mesh%nz(:,ie),free1,.false.,par%set_pml)
                        end if
                        c=1
                        do i=1,nMB
                            do j=1,3
                                a_rQ(iv,c) = a_rQ(iv,c) + matmul(mesh%lift,mesh%fscale(:,ie)*mesh%wl(i,ie)*aflux(:,j)/2.0)
                                a_rQ(iv,c) = a_rQ(iv,c) - mesh%wl(i,ie) * thetam(iv,j,i)
                                c=c+1
                            end do
                        end do

                        do j=1,nMB
                            rq(iv,1) = rq(iv,1) +( - (mesh%lambdau(ie)*mesh%ylambda(j,ie) + 2*mesh%muu(ie)*mesh%ymu(j,ie))*thetam(iv,1,j) - mesh%lambdau(ie) * mesh%ylambda(j,ie) * thetam(iv,2,j))
                            rq(iv,2) = rq(iv,2) +( - mesh%lambdau(ie)*mesh%ylambda(j,ie)*thetam(iv,1,j) - (mesh%lambdau(ie)* mesh%ylambda(j,ie) + 2*mesh%muu(ie)*mesh%ymu(j,ie)) * thetam(iv,2,j))
                            rq(iv,3) = rq(iv,3) +( - 2*mesh%muu(ie)*mesh%ymu(j,ie) * thetam(iv,3,j))
                        end do
                    end if

                    ! lift surface integral
                    if (par%strongform) then
                        do i=1,5
                            rq(iv,i) = rq(iv,i) + matmul(mesh%lift,mesh%fscale(:,ie)*flux(:,i)/2.0)
                        end do
                    else
                        do i=1,5
                            rq(iv,i) = rq(iv,i) - matmul(mesh%lift,mesh%fscale(:,ie)*flux(:,i)/2.0)
                        end do
                    endif

                    ! do runge kutta
                    ! RK 2 and Euler
                    if (wrk==1) then
                        if (irk==1) then
                            qn(iv,:)=q(iv,:)
                            if(par%set_pml) then
                                fprimen(iv,:)=fprime(iv,:)
                                gprimen(iv,:)=gprime(iv,:)
                            end if
                            if(par%attenuation) then
                                thetan(iv,:,:) = theta(iv,:,:)
                            end if
                            ! if not source
                            if ( (.not.par%adjoint .and. mesh%has_src .and. (srcelemv(ie)>0)).and.whichsource==0) then
                                j=srcelemv(ie)
                                q(iv,4)=q(iv,4) + dt*(rq(iv,4)+srcArray(:,1,j))
                                q(iv,5)=q(iv,5) + dt*(rq(iv,5)+srcArray(:,2,j))
                                q(iv,1:3)=q(iv,1:3) + dt*rq(iv,1:3)
                            else if ( (.not.par%adjoint .and. mesh%has_src .and. (srcelemv(ie)>0)).and.whichsource==1) then
                                j=srcelemv(ie)
                                q(iv,1:3) = q(iv,1:3)+dt*(rq(iv,1:3) + srcArrayM(:,1:3,j))
                                q(iv,4:5)=q(iv,4:5) + dt*rq(iv,4:5)
                            else if (par%adjoint .and. (mesh%has_rec .and. (recelemv(ie) > 0))) then
                                j=recelemv(ie)
                                q(iv,4)=q(iv,4) + dt*(rq(iv,4)+srcArray(:,1,j))
                                q(iv,5)=q(iv,5) + dt*(rq(iv,5)+srcArray(:,2,j))
                                q(iv,1:3)=q(iv,1:3) + dt*rq(iv,1:3)
                            else
                                q(iv,1:5)=q(iv,1:5) + dt*rq(iv,1:5)
                            end if
                            !pml
                            if (pmlcheck(ie)) then
                                do i=1,5
                                    fprime(iv,i) = fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)*(ftemp(:,i)))
                                    gprime(iv,i) = gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)*(gtemp(:,i)))
                                end do
                            end if
                            ! attenuation
                            if(par%attenuation) then
                                c = 1
                                do i=1,nMB
                                    do j=1,3
                                        theta(iv,j,i) = theta(iv,j,i) + dt*(a_rQ(iv,c))
                                        c = c + 1
                                    end do
                                end do
                            end if
                        end if !rk==1
                        if (irk==2) then
                            if ( (.not.par%adjoint .and. mesh%has_src .and. (srcelemv(ie)>0)).and.whichsource==0) then
                                j=srcelemv(ie)
                                q(iv,1:3)=onehalf*qn(iv,1:3)+onehalf*(q(iv,1:3) + dt*rq(iv,1:3))
                                q(iv,4)=onehalf*qn(iv,4)+onehalf*(q(iv,4) + dt*(rq(iv,4)+srcArray(:,1,j)))
                                q(iv,5)=onehalf*qn(iv,5)+onehalf*(q(iv,5) + dt*(rq(iv,5)+srcArray(:,2,j)))
                            else if ( (.not.par%adjoint .and. mesh%has_src .and. (srcelemv(ie)>0)).and.whichsource==1) then
                                j=srcelemv(ie)
                                q(iv,1:3)=onehalf*qn(iv,1:3)+onehalf*(q(iv,1:3) + dt*(rq(iv,1:3) +srcArrayM(:,1:3,j) ))
                                q(iv,4:5)=onehalf*qn(iv,4:5)+onehalf*(q(iv,4:5) + dt*rq(iv,4:5))
                            else if (par%adjoint .and. (mesh%has_rec .and. (recelemv(ie) > 0))) then
                                j=recelemv(ie)
                                q(iv,1:3)=onehalf*qn(iv,1:3)+onehalf*(q(iv,1:3) + dt*rq(iv,1:3))
                                q(iv,4)=onehalf*qn(iv,4)+onehalf*(q(iv,4) + dt*(rq(iv,4)+srcArray(:,1,j)))
                                q(iv,5)=onehalf*qn(iv,5)+onehalf*(q(iv,5) + dt*(rq(iv,5)+srcArray(:,2,j)))
                            else
                                q(iv,1:5)=onehalf*qn(iv,1:5)+onehalf*(q(iv,1:5) + dt*rq(iv,1:5))
                            end if

                            if (pmlcheck(ie)) then
                                do i=1,5
                                    fprime(iv,i) = onehalf* fprimen(iv,i) + onehalf*(fprime(iv,i) + &
                                                   dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)*(ftemp(:,i))))
                                    gprime(iv,i) = onehalf* gprimen(iv,i) + onehalf*(gprime(iv,i) + &
                                                   dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)*(gtemp(:,i))))
                                end do
                            end if
                            ! attenuation
                            if(par%attenuation) then
                                c = 1
                                do i=1,nMB
                                    do j=1,3
                                        theta(iv,j,i) = onehalf*thetan(iv,j,i) + onehalf*(theta(iv,j,i) + dt*(a_rQ(iv,c)))
                                        c = c + 1
                                    end do
                                end do
                            end if
                        end if !rk==2
                    ! RK 3
                    else if ( wrk==2) then
                        if (irk==1) then
                            qn(iv,:)=q(iv,:)
                            ! if not source
                            if ( (mesh%has_src .and. (srcelemv(ie)>0)).and.whichsource==0) then
                                j=srcelemv(ie)
                                q(iv,4) = q(iv,4) + dt*(rq(iv,4)+srcArray(:,1,j))
                                q(iv,5) = q(iv,5) + dt*(rq(iv,5)+srcArray(:,2,j))
                            else
                                q(iv,4) = q(iv,4) + dt*rq(iv,4)
                                q(iv,5) = q(iv,5) + dt*rq(iv,5)
                            end if
                            q(iv,1) = q(iv,1) + dt*rq(iv,1)
                            q(iv,2) = q(iv,2) + dt*rq(iv,2)
                            q(iv,3) = q(iv,3) + dt*rq(iv,3)
                        end if !rk==1
                        if (irk==2) then
                            if ( (mesh%has_src .and. (srcelemv(ie)>0)).and.whichsource==0) then
                                j=srcelemv(ie)
                                q(iv,4) = threefor*qn(iv,4)+onefor*(q(iv,4) + dt*(rq(iv,4)+srcArray(:,1,j)))
                                q(iv,5) = threefor*qn(iv,5)+onefor*(q(iv,5) + dt*(rq(iv,5)+srcArray(:,2,j)))
                            else
                                q(iv,4) = threefor*qn(iv,4)+onefor*(q(iv,4) + dt*rq(iv,4))
                                q(iv,5) = threefor*qn(iv,5)+onefor*(q(iv,5) + dt*rq(iv,5))
                            end if
                                q(iv,1) = threefor*qn(iv,1)+onefor*(q(iv,1) + dt*rq(iv,1))
                                q(iv,2) = threefor*qn(iv,2)+onefor*(q(iv,2) + dt*rq(iv,2))
                                q(iv,3) = threefor*qn(iv,3)+onefor*(q(iv,3) + dt*rq(iv,3))
                        end if !rk==2
                        if (irk==3) then
                            if ( (mesh%has_src .and. (srcelemv(ie)>0)).and.whichsource==0) then
                                j=srcelemv(ie)
                              q(iv,4) = onethree*qn(iv,4) + twothree*(q(iv,4) + dt*(rq(iv,4)+srcArray(:,1,j)))
                              q(iv,5) = onethree*qn(iv,5) + twothree*(q(iv,5) + dt*(rq(iv,5)+srcArray(:,2,j)))
                           else
                              q(iv,4) = onethree*qn(iv,4) + twothree*(q(iv,4) + dt*rq(iv,4))
                              q(iv,5) = onethree*qn(iv,5) + twothree*(q(iv,5) + dt*rq(iv,5))
                           end if
                           q(iv,1) = onethree*qn(iv,1) + twothree*(q(iv,1) + dt*rq(iv,1))
                           q(iv,2) = onethree*qn(iv,2) + twothree*(q(iv,2) + dt*rq(iv,2))
                           q(iv,3) = onethree*qn(iv,3) + twothree*(q(iv,3) + dt*rq(iv,3))
                        end if !rk==3
                    end if

                    !Displacement and uplot so norm u will be calculuated
                    if (irk==nrk) then
                        ux(iv) = ux(iv) + dt*q(iv,4)
                        uz(iv) = uz(iv) + dt*q(iv,5)
                        uplot(iv) = sqrt(ux(iv)**2+uz(iv)**2)
                        if (movie%movie) then
                            if (movie%velocity) then
                                !Calculate the norm of the velocity
                                vplot(iv) = sqrt(q(iv,1)**2+q(iv,2)**2)
                            end if
                        end if

                        !Acceleration
                        ax(iv) = rq(iv,4)
                        az(iv) = rq(iv,5)
                        ! LL LL strain
                        e(iv,1) = (-mesh%lambda(ie) * q(iv,2) + &
                                  (2 * mesh%mu(ie) - mesh%lambda(ie)) * q(iv,1))/ (4*mesh%mu(ie)**2 + 4 * mesh%lambda(ie)*mesh%mu(ie))
                        e(iv,2) = q(iv,3)/(2 * mesh%mu(ie))
                        e(iv,3) = ((2 * mesh%mu(ie) + mesh%lambda(ie)) * q(iv,2) - &
                                  mesh%lambda(ie)*q(iv,1))/(4*mesh%mu(ie)**2 + 4*mesh%lambda(ie)*mesh%mu(ie))
                    endif
                end do !element loop
                call sync_mpi()
            end do ! RK-LOOP

            do ie=1,mesh%nelem
                iv=mesh%ibool(:,ie)
                do i=1,Np
                    energy_kin(it) = energy_kin(it) + mesh%rho(ie) * (q(iv(i),4)**2 + q(iv(i),5)**2)
                    energy_pot(it) = energy_pot(it) + 0.5 * (q(iv(i),1)*e(iv(i),1) + q(iv(i),2)*e(iv(i),2) + q(iv(i),3)*e(iv(i),3))
                end do
            end do
            energy(it) = energy_kin(it) + energy_pot(it)
            call sum_all_real(energy(it),all_energy(it))

            ! check if energy grows to much, its done with sta/lta trigger to switch off the pml if nessersary
            if(par%use_trigger) then
                if(par%set_pml) then
                    if(pmlcrit) then
                        if( (it>timecrit) ) then
                            avg_energy1=0
                            do i=1,par%avg_window1
                                avg_energy1=avg_energy1+all_energy(it-i)
                            end do
                            avg_energy1=avg_energy1/par%avg_window1

                            avg_energy2=0
                            do i=1,par%avg_window2
                                avg_energy2=avg_energy2+all_energy(it-i)
                            end do
                            avg_energy2=avg_energy2/par%avg_window2

                            sta_lta=avg_energy2/avg_energy1

                            if( abs(sta_lta -1) > par%sta_lta_trigger) then
                                write(*,*)
                                write(*,*) " !!!!!!!!!! warning, energy grows to much and pml is switched off !!!!!!"
                                write(*,*)
                                pmlcrit=.false.
                                pmlcheck(:)=.false.
                            end if
                        end if
                    end if
                end if
            end if
            !LL LL
            ! adjoint test
            if (par%save_forward) then
                if (mod(it,par%adj_nstep)==0) then
                    write(86, rec = isave ) ax
                    write(87, rec = isave ) az
                    write(88, rec = isave ) q(:,1)
                    write(89, rec = isave ) q(:,2)
                    write(90, rec = isave ) q(:,3)
                    isave = isave +1
                end if
            endif

            if (par%adjoint) then
                if (mod(it,par%adj_nstep)==0) then
                    read(86, rec = (isave)) adj_ax
                    read(87, rec = (isave)) adj_az
                    read(88, rec = (isave)) adj_sxx
                    read(89, rec = (isave)) adj_szz
                    read(90, rec = (isave)) adj_sxz
                    isave = isave - 1
                end if

                do ie=1,mesh%nelem
                    iv=mesh%ibool(:,ie)
                    adj_kernel_lambda(iv) = adj_kernel_lambda(iv) -&
                         ( ( ( q(iv,1) + q(iv,2) )*( adj_sxx(iv) + adj_szz(iv) ) ) / (4 * (mesh%lambda(ie) + mesh%mu(ie))**2) )*dt
                    adj_kernel_mu(iv) = adj_kernel_mu(iv) - &
                         ( (q(iv,3)*adj_sxz(iv))/mesh%mu(ie)**2 +&
                         1/4 * ( ( ( q(iv,1) + q(iv,2) )*( adj_sxx(iv) + adj_szz(iv) ) ) / (mesh%lambda(ie) + mesh%mu(ie))**2 +&
                         ( ( q(iv,1) - q(iv,2) )*( adj_sxx(iv) - adj_szz(iv) ) ) / mesh%mu(ie)**2 ) )*dt
                    adj_kernel_rho(iv) = adj_kernel_rho(iv) -&
                         ( ux(iv) * adj_ax(iv) + uz(iv) * adj_az(iv) ) *dt
                end do
            endif

            if (mod(it,movie%frame) == 0 .or. it == par%nt) then
                if (myrank == 0 ) then
                    call  system_clock(t2, rate,max) !get timer
                    dt1=real(t2 - t1, kind=8)/real(rate, kind=8)
                    dt_ges=dt_ges+dt1
                end if
                call maxval_real(maxval(uplot),maxu)
                call maxval_real(maxval(vplot),maxv)

                if (movie%movie) then
                    if (movie%velocity) then
                        call writePlotBinaries(vplot, "normV", myrank, it)
                        call writePlotBinaries(q(:,4), "vx", myrank, it)
                        call writePlotBinaries(q(:,5), "vz", myrank, it)
                    end if
                    if (movie%stress) then
                        call writePlotBinaries(q(:,1), "sigmaXX", myrank, it)
                        call writePlotBinaries(q(:,2), "sigmaZZ", myrank, it)
                        call writePlotBinaries(q(:,3), "sigmaXZ", myrank, it)
                    end if
                    if (movie%displacement) then
                        call writePlotBinaries(uplot, "normU", myrank, it)
                        call writePlotBinaries(ux, "ux", myrank, it)
                        call writePlotBinaries(uz, "uz", myrank, it)
                    end if
                endif

                if (myrank == 0) then
                    write(filename,"('out/energytemp')")
                    call plotPoints2D(stf(:),all_energy,filename)
                end if

                if (par%log .and. (myrank == 0)) then
                    write(*,*) "--------------------------------------------------------------------------------------"
                    write(*,*) " timestep :", it, "norm u :",maxu
                    write(*,*) movie%frame, "timesteps take ", dt1, " seconds and ", dt_ges," seconds so far."
                    write(*,*) "--------------------------------------------------------------------------------------"
                end if
                if (myrank == 0) then
                    call  system_clock(t1, rate,max) !get timer
                end if
            end if
        end do! timeloop

        if (par%adjoint) then
            do ie=1,mesh%nelem
                iv=mesh%ibool(:,ie)
                adj_kernel_vp(iv) = 2 * mesh%rho(ie) * mesh%vp(ie) * adj_kernel_lambda(iv)
                adj_kernel_vs(iv) = -4 * mesh%rho(ie) * mesh%vs(ie) * adj_kernel_lambda(iv) + 2 * mesh%rho(ie) * mesh%vs(ie) * adj_kernel_mu(iv)
                adj_kernel_rhop(iv) = (mesh%vp(ie)**2 - 2*mesh%vs(ie)**2) * adj_kernel_lambda(iv) + mesh%vs(ie)**2 * adj_kernel_mu(iv) + adj_kernel_rho(iv)
            end do
        endif

        write(filename,"('out/model',i6.6,'.dat')") myrank+1
        open(unit = 76, file = trim(filename),status = 'unknown')
        tmp = 0
        do ie=1,mesh%nelem ! element loop
            do j=1,Np
                tmp = tmp + 1
                i=mesh%ibool(j,ie)
                write(76,*)  tmp, mesh%vx(i), mesh%vz(i), mesh%rho(ie), mesh%vp(ie), mesh%vs(ie)
            end do
        end do
        close(76)

        if (par%save_forward) then
            close(86)
            close(87)
            close(88)
            close(89)
        end if

        if (par%adjoint) then
            close(86)
            close(87)
            close(88)
            close(89)

            write(filename,'(a,i6.6,a)') 'proc',myrank+1,'_rhop_alpha_beta_kernel.dat'
            open(unit = 86, file = 'out/'//filename,status = 'unknown')
            do ie=1,mesh%nelem ! element loop
                do j=1,Np
                    i=mesh%ibool(j,ie)
                    write(86,'(5e15.5e4)') mesh%vx(i), mesh%vz(i), adj_kernel_rhop(i), adj_kernel_vp(i), adj_kernel_vs(i)
                end do
            end do
            close(86)
            write(filename,'(a,i6.6,a)') 'proc',myrank+1,'_rho_lambda_mu_kernel.dat'
            open(unit = 86, file = 'out/'//filename,status = 'unknown')
            do ie=1,mesh%nelem ! element loop
                do j=1,Np
                    i=mesh%ibool(j,ie)
                    write(86,'(5e15.5e4)') mesh%vx(i), mesh%vz(i), adj_kernel_rho(i), adj_kernel_lambda(i), adj_kernel_mu(i)
                end do
            end do
            close(86)
        end if

        ! save seismos
        if (.not.par%adjoint) then
            if (mesh%has_rec) then
                do r=1,rec%nrec
                    write(filename,"('out/seismo.x.',i7.7,'.sdu')") rec%recnr(r)
                    call plotPoints2D(stf(:)-plott0,plotux(r,:),filename)
                    write(filename,"('out/seismo.z.',i7.7,'.sdu')") rec%recnr(r)
                    call plotPoints2D(stf(:)-plott0,plotuz(r,:),filename)
                    write(filename,"('out/seismo.x.',i7.7,'.sdv')") rec%recnr(r)
                    call plotPoints2D(stf(:)-plott0,plotv(r,:),filename)
                    write(filename,"('out/seismo.z.',i7.7,'.sdv')") rec%recnr(r)
                    call plotPoints2D(stf(:)-plott0,plotw(r,:),filename)
                    write(filename,"('out/seismo.x.',i7.7,'.sda')") rec%recnr(r)
                    call plotPoints2D(stf(:)-plott0,plotax(r,:),filename)
                    write(filename,"('out/seismo.z.',i7.7,'.sda')") rec%recnr(r)
                    call plotPoints2D(stf(:)-plott0,plotaz(r,:),filename)
                end do
            end if
        endif

        if (myrank==0) then
            write(filename,"('out/energy')")
            call plotPoints2D(stf(:),all_energy,filename)
        end if

        if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log.and.myrank==0) write(*,*) " done timeloop after", dt_ges," seconds"
        if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"

        close(27)

        deallocate(ux,uz,ax,az,uplot,vplot)
        deallocate(rQ,Q,Qn,Qm,e)
        deallocate(energy, all_energy, energy_kin, energy_pot)

        allocate(ux(mesh%nglob),uz(mesh%nglob),ax(mesh%nglob),az(mesh%nglob),uplot(mesh%nglob))
        allocate(rQ(mesh%nglob,5),Q(mesh%nglob,5),Qn(mesh%nglob,5),Qm(mesh%nglob,5))
        allocate(e(mesh%nglob,3))
        deallocate(a_rQ,theta,thetan,thetam)
        deallocate(stf)


        deallocate(srcelemV)
        deallocate(recelemV)
        if (.not.par%adjoint) then
            if (mesh%has_src) then
                deallocate(srcInt,srcTemp)
                deallocate(t0)
                deallocate(srcArray)
                deallocate(srcArrayM)
                deallocate(plotstf)
                deallocate(plotDiffstf)
            end if
            if (mesh%has_rec) then
                deallocate(recInt,recTemp)
                deallocate(plotv,plotw,plotux,plotuz,plotax,plotaz)
            end if
        else
            if (mesh%has_rec) then
                deallocate(recInt,recTemp)
                deallocate(adj_src_x,adj_src_z)
                deallocate(srcArray)
            end if
        endif

        deallocate(q_send)
        deallocate(q_rec)
        deallocate(qi)
        deallocate(qi_test)
        deallocate(req1)
        deallocate(pmlcheck)

        if(par%set_pml) then
            deallocate(fprime)
            deallocate(gprime)
            deallocate(fprimen,fprimem)
            deallocate(gprimen,gprimem)
            deallocate(ddx, ddz)
            deallocate(alphax, alphaz)
            deallocate(kx,kz)
        end if
        if (par%adjoint) then
            deallocate(adj_ax,adj_az)
            deallocate(adj_sxx,adj_szz, adj_sxz)
            deallocate(adj_kernel_lambda,adj_kernel_mu,adj_kernel_rho)
            deallocate(adj_kernel_vp,adj_kernel_vs,adj_kernel_rhop)
        endif
        ! stop MPI
        call finalize_mpi()
    end subroutine timeloop2d
end module timeloopMod
