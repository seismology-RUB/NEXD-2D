!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2019 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2019 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2019 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
!   Copyright 2016 Elena Busch (Rice University, USA)
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
    use slipInterfaceMod
    use rhsSlipInterface
    use rhsElasticMod
    use timestampMod
    use derMod
    use, intrinsic :: iso_fortran_env

    implicit none

    contains

    subroutine timeloop2d(par,mesh,src,rec,lsi, lsi_spec, lsipar, etolsi, movie, myrank, errmsg)
        type(error_message) :: errmsg
        type(parameterVar) :: par
        type(meshVar) :: mesh
        type(srcVar) ::src
        type(recVar) ::rec
        type(lsiVar), dimension(:) :: lsi
        type(interface_spec), dimension(:) :: lsi_spec
        type(elementToLSI), dimension(:) :: etolsi
        type(lsi_parameter) :: lsipar
        type(movie_parameter) :: movie
        ! time variables
        real(kind=CUSTOM_REAL) :: dt
        real(kind=CUSTOM_REAL) :: f0,f0tmp
        real(kind=CUSTOM_REAL), pointer, dimension(:) :: t0
        real(kind=CUSTOM_REAL) :: t0max, t0maxtmp
        real(kind=CUSTOM_REAL) :: time
        real(kind=CUSTOM_REAL) :: simt0
        real(kind=CUSTOM_REAL) :: fcrit,avg_energy1,avg_energy2,sta_lta
        integer :: timecrit
        logical :: pmlcrit = .true.
        ! free
        real(kind=custom_real), dimension(:,:), allocatable :: free
        ! indices
        integer :: i,j,r
        integer :: iglob
        integer, dimension(Np) :: iv
        integer :: ie, it, is
        ! RK
        integer :: irk, nrk, wrk
        real(kind=custom_real), dimension(:,:), allocatable :: resQ, resU !Runge-Kutta residual for 4th order rk
        ! rec variabels
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: recInt,recTemp
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: plotv,plotw,plotux,plotuz,plotax,plotaz,plotp1,plotp2,plotv1x,plotv1z,plotv2x,plotv2z,plot_r,plot_t
        real(kind=CUSTOM_REAL), dimension(1) :: r_v,s_v
        real(kind=CUSTOM_REAL) :: ux_temp, uz_temp, vx_temp, vz_temp, ax_temp, az_temp, angle_temp
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: energy,all_energy, energy_kin, energy_pot, np_zeros
        integer, dimension(:), pointer :: recelemv
        integer :: nsamples, isample
        ! source
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: srcInt,srcTemp
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotstf
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotDiffstf
        real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: srcArray
        real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: srcArrayM
        integer, dimension(:), pointer :: srcelemv
        real(kind=CUSTOM_REAL), dimension(2,2) :: M,Ms,Ma
        real(kind=CUSTOM_REAL) :: rk_time, stf_val, stf_val_diff
        ! PML
        logical, dimension(:), allocatable :: pmlcheck
        integer, dimension(:,:), allocatable :: pmlloc
        real(kind=CUSTOM_REAL) :: amax
        real(kind=CUSTOM_REAL) :: xmean, zmean
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ddx,ddz,alphax,alphaz,kx,kz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: fprime, gprime
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: resFprime, resGprime                 !runge Kutte 4th order residual storage
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: fprimen, fprimem, gprimen, gprimem

        ! fields
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ux,uz,ax,az, uplot, vplot, v1plot, v2plot
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rQ,Q,Qn,Qm, rQm,e,u
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ftemp,gtemp,htemp,qtemp
        real(kind=CUSTOM_REAL), dimension(Np) :: dFdr,dFds,dGdr,dGds
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: flux
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: invmass, vdmTinv
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: stf
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: anelasticvar
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: elasticfluxvar
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: APA, aAPA
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: T, invT, VT, VTfree, aVT, aVTfree, at

        !output
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: div
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: curl

        !logical :: movie
        logical :: file_exists
        character(len=80) ::filename

        ! usefulls
        real(kind=CUSTOM_REAL) :: onehalf = 1./2.
        real(kind=CUSTOM_REAL) :: onethree = 1./3.
        real(kind=CUSTOM_REAL) :: twothree = 2./3.
        real(kind=CUSTOM_REAL) :: onefor = 1./4.
        real(kind=CUSTOM_REAL) :: threefor = 3./4.

        ! timer
        type(timestampVar) :: timestamp
        real(kind = custom_real) :: localtime
        ! MPI
        integer :: myrank
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_send
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_rec
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi_test
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rq_send
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rq_rec
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rqi
        integer :: c,k
        integer :: dest
        integer :: req, req_r, tag, tag_r
        integer ,dimension(:), allocatable:: req1
        real(kind=CUSTOM_REAL) :: maxv
        real(kind=CUSTOM_REAL) :: maxu

        ! poroelasticity
        integer :: dimens

        ! attenuation
        real(kind=CUSTOM_REAL), dimension(Np,3*nMB) :: a_ftemp,a_gtemp
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: theta,thetam,thetan, resTheta, int_theta
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: a_rQ
        real(kind=CUSTOM_REAL), dimension(Ngll*3,3) :: aflux
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: a_sigma_tmp

        !LSI
        integer :: ilsi
        ! error message
        character(len=10) :: myname = "timeloop2d"



        if (par%autodt) then
            dt = mesh%dtfactor*par%cfl
        else
            dt = par%dt
        endif
        if (par%autont) then
            par%nt = int(par%t_total/dt) + 1
        endif

        call addTrace(errmsg, myname)

        if (mesh%has_src) then
            ! load sources
            write(filename,"('out/srcVar',i6.6)") myrank+1
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readSrcVar(src,filename)
            else
                call add(errmsg, 2, "File does not exist!", myname, filename)
                call print(errmsg)
                call stop_mpi()
            end if
            if (myrank == 0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            do i = 1, size(src%srcelem)
                write(*,'(a41, i5, a13, i5, a16)') "|                   Found source in rank ", myrank, " and element ", src%srcelem(i), "               |"
            enddo
            f0tmp = maxval(src%srcf0)
        else
            f0tmp = 0.0
        endif
        ! get the maximum f0
        call maxval_real_all(f0tmp,f0,CUSTOM_REAL)

        ! load receivers
        if (mesh%has_rec) then
            write(filename,"('out/recVar',i6.6)") myrank+1
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readRecVar(rec,filename)
            else
                call add(errmsg, 2, "Error in recVar, file does not exist!", myname, filename)
                call print(errmsg)
                call stop_mpi()
            end if
            if (myrank == 0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            do i = 1, size(rec%recelem)
                write(*,'(a41, i5, a13, i5, a16)') "|                 Found receiver in rank ", myrank, " and element ", rec%recelem(i), "               |"
            enddo
        end if

        call sync_mpi()
        if (myrank == 0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
        call sync_mpi()

        !if (par%debug) write(*,*) " Hello from task ",myrank+1, " / ", world_size, "with ", mesh%pinterfaces , " shared interfaces"
        ! allocate fields
        dimens = 5
        if (par%poroelastic) then
            dimens = dimens + 3*mesh%nfluids
        endif
        allocate(ux(mesh%nglob),uz(mesh%nglob),ax(mesh%nglob),az(mesh%nglob),uplot(mesh%nglob),vplot(mesh%nglob),v1plot(mesh%nglob),v2plot(mesh%nglob))
        allocate(rQ(mesh%nglob,dimens),Q(mesh%nglob,dimens),Qn(mesh%nglob,dimens),Qm(mesh%nglob,dimens), rQm(mesh%nglob, dimens), resQ(mesh%nglob,dimens), resU(mesh%nglob, 2), u(mesh%nglob, 2))
        allocate(ftemp(Np,dimens),gtemp(Np,dimens),htemp(Np,dimens),qtemp(Np,dimens))
        allocate(e(mesh%nglob,3))
        allocate(a_rQ(mesh%nglob,3*nMB),theta(mesh%nglob,3,nMB),thetan(mesh%nglob,3,nMB),thetam(mesh%nglob,3,nMB), resTheta(mesh%nglob,3,nMB))
        allocate(int_theta(mesh%nglob, 3, nMB))
        allocate(flux(3*NGLL,dimens))
        !allocate(plotstf(2,par%nt))
        allocate(stf(par%nt))
        allocate(energy(par%nt),energy_kin(par%nt), energy_pot(par%nt))
        allocate(all_energy(par%nt), a_sigma_tmp(Np, 3))
        energy     = 0.0
        energy_kin = 0.0
        energy_pot = 0.0
        all_energy = 0.0

        if (mesh%has_src) then
            allocate(srcInt(Np,src%nsrc),srcTemp(Np,src%nsrc))
            allocate(t0(src%nsrc))
            if (par%poroelastic) then
                    allocate(srcArray(Np,2+2*mesh%nfluids,src%nsrc))
                else
                    allocate(srcArray(Np,2,src%nsrc))
                endif
            allocate(srcArrayM(Np,3,src%nsrc))
            allocate(plotstf(par%nt,2,src%nsrc))
            allocate(plotDiffstf(par%nt,2,src%nsrc))
        end if
        if (mesh%has_rec) then
            allocate(recInt(Np,rec%nrec),recTemp(Np,rec%nrec))
            nsamples = int(par%nt/par%subsampling_factor)
            allocate(plotv(rec%nrec,nsamples), plotw(rec%nrec,nsamples))
            allocate(plotux(rec%nrec,nsamples), plotuz(rec%nrec,nsamples))
            allocate(plotax(rec%nrec,nsamples), plotaz(rec%nrec,nsamples))
            allocate(plot_r(rec%nrec,nsamples), plot_t(rec%nrec,nsamples))
            if (par%poroelastic .and. mesh%nfluids >= 1) then
                allocate(plotp1(rec%nrec,nsamples))
                allocate(plotv1x(rec%nrec,nsamples))
                allocate(plotv1z(rec%nrec,nsamples))
                if (mesh%nfluids == 2) then
                    allocate(plotp2(rec%nrec,nsamples))
                    allocate(plotv2x(rec%nrec,nsamples))
                    allocate(plotv2z(rec%nrec,nsamples))
                endif
            endif
        end if

        ! mpi
        allocate(q_send(NGLL*mesh%mpi_ne*dimens,mesh%mpi_nn))
        allocate(q_rec(NGLL*mesh%mpi_ne*dimens,mesh%mpi_nn))
        allocate(qi(NGLL,dimens,mesh%mpi_ne,mesh%mpi_nn))
        allocate(qi_test(NGLL,dimens,mesh%mpi_ne,mesh%mpi_nn))
        allocate(req1(mesh%mpi_nnmax))
        allocate(rq_send(NGLL*mesh%mpi_ne*dimens,mesh%mpi_nn))
        allocate(rq_rec(NGLL*mesh%mpi_ne*dimens,mesh%mpi_nn))
        allocate(rqi(NGLL,dimens,mesh%mpi_ne,mesh%mpi_nn))
        ! pml
        allocate(pmlcheck(mesh%nelem))
        pmlcheck=.false.

        if(par%set_pml) then
            allocate(fprime(mesh%nglob,dimens))
            allocate(gprime(mesh%nglob,dimens))
            allocate(resFprime(mesh%nglob,dimens))
            allocate(resGprime(mesh%nglob,dimens))
            allocate(fprimen(mesh%nglob,dimens),fprimem(mesh%nglob,dimens))
            allocate(gprimen(mesh%nglob,dimens),gprimem(mesh%nglob,dimens))
            allocate(np_zeros(mesh%nglob))
        end if
        if (mesh%has_rec) then
            plotux = 0.0
            plotuz = 0.0
            plotv  = 0.0
            plotw  = 0.0
            plotax = 0.0
            plotaz = 0.0
            if (par%poroelastic .and. mesh%nfluids >= 1) then
                plotp1  = 0.0
                plotv1x = 0.0
                plotv1z = 0.0
                if (mesh%nfluids == 2) then
                    plotp2  = 0.0
                    plotv2x = 0.0
                    plotv2z = 0.0
                endif
            endif
        end if

        q      = 1e-24!eps
        e      = 1e-24!eps
        rQ     = 1e-24!eps
        Qn     = 1e-24!eps
        resQ   = 1e-24!eps
        if (par%attenuation) then
            resTheta = 1e-24!eps !rk4
            a_rQ   = 1e-24!eps
            theta  = 1e-24!eps
            thetan = 1e-24!eps
            int_theta=1e-24!eps
        endif

        ux     = 0.0
        uz     = 0.0
        ax     = 0.0
        az     = 0.0
        uplot  = 0.0
        vplot  = 0.0
        v1plot = 0.0
        v2plot = 0.0

        if (par%use_trigger) then
            fcrit    = 5/f0
            timecrit = int(((1.2/f0)*3)/dt)
            ! sta lta properties for avoiding pml instabilities
            if (par%log.and.myrank==0) then
                write(*, '(a80)') "|                           sta-lta trigger properties                         |"
                write(*,'(a40, i7, a33)')   "|                             timecrit: ",timecrit, "                                |"
                write(*,'(a40, i7, a33)')   "|                         avg_window1 : ",par%avg_window1, "                                |"
                write(*,'(a40, i7, a33)')   "|                          avg_window2: ",par%avg_window2, "                                |"
                write(*,'(a40, f7.2, a33)') "|                      sta_lta_trigger: ",par%sta_lta_trigger, "                                |"
            end if
        endif

        nrk = 0
        wrk = 0
        if (par%timeint == 1) then
            if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            if (par%log.and.myrank==0) write(*,'(a80)') "|                                   use euler                                  |"
            if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            nrk = 1 !euler
            wrk = 1
        else if (par%timeint == 2) then
            if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            if (par%log.and.myrank==0) write(*,'(a80)') "|                       use tvd runge kutta second order                       |"
            if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            nrk = 2 ! rk2
            wrk = 1
        else if (par%timeint == 3) then
            if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            if (par%log.and.myrank==0) write(*,'(a80)') "|                       use tvd runge kutta third order                        |"
            if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            nrk = 3 ! rk3
            wrk = 2
        else if (par%timeint == 4) then
            if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            if (par%log.and.myrank==0) write(*,'(a80)') "|                        use runge kutta fourth order                          |"
            if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            nrk = 5 ! rk4
            wrk = 3
        end if

        if(par%set_pml) then
            if (par%log.and.myrank==0) then
                write(*,'(a80)') "|                        using NPML boundary conditions                        |"
                write(*,'(a40, f5.2, a35)') "|                                 xmin: ", mesh%pmlxmin, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 xmax: ", mesh%pmlxmax, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 zmin: ", mesh%pmlzmin, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 zmax: ", mesh%pmlzmax, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                            pml delta: ", par%pml_delta, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                   rc: ", par%pml_rc, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 kamx: ", par%pml_kmax, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 afac: ", par%pml_afac, "                                  |"
                write(*,'(a37,f8.2,a4,f8.2,a23)') '|                 x axis ranges from ', mesh%pmlxmin, ' to ', mesh%pmlxmax, ' |'
                write(*,'(a37,f8.2,a4,f8.2,a23)') '|                 z axis ranges from ', mesh%pmlzmin, ' to ', mesh%pmlzmax, ' |'
                write(*,"(a80)") "|------------------------------------------------------------------------------|"
            end if
            amax=par%pml_afac*pi*f0

            allocate(ddx(mesh%nglob), ddz(mesh%nglob))
            allocate(alphax(mesh%nglob), alphaz(mesh%nglob))
            allocate(kx(mesh%nglob),kz(mesh%nglob))
            allocate(pmlloc(mesh%nelem,2))
            ddx=0
            ddz=0
            kx=1
            kz=1
            alphax=0.
            alphaz=0.
            pmlloc=0
            pmlcheck=.false.
            do ie=1, mesh%nelem
                iv=mesh%ibool(:,ie)
                if (mesh%pml(ie)>0) then
                    pmlcheck(ie)=.true.      ! flag element as PML element
                    xmean=0.
                    zmean=0.
                    do j=1,NP               ! get center of element to determine its position
                        xmean=xmean+mesh%vx(iv(j))/Np
                        zmean=zmean+mesh%vz(iv(j))/Np
                    enddo
                    if (xmean<mesh%pmlxmin+par%pml_delta) pmlloc(ie,1)=-1     ! define to which PML (x,z) an element belongs (possibly multiple)
                    if (xmean>mesh%pmlxmax-par%pml_delta) pmlloc(ie,1)=1
                    if (zmean<mesh%pmlzmin+par%pml_delta) pmlloc(ie,2)=-1
                    if (zmean>mesh%pmlzmax+par%pml_delta) pmlloc(ie,2)=1
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(ie),mesh%pmlxmin,mesh%pmlxmax,mesh%pmlzmin,mesh%pmlzmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,pmlloc(ie,:))    ! get damping profiles
                endif
            enddo
            fprime = 0.
            gprime = 0.
            if (wrk == 3) then !rk4
                resFprime = 0.
                resGprime = 0.
            end if
        end if

        ! free surface conditions
        allocate(free(dimens,dimens))
        free=0.0
        free(1,1)=-2
        free(2,2)=0
        free(3,3)=-2
        free(4,4)=0
        free(5,5)=0
        if (par%poroelastic .and. mesh%nfluids >= 1) then
            free( 6, 6) = -2
            free( 7, 7) =  0
            free( 8, 8) =  0
            if (mesh%nfluids == 2) then
                free( 9, 9) = -2
                free(10,10) =  0
                free(11,11) =  0
            endif
        endif

        ftemp=0.
        gtemp=0.
        if (par%attenuation) then
            a_ftemp=0.
            a_gtemp=0.
        endif

        allocate(elasticfluxvar(mesh%nelem,4))
        allocate(APA(mesh%nelem,dimens,dimens))
        allocate(T(mesh%nelem,3,dimens,dimens), invT(mesh%nelem,3,dimens,dimens))
        allocate(VT(mesh%nelem,3,dimens,dimens), VTfree(mesh%nelem,3,dimens,dimens))
        if (par%attenuation) then
            allocate(anelasticvar(mesh%nelem*3*nMB))
            allocate(aVT(mesh%nelem,3,3,5), aVTfree(mesh%nelem,3,3,5),aT(mesh%nelem,3,3,3),aAPA(mesh%nelem,3,5))
        endif

        do ie=1,mesh%nelem
            elasticfluxvar(ie,1) = -(mesh%lambda(ie)+2*mesh%mu(ie))
            elasticfluxvar(ie,2) = -mesh%lambda(ie)
            elasticfluxvar(ie,3) = -mesh%mu(ie)
            elasticfluxvar(ie,4) = -1.0/mesh%rho(ie)
            APA(ie,:,:) = getAPA(mesh%vp(ie),mesh%vs(ie),mesh%rho(ie),mesh%lambda(ie),mesh%mu(ie))
            if (par%attenuation) aAPA(ie,:,:) = getAnelasticAPA(mesh%vp(ie),mesh%vs(ie),mesh%rho(ie))
            do is=1,3
                T(ie,is,:,:) =getT(mesh%nx(is*NGLL,ie),mesh%nz(is*NGLL,ie),dimens)
                invT(ie,is,:,:) =getinvT(mesh%nx(is*NGLL,ie),mesh%nz(is*NGLL,ie),dimens)
                VT(ie,is,:,:) = matmul(T(ie,is,:,:),matmul(APA(ie,:,:),invT(ie,is,:,:)))
                VTfree(ie,is,:,:) = matmul(T(ie,is,:,:),matmul(matmul(APA(ie,:,:),free),invT(ie,is,:,:)))
            enddo
            if (par%attenuation) then
                do is=1,3
                    aT(ie,is,:,:) = T(ie,is,1:3,1:3)
                    aVT(ie,is,:,:) = matmul(aT(ie,is,:,:),matmul(aAPA(ie,:,:),invT(ie,is,:,:)))
                    aVTfree(ie,is,:,:)= matmul(aT(ie,is,:,:),matmul(matmul(aAPA(ie,:,:),free),invT(ie,is,:,:)))
                enddo
            endif
        enddo
        if (par%attenuation) then
            do ie=1,mesh%nelem
                do i=1,nMB
                    anelasticvar(3*nMB*(ie-1)+3*(i-1)+1) = -mesh%lambda(ie)*mesh%ylambda(i,ie)
                    anelasticvar(3*nMB*(ie-1)+3*(i-1)+2) = -2*mesh%mu(ie)*mesh%ymu(i,ie)
                    anelasticvar(3*nMB*(ie-1)+3*(i-1)+3) = anelasticvar(3*nMB*(ie-1)+3*(i-1)+1)+anelasticvar(3*nMB*(ie-1)+3*(i-1)+2)
                enddo
            enddo
        endif

        vdmTinv=mesh%vdm
        vdmTinv=transpose(vdmTinv)
        call invert(vdmTinv, errmsg)
        allocate(srcelemV(mesh%nelem))
        allocate(recelemV(mesh%nelem))
        srcelemv=0
        recelemv=0

        ms = 0.
        if (mesh%has_src) then
            do i=1,src%nsrc
                r_v(1) = src%srcrs(1,i)
                s_v(1) = src%srcrs(2,i)
                call vdm2D(srcTemp(:,i),r_v,s_v)
                srcInt(:,i)=matmul(vdmTinv,srcTemp(:,i))

                ! set t0 to a save value
                t0(i)=1.2/src%srcf0(i)
                srcelemv(src%srcelem(i)) = i

                if (src%srctype(i) == 1) then  ! momenttensor
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
            t0maxtmp = maxval(t0)
        end if

        call sync_mpi()
        call maxval_real_all(t0maxtmp,t0max,CUSTOM_REAL)

        if (.not. par%shift_sources) then
            if (mesh%has_src) then
                do i=1,src%nsrc
                    t0(i)=0.
                enddo
            endif

            if (-t0max < par%simt0) then
                simt0 = -t0max
                if (myrank == 0) then
                    call add(errmsg, 1, "Source(s) started earlier than chosen simt0, so we changed simt0!", myname)
                end if
                if (par%log .and. myrank == 0) then
                    write(*,'(a80)') "|        Source(s) start earlier than chosen simt0, so we change simt0!        |"
                    write(*,'(a40, f10.7, a30)') "|                       original simt0: ", par%simt0, "                             |"
                    write(*,'(a40, f10.7, a30)') "|                            new simt0: ", simt0, "                             |"
                    write(*,'(a80)') "|------------------------------------------------------------------------------|"
                end if
            else
                simt0 = par%simt0
            endif
        else
            simt0 = par%simt0
        endif

        if (par%autoshift) then
            if (par%shift_sources) then
                ! This shifts the seismograms so that t=0 actually is the maximum of the wavelet with the smallest frequency (Ricker, Gaussian).
                par%plott0 = t0max ! 1.2/f0
            else
                par%plott0 = 0.
            end if
            ! else par%plott0 remains the same as provided in the parfile
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
        if (mesh%has_src) then
            do i=1,src%nsrc
                select case (src%srcstf(i))
                    case (1) !GAUSS
                        do it=1,par%nt
                            time = (float(it)-1.)*dt + simt0
                            select case (src%srctype(i))
                                case (0)  ! single force
                                    plotstf(it,1,i) = time
                                    plotstf(it,2,i) = -stfGauss(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                case (1)  ! moment tensor
                                    plotDiffstf(it,2,i) = -stfDiffGauss(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    plotDiffstf(it,1,i) = time
                            end select
                        enddo
                    case (2) !RICKER
                        do it=1,par%nt
                            time = (float(it)-1.)*dt + simt0
                            select case (src%srctype(i))
                                case (0)  ! single force
                                    plotstf(it,1,i) = time
                                    plotstf(it,2,i) = -stfRicker(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                case (1)  ! moment tensor
                                    plotDiffstf(it,1,i) = time
                                    plotDiffstf(it,2,i) = -stfDiffRicker(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                            end select
                        enddo
                    case (3) !SIN^3
                        do it=1,par%nt
                            time = (float(it)-1.)*dt + simt0
                            select case (src%srctype(i))
                                case (0)  ! single force
                                    plotstf(it,1,i) = time
                                    plotstf(it,2,i) = -stfSin3(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                case (1)  ! moment tensor
                                    plotDiffstf(it,1,i) = time
                                    plotDiffstf(it,2,i) = -stfDiffSin3(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                            end select
                        enddo
                    case (4) !EXTERNAL
                        filename=src%extwavelet(i)
                        select case (src%srctype(i))
                            case (0)  ! single force
                                call stfExternal(plotstf(:,2,i),plotstf(:,1,i),dt,par%nt,.true.,6,3,trim(filename),0, errmsg)
                            case (1)  ! moment tensor
                                call stfExternal(plotDiffstf(:,2,i),plotDiffstf(:,1,i),dt,par%nt,.true.,6,3,trim(filename),1, errmsg)
                        end select
                    case default
                        call add(errmsg, 2, "Chose a valid source time function. Available functions are listed in the 'source' parameter file.", myname, "data/source")
                        call print(errmsg)
                        call stop_mpi()
                end select

                ! write stf
                select case (src%srctype(i))
                    case (0)  ! single force
                        write(filename,"('out/stf',i6.6, '_rank', i7.7)") i, myrank+1
                        write(*,'(a40, a25, a15)')  "|     write sourcetimefunction to file: ", filename, "              |"
                        write(*,"(a80)") "|------------------------------------------------------------------------------|"
                        open(unit=27,file=trim(filename),status='unknown')
                        do it=1,par%nt
                            write(27,*) plotstf(it,1,i),plotstf(it,2,i)
                        end do
                        close(27)
                    case (1)  ! moment tensor
                        write(filename,"('out/stfdiff',i6.6, '_rank', i7.7)") i, myrank+1
                        write(*,'(a40, a30, a10)')  "|     write sourcetimefunction to file: ", filename, "         |"
                        write(*,"(a80)") "|------------------------------------------------------------------------------|"
                        open(unit=27,file=trim(filename),status='unknown')
                        do it=1,par%nt
                            write(27,*) plotDiffstf(it,1,i),plotDiffstf(it,2,i)
                        end do
                        close(27)
                end select
            end do
        end if

        if (lsipar%lsi .and. par%log .and. myrank == 0) then
            write(*,"(a80)") "|                                initialize LSI                                |"
            write(*,"(a80)") "|------------------------------------------------------------------------------|"
        end if

        ! mpi comunication
        ! build send buffer
        qm = q
        tag = 0
        q_send = 0.
        do i=1,mesh%mpi_nn
            do ie=1,mesh%mpi_ne ! loop over interface elements
                do k=1,dimens
                    do j=1,NGLL
                        if ( mesh%mpi_connection(i,ie,1) >0) then
                            q_send((ie-1)*dimens*NGLL + (k-1)*NGLL + j,i) = &
                                qm(mesh%ibt(j,mesh%mpi_connection(i,ie,2),mesh%mpi_connection(i,ie,1)),k)
                        end if
                    end do
                end do
            end do
        end do ! all interfaces

        !send and rec
        do i=1,mesh%mpi_nn
            dest=mesh%mpi_neighbor(i)-1
            call isendV_real(q_send(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag,req,CUSTOM_REAL)
            call irecV_real(q_rec(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag,req_r,CUSTOM_REAL)
            call wait_req(req)
            call wait_req(req_r)
        end do

        !unpack mpi buffers
        do i=1,mesh%mpi_nn
            c = 1
            do ie=1,mesh%mpi_ne
                do k=1,dimens
                    do j=1,NGLL
                        if ( mesh%mpi_connection(i,ie,1) > 0) then
                            qi(j,k,ie,i) = q_rec(c,i)
                        end if
                        c = c+1
                    end do
                end do
            end do
        end do

        if (lsipar%lsi) then !n > 0
            call initialSlipInterfaceActivity(mesh, lsi, lsi_spec, lsipar, etolsi, q, qi, errmsg)
        end if

        call sync_mpi()

        invmass=matmul(mesh%vdm,transpose(mesh%vdm))
        if (par%log.and.myrank==0) then
            call initialTimestamp(timestamp)
            write(*,"(a80)") "|                             starting timeloop...                             |"
            write(*, "(a40, es12.5, a28)") "|                                   dt: ", dt, "                           |"
            write(*, "(a40, i12, a28)") "|                                   nt: ", par%nt, "                           |"
            write(*, "(a40, es12.5, a28)") "|                 total simulated time: ", dt*par%nt, "                           |"
            write(*,"(a80)") "|------------------------------------------------------------------------------|"
        end if
        ! ---------------------------------------------------------------------------------------------
        ! ------------------------------------timeloop-------------------------------------------------
        ! ---------------------------------------------------------------------------------------------
        do it=1,par%nt ! timeloop
            if (myrank == 0) then
                localtime = (float(it)-1.) * dt + simt0
            end if
            rk_time=0.
            stf(it) = (float(it)-1.) * dt + simt0
            ! interpolate receivers
            if (mesh%has_rec .and. mod(it,par%subsampling_factor) == 0) then
                isample = int(it/par%subsampling_factor)
                do r=1,rec%nrec
                    iv = mesh%ibool(:,rec%recelem(r))
                    ! Calculate divergence and curl of the velocity component in order to create separate seismograms for
                    ! Radial and tangential component
                    if (par%curl) then
                        call curl2d(q(iv,4),q(iv,5),mesh%rx(iv),mesh%sx(iv),mesh%rz(iv),mesh%sz(iv),mesh%Dr,mesh%Ds,curl)
                        do j = 1, Np
                            plot_t(r, isample) = plot_t(r, isample) + curl(j,2)*recint(j,r)
                        end do
                    end if
                    if (par%div) then
                        call div2d(div,q(iv,4),q(iv,5),mesh%Dr,mesh%Ds,mesh%rx(iv),mesh%sx(iv),mesh%rz(iv),mesh%sz(iv))
                        do j = 1, Np
                            plot_r(r, isample) = plot_r(r, isample) + div(j)*recint(j,r)
                        end do
                    end if

                    do j=1,Np
                        iglob=mesh%ibool(j,rec%recelem(r))
                        plotux(r,isample) = plotux(r,isample) + ux(iglob)  * recint(j,r)
                        plotuz(r,isample) = plotuz(r,isample) + uz(iglob)  * recint(j,r)
                        plotv(r,isample)  = plotv(r,isample)  + q(iglob,4) * recint(j,r)
                        plotw(r,isample)  = plotw(r,isample)  + q(iglob,5) * recint(j,r)
                        plotax(r,isample) = plotax(r,isample) + ax(iglob)  * recint(j,r)
                        plotaz(r,isample) = plotaz(r,isample) + az(iglob)  * recint(j,r)
                        if (par%poroelastic .and. mesh%nfluids >= 1) then
                            plotp1(r,isample)  = plotp1(r,isample) + q(iglob,6) * recint(j,r)
                            plotv1x(r,isample) = plotv1x(r,isample) + q(iglob,7) * recint(j,r)
                            plotv1z(r,isample) = plotv1z(r,isample) + q(iglob,8) * recint(j,r)
                            if (mesh%nfluids == 2) then
                                plotp2(r,isample)  = plotp2(r,isample) + q(iglob,9) * recint(j,r)
                                plotv1x(r,isample) = plotv1x(r,isample) + q(iglob,7) * recint(j,r)
                                plotv1z(r,isample) = plotv1z(r,isample) + q(iglob,8) * recint(j,r)
                            endif
                        endif
                    end do
                    ! rotate receivers
                    angle_temp = (rec%rec_ang(r)*PI)/180.
                    ux_temp = cos(angle_temp) * plotux(r,isample) + sin(angle_temp) * plotuz(r,isample)
                    uz_temp = -sin(angle_temp) * plotux(r,isample) + cos(angle_temp) * plotuz(r,isample)
                    vx_temp = cos(angle_temp) * plotv(r,isample) + sin(angle_temp) * plotw(r,isample)
                    vz_temp = -sin(angle_temp) * plotv(r,isample) + cos(angle_temp) * plotw(r,isample)
                    ax_temp = cos(angle_temp) * plotax(r,isample) + sin(angle_temp) * plotaz(r,isample)
                    az_temp = -sin(angle_temp) * plotax(r,isample) + cos(angle_temp) * plotaz(r,isample)
                    plotux(r,isample) = ux_temp
                    plotuz(r,isample) = uz_temp
                    plotv(r,isample) = vx_temp
                    plotw(r,isample) = vz_temp
                    plotax(r,isample) = ax_temp
                    plotaz(r,isample) = az_temp
                end do
            end if

            tag   = 0
            tag_r = 0
            stf_val = 0.
            stf_val_diff = 0.
            do irk=1,nrk ! runge-kutta loop
                time = (float(it) - 1. + rk_time) * dt + simt0
                if (mesh%has_src) then
                    do i=1, src%nsrc
                        select case (src%srctype(i))
                            case (0)  ! single force
                                select case (src%srcstf(i))
                                    case (1)  ! GAUSS
                                        stf_val = -stfGauss(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    case (2)  ! RICKER
                                        stf_val = -stfRicker(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    case (3)  ! SIN^3
                                        stf_val = -stfSin3(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    case (4)  ! EXTERNAL
                                        if (it<par%nt) then
                                            stf_val=(1.-rk_time)*plotstf(it,2,i)+rk_time*plotstf(it+1,2,i)
                                        else
                                            stf_val=(1.+rk_time)*plotstf(it,2,i)-rk_time*plotstf(it-1,2,i)
                                        endif
                                end select

                                if (par%poroelastic) then
                                    srcArray(:,1,i)=srcInt(:,i)*(-sin((src%srcangle_force(i)*PI)/180.)*stf_val)*(-mesh%A(4,1,src%srcelem(i)))
                                    srcArray(:,2,i)=srcInt(:,i)*(+cos((src%srcangle_force(i)*PI)/180.)*stf_val)*(-mesh%B(5,2,src%srcelem(i)))
                                    if (mesh%nfluids >= 1) then
                                        srcArray(:,3,i)=srcInt(:,i)*(-sin((src%srcangle_force(i)*PI)/180.)*stf_val)*(-mesh%A(7,1,src%srcelem(i)))
                                        srcArray(:,4,i)=srcInt(:,i)*(+cos((src%srcangle_force(i)*PI)/180.)*stf_val)*(-mesh%B(8,2,src%srcelem(i)))
                                        if (mesh%nfluids == 2) then
                                            srcArray(:,5,i)=srcInt(:,i)*(-sin((src%srcangle_force(i)*PI)/180.)*stf_val)*(-mesh%A(10,1,src%srcelem(i)))
                                            srcArray(:,6,i)=srcInt(:,i)*(+cos((src%srcangle_force(i)*PI)/180.)*stf_val)*(-mesh%B(11,2,src%srcelem(i)))
                                        endif
                                    endif
                                else
                                    srcArray(:,1,i)=srcInt(:,i)*(-sin((src%srcangle_force(i)*PI)/180.)*stf_val)/mesh%rho(src%srcelem(i))
                                    srcArray(:,2,i)=srcInt(:,i)*(+cos((src%srcangle_force(i)*PI)/180.)*stf_val)/mesh%rho(src%srcelem(i))
                                endif

                                iv=mesh%ibool(:,src%srcelem(i))
                                srcArray(:,1,i)=matmul(invmass,srcArray(:,1,i))/mesh%jacobian(iv)
                                srcArray(:,2,i)=matmul(invmass,srcArray(:,2,i))/mesh%jacobian(iv)
                                if (par%poroelastic .and. mesh%nfluids >= 1) then
                                    srcArray(:,3,i)=matmul(invmass,srcArray(:,3,i))/mesh%jacobian(iv)
                                    srcArray(:,4,i)=matmul(invmass,srcArray(:,4,i))/mesh%jacobian(iv)
                                    if (mesh%nfluids == 2) then
                                        srcArray(:,5,i)=matmul(invmass,srcArray(:,5,i))/mesh%jacobian(iv)
                                        srcArray(:,6,i)=matmul(invmass,srcArray(:,6,i))/mesh%jacobian(iv)
                                    endif
                                endif

                            case (1)  ! moment tensor
                                select case (src%srcstf(i))
                                    case (1)  ! GAUSS
                                        stf_val_diff = -stfDiffGauss(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    case (2)  ! RICKER
                                        stf_val_diff = -stfDiffRicker(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    case (3)  ! SIN^3
                                        stf_val_diff = -stfDiffSin3(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    case (4)  ! EXTERNAL
                                        if (it<par%nt) then
                                            stf_val_diff=(1.-rk_time)*plotDiffstf(it,2,i)+rk_time*plotDiffstf(it+1,2,i)
                                        else
                                            stf_val_diff=(1.+rk_time)*plotDiffstf(it,2,i)-rk_time*plotDiffstf(it-1,2,i)
                                        endif
                                end select

                                srcArrayM(:,1,i)=srcInt(:,i)*(ms(1,1)*stf_val_diff)
                                srcArrayM(:,2,i)=srcInt(:,i)*(ms(2,2)*stf_val_diff)
                                srcArrayM(:,3,i)=srcInt(:,i)*(ms(1,2)*stf_val_diff)

                                iv=mesh%ibool(:,src%srcelem(i))
                                srcArrayM(:,1,i)=matmul(invmass,srcArrayM(:,1,i))/mesh%jacobian(iv)
                                srcArrayM(:,2,i)=matmul(invmass,srcArrayM(:,2,i))/mesh%jacobian(iv)
                                srcArrayM(:,3,i)=matmul(invmass,srcArrayM(:,3,i))/mesh%jacobian(iv)
                        end select
                    end do
                end if
                qm=q
                rQm = rQ
                if (par%attenuation) then
                    thetam=theta
                end if
                if(par%set_pml) then
                    fprimem=fprime
                    gprimem=gprime
                end if
                q_send(:,:) = 0
                rq_send(:,:) = 0
                ! mpi comunication
                ! build send buffer
                do i=1,mesh%mpi_nn
                    do ie=1,mesh%mpi_ne ! loop over interface elements
                        do k=1,dimens
                            do j=1,NGLL
                                if ( mesh%mpi_connection(i,ie,1) >0) then
                                    q_send((ie-1)*dimens*NGLL + (k-1)*NGLL + j,i) = &
                                        qm(mesh%ibt(j,mesh%mpi_connection(i,ie,2),mesh%mpi_connection(i,ie,1)),k)
                                    rq_send((ie-1)*dimens*NGLL + (k-1)*NGLL + j,i) = &
                                        rqm(mesh%ibt(j,mesh%mpi_connection(i,ie,2),mesh%mpi_connection(i,ie,1)),k)
                                end if
                            end do
                        end do
                    end do
                end do ! all interfaces

                ! send and rec
                do i=1,mesh%mpi_nn
                    dest=mesh%mpi_neighbor(i)-1
                    call isendV_real(q_send(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag,req,CUSTOM_REAL)
                    call irecV_real(q_rec(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag,req_r,CUSTOM_REAL)
                    call wait_req(req)
                    call wait_req(req_r)
                    call isendV_real(rq_send(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag_r,req,CUSTOM_REAL)
                    call irecV_real(rq_rec(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag_r,req_r,CUSTOM_REAL)
                    call wait_req(req)
                    call wait_req(req_r)
                end do

                ! unpack mpi buffers
                do i=1,mesh%mpi_nn
                    c = 1
                    do ie=1,mesh%mpi_ne
                        do k=1,dimens
                            do j=1,NGLL
                                if ( mesh%mpi_connection(i,ie,1) > 0) then
                                    qi(j,k,ie,i) = q_rec(c,i)
                                    rqi(j,k,ie,i) = rq_rec(c,i)
                                end if
                                c = c+1
                            end do
                        end do
                    end do
                end do

                !Start calculation of the right-hand-side of the dg equation
                if (lsipar%lsi) then !n > 0
                    ilsi = 1
                else
                    ilsi = 0
                end if
                do ie=1,mesh%nelem ! element loop
                    iv=mesh%ibool(:,ie)
                    select case(par%fluxtype)
                        case (0)
                            call rhsElastic(ie, iv, ilsi, mesh, elasticfluxvar(ie,:), lsi, lsipar, etolsi, free, qm, qi, rQ, ftemp, gtemp, fprimem, gprimem, pmlcheck, kx, kz)
                        case (1)
                            qtemp=qm(iv,:)

                            ! get fluxes
                            if (par%poroelastic) then
                                call poroelasticFluxes(qtemp,mesh%A(:,:,ie),mesh%B(:,:,ie),mesh%E(:,:,ie),ftemp,gtemp,htemp)
                            else
                                call elasticFluxes(qtemp,elasticfluxvar(ie,:),ftemp,gtemp)
                            endif

                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    ftemp(:,i) = (fprimem(iv,i) + ftemp(:,i))/kx(iv)
                                    gtemp(:,i) = (gprimem(iv,i) + gtemp(:,i))/kz(iv)
                                enddo
                            end if
                            ! comp strong deriverate
                            do i=1,dimens
                                dFdr = matmul(mesh%Dr,ftemp(:,i))
                                dFds = matmul(mesh%Ds,ftemp(:,i))
                                dGdr = matmul(mesh%Dr,gtemp(:,i))
                                dGds = matmul(mesh%Ds,gtemp(:,i))
                                rQ(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds) + (mesh%rz(iv)*dGdr + mesh%sz(iv)*dGds)
                                if (par%poroelastic) then
                                    rQ(iv,i) = rQ(iv,i) + htemp(:,i)
                                endif
                            end do

                            if (par%attenuation) then
                                call anelasticFluxes(qtemp,mesh%wl(:,ie),a_ftemp,a_gtemp)
                                do i=1,3*nMB
                                    dFdr = matmul(mesh%Dr,a_ftemp(:,i))
                                    dFds = matmul(mesh%Ds,a_ftemp(:,i))
                                    dGdr = matmul(mesh%Dr,a_gtemp(:,i))
                                    dGds = matmul(mesh%Ds,a_gtemp(:,i))
                                    a_rQ(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds) + (mesh%rz(iv)*dGdr + mesh%sz(iv)*dGds)
                                end do
                            end if

                            ! compute fluxes on the surface
                            if (par%poroelastic) then
                                call computeExactRiemannSF(flux,qm,qi,mesh%neighbor(:,ie),mesh%AP(:,:,ie),&
                                    mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie),&
                                    mesh%ibt(:,:,ie),mesh%ibn(:,:,ie),mesh%nx(:,ie),mesh%nz(:,ie),free,dimens)
                            else
                                call computeExactRiemannSF(flux,qm,qi,mesh%neighbor(:,ie),VT(ie,:,:,:), VTfree(ie,:,:,:),&
                                    mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie),&
                                    mesh%ibt(:,:,ie),mesh%ibn(:,:,ie),5)
                            endif

                            if(par%attenuation) then
                                call computeExactRiemannSF(aflux,qm,qi,mesh%neighbor(:,ie),aVT(ie,:,:,:), aVTfree(ie,:,:,:),&
                                mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie),&
                                mesh%ibt(:,:,ie),mesh%ibn(:,:,ie),3)
                                c=1
                                do i=1,nMB
                                    do j=1,3
                                        a_rQ(iv,c) = a_rQ(iv,c) + matmul(mesh%lift,mesh%fscale(:,ie)*mesh%wl(i,ie)*aflux(:,j)/2.0)
                                        a_rQ(iv,c) = a_rQ(iv,c) - mesh%wl(i,ie) * thetam(iv,j,i)
                                        c=c+1
                                    end do
                                end do

                                do j=1,nMB
                                    rq(iv,1) = rq(iv,1) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+3) * thetam(iv,1,j) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+1) * thetam(iv,2,j)
                                    rq(iv,2) = rq(iv,2) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+1) * thetam(iv,1,j) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+3) * thetam(iv,2,j)
                                    rq(iv,3) = rq(iv,3) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+2) * thetam(iv,3,j)
                                end do
                            end if
                            ! lift surface integral
                            do i=1,dimens
                                rq(iv,i) = rq(iv,i) + matmul(mesh%lift,mesh%fscale(:,ie)*flux(:,i)/2.0)
                            end do
                            !End of calculation of the rhs
                        case default
                            call add(errmsg,2,'No proper fluxtype entered!',myname)
                            return
                    end select

                    do i=1,src%nsrc
                        if (src%srcelem(i)==ie) then      ! add source term
                            if (src%srctype(i) == 0) then
                                !single force
                                rq(iv,4:5)=rq(iv,4:5)+srcArray(:,1:2,i)
                                if (par%poroelastic .and. mesh%nfluids >= 1) then
                                    rq(iv,7:8)=rq(iv,7:8)+srcArray(:,3:4,i)
                                    if (mesh%nfluids == 2) then
                                        rq(iv,10:11)=rq(iv,10:11)+srcArray(:,5:6,i)
                                    endif
                                endif
                            elseif (src%srctype(i) == 1) then
                                rq(iv,1:3)=rq(iv,1:3)+srcArrayM(:,1:3,i)
                            endif
                        endif
                    enddo

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
                            if (lsipar%lsi) then !n > 0
                                if (etolsi(ie)%isLSI) then
                                    lsi(ilsi)%resS_n = lsi(ilsi)%S_n
                                    lsi(ilsi)%resS_t = lsi(ilsi)%S_t
                                endif
                            endif

                            q(iv,:) = q(iv,:) + dt*rq(iv,:)

                            !pml
                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    fprime(iv,i) = fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
                                    gprime(iv,i) = gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
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
                            if (lsipar%lsi) then !n > 0
                                if (etolsi(ie)%isLSI) then
                                    call activityFluxSlipInterface(ie, mesh, lsi(ilsi), lsipar, lsi_spec, etolsi, rQm, rQi)
                                    if (lsipar%normal) then
                                        lsi(ilsi)%S_n = lsi(ilsi)%S_n + dt*lsi(ilsi)%rhsS_n
                                    end if
                                    if (lsipar%tangential) then
                                        lsi(ilsi)%S_t = lsi(ilsi)%S_t + dt*lsi(ilsi)%rhsS_t
                                    end if
                                    ilsi = ilsi + 1
                                end if
                            endif
                            rk_time=1.
                        end if !rk==1
                        if (irk==2) then
                            q(iv,1:dimens) = onehalf*(qn(iv,:) + q(iv,:) + dt*rq(iv,:))

                            !pml
                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    fprime(iv,i) = onehalf* (fprimen(iv,i) + fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i))))
                                    gprime(iv,i) = onehalf* (gprimen(iv,i) + gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i))))
                                end do
                            end if
                            ! attenuation
                            if(par%attenuation) then
                                c = 1
                                do i=1,nMB
                                    do j=1,3
                                        theta(iv,j,i) = onehalf*(thetan(iv,j,i) + theta(iv,j,i) + dt*(a_rQ(iv,c)))
                                        c = c + 1
                                    end do
                                end do
                            end if
                            !LSI
                            if (lsipar%lsi) then !n > 0
                                if (etolsi(ie)%isLSI) then
                                    call activityFluxSlipInterface(ie, mesh, lsi(ilsi), lsipar, lsi_spec, etolsi, rQm, rQi)
                                    if (lsipar%normal) then
                                        lsi(ilsi)%S_n = onehalf*(lsi(ilsi)%resS_n + lsi(ilsi)%S_n + dt*lsi(ilsi)%rhsS_n)
                                    end if
                                    if (lsipar%tangential) then
                                        lsi(ilsi)%S_t = onehalf*(lsi(ilsi)%resS_t + lsi(ilsi)%S_t + dt*lsi(ilsi)%rhsS_t)
                                    end if
                                    ilsi = ilsi + 1
                                end if
                            endif
                        end if !rk==2
                    !RK 3
                    else if (wrk==2) then
                        if (irk==1) then
                            qn(iv,:)=q(iv,:)
                            if(par%set_pml) then
                                fprimen(iv,:)=fprime(iv,:)
                                gprimen(iv,:)=gprime(iv,:)
                            end if
                            if(par%attenuation) then
                                thetan(iv,:,:) = theta(iv,:,:)
                            end if
                            if (lsipar%lsi) then
                                if (etolsi(ie)%isLSI) then
                                    lsi(ilsi)%resS_n = lsi(ilsi)%S_n
                                    lsi(ilsi)%resS_t = lsi(ilsi)%S_t
                                endif
                            endif

                            q(iv,:) = q(iv,:) + dt*rq(iv,:)
                            !PML
                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    fprime(iv,i) = fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
                                    gprime(iv,i) = gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
                                end do
                            end if
                            !Attenuation
                            if(par%attenuation) then
                                c = 1
                                do i=1,nMB
                                    do j=1,3
                                        theta(iv,j,i) = theta(iv,j,i) + dt*a_rQ(iv,c)
                                        c = c + 1
                                    end do
                                end do
                            end if
                            !LSI
                            if (lsipar%lsi) then !n > 0
                                if (etolsi(ie)%isLSI) then
                                    call activityFluxSlipInterface(ie, mesh, lsi(ilsi), lsipar, lsi_spec, etolsi, rQm, rQi)
                                    if (lsipar%normal) then
                                        lsi(ilsi)%S_n = lsi(ilsi)%S_n + dt*lsi(ilsi)%rhsS_n
                                    end if
                                    if (lsipar%tangential) then
                                        lsi(ilsi)%S_t = lsi(ilsi)%S_t + dt*lsi(ilsi)%rhsS_t
                                    end if
                                    ilsi = ilsi + 1
                                end if
                            endif
                            rk_time=1.
                        end if !rk==1
                        if (irk==2) then
                            q(iv,:) = threefor*qn(iv,:) + onefor*(q(iv,:) + dt*rq(iv,:))
                            !PML
                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    fprime(iv,i) = threefor*fprimen(iv,i) + onefor*(fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i))))
                                    gprime(iv,i) = threefor*gprimen(iv,i) + onefor*(gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i))))
                                end do
                            end if
                            !Attenuation
                            if(par%attenuation) then
                                c = 1
                                do i=1,nMB
                                    do j=1,3
                                        theta(iv,j,i) = threefor*thetan(iv,j,i) + onefor*(theta(iv,j,i) + dt*a_rQ(iv,c))
                                        c = c + 1
                                    end do
                                end do
                            end if
                            !LSI
                            if (lsipar%lsi) then
                                if (etolsi(ie)%isLSI) then
                                    call activityFluxSlipInterface(ie, mesh, lsi(ilsi), lsipar, lsi_spec, etolsi, rQm, rQi)
                                    if (lsipar%normal) then
                                        lsi(ilsi)%S_n = threefor*lsi(ilsi)%resS_n + onefor*(lsi(ilsi)%S_n + dt*lsi(ilsi)%rhsS_n)
                                    end if
                                    if (lsipar%tangential) then
                                        lsi(ilsi)%S_t = threefor*lsi(ilsi)%resS_t + onefor*(lsi(ilsi)%S_t + dt*lsi(ilsi)%rhsS_t)
                                    end if
                                    ilsi = ilsi + 1
                                end if
                            endif
                            rk_time=0.5
                        end if !rk==2
                        if (irk==3) then
                            q(iv,:) = onethree*qn(iv,:) + twothree*(q(iv,:) + dt*rq(iv,:))
                            !PML
                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    fprime(iv,i) = onethree*fprimen(iv,i) + twothree*(fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i))))
                                    gprime(iv,i) = onethree*gprimen(iv,i) + twothree*(gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i))))
                                end do
                            end if
                            !Attenuation
                            if(par%attenuation) then
                                c = 1
                                do i=1,nMB
                                    do j=1,3
                                        theta(iv,j,i) = onethree*thetan(iv,j,i) + twothree*(theta(iv,j,i) + dt*a_rQ(iv,c))
                                        c = c + 1
                                    end do
                                end do
                            end if
                            !LSI
                            if (lsipar%lsi) then
                                if (etolsi(ie)%isLSI) then
                                    call activityFluxSlipInterface(ie, mesh, lsi(ilsi), lsipar, lsi_spec, etolsi, rQm, rQi)
                                    if (lsipar%normal) then
                                        lsi(ilsi)%S_n = onethree*lsi(ilsi)%resS_n + twothree*(lsi(ilsi)%S_n + dt*lsi(ilsi)%rhsS_n)
                                    end if
                                    if (lsipar%tangential) then
                                        lsi(ilsi)%S_t = onethree*lsi(ilsi)%resS_t + twothree*(lsi(ilsi)%S_t + dt*lsi(ilsi)%rhsS_t)
                                    end if
                                    ilsi = ilsi + 1
                                end if
                            endif
                        end if !rk==3
                    ! RK 4
                    else if (wrk == 3) then
                        resQ(iv, :) = rk4a(irk)*resQ(iv, :) + dt*rq(iv,:)
                        q(iv, :) = q(iv, :) + rk4b(irk)*resQ(iv, :)
                        !pml
                        if (pmlcheck(ie)) then
                            do i=1,dimens
                                resFprime(iv,i) = rk4a(irk)*resFprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
                                resGprime(iv,i) = rk4a(irk)*resGprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
                                fprime(iv,i)    = fprime(iv,i) + rk4b(irk)*resFprime(iv,i)
                                gprime(iv,i)    = gprime(iv,i) + rk4b(irk)*resGprime(iv,i)
                            end do
                        end if
                        ! attenuation
                        if(par%attenuation) then
                            c = 1
                            do i=1,nMB
                                do j=1,3
                                    resTheta(iv, j, i) = rk4a(irk)*resTheta(iv, j, i) + dt*a_rQ(iv,c)
                                    theta(iv, j, i)    = theta(iv, j, i) + rk4b(irk)*resTheta(iv, j, i)
                                    c = c + 1
                                end do
                            end do
                        end if
                        if (lsipar%lsi) then !n > 0
                            if (etolsi(ie)%isLSI) then
                                call activityFluxSlipInterface(ie, mesh, lsi(ilsi), lsipar, lsi_spec, etolsi, rQm, rQi)
                                if (lsipar%normal) then
                                    lsi(ilsi)%resS_n = rk4a(irk)*lsi(ilsi)%resS_n + dt*lsi(ilsi)%rhsS_n
                                    lsi(ilsi)%S_n    = lsi(ilsi)%S_n + rk4b(irk)*lsi(ilsi)%resS_n
                                end if
                                if (lsipar%tangential) then
                                    lsi(ilsi)%resS_t = rk4a(irk)*lsi(ilsi)%resS_t + dt*lsi(ilsi)%rhsS_t
                                    lsi(ilsi)%S_t    = lsi(ilsi)%S_t + rk4b(irk)*lsi(ilsi)%resS_t
                                end if
                                ilsi = ilsi + 1
                            end if
                        endif
                        if (irk < nrk) rk_time=rk4c(irk+1)
                    end if !end rk4
                    if (irk == nrk) then
                        !integrate fields with euler, for plotting
                        !This is done inside of the rk-loop so it is done for each iteration of the rk-loop.
                        !In case of rk2 it is done twice and in case of rk4 it would be done 4 times.
                        !Displacement
                        ux(iv) = ux(iv) + dt*q(iv,4)
                        uz(iv) = uz(iv) + dt*q(iv,5)
                        uplot(iv) = sqrt(ux(iv)**2+uz(iv)**2)
                        if (movie%movie) then
                            if (movie%velocity) then
                                !Calculate the norm of the velocity
                                vplot(iv) = sqrt(q(iv,4)**2+q(iv,5)**2)
                            endif
                            if (par%poroelastic .and. mesh%nfluids >= 1 .and. movie%v1) then
                                !Calculate the norm of the velocity
                                v1plot(iv) = sqrt(q(iv,7)**2+q(iv,8)**2)
                            endif
                            if (par%poroelastic .and. mesh%nfluids == 2 .and. movie%v2) then
                                !Calculate the norm of the velocity
                                v2plot(iv) = sqrt(q(iv,10)**2+q(iv,11)**2)
                            endif
                        endif

                        !Acceleration
                        ax(iv) = rq(iv,4)
                        az(iv) = rq(iv,5)
                        ! Strain
                        if (par%attenuation) then
                            do i=1,nMB
                                int_theta(iv,:,i)= int_theta(iv,:,i) + dt * theta(iv,:,i)
                            enddo
                        endif
                        a_sigma_tmp(:,1)=q(iv,1)
                        a_sigma_tmp(:,2)=q(iv,2)
                        a_sigma_tmp(:,3)=q(iv,3)

                        if (par%attenuation) then
                            do i=1,nMB
                                a_sigma_tmp(:,1)=a_sigma_tmp(:,1)+mesh%lambda(ie)*mesh%ylambda(i,ie)*(int_theta(iv,1,i)+int_theta(iv,2,i))+2*mesh%mu(ie)*mesh%ymu(i,ie)*int_theta(iv,1,i)
                                a_sigma_tmp(:,2)=a_sigma_tmp(:,2)+mesh%lambda(ie)*mesh%ylambda(i,ie)*(int_theta(iv,1,i)+int_theta(iv,2,i))+2*mesh%mu(ie)*mesh%ymu(i,ie)*int_theta(iv,2,i)
                                a_sigma_tmp(:,3)=a_sigma_tmp(:,3)+2*mesh%mu(ie)*mesh%ymu(i,ie)*int_theta(iv,3,i)
                            enddo
                        endif

                        if (.not. par%poroelastic) then
                            e(iv,1) = (-mesh%lambda(ie) * a_sigma_tmp(:,2) + (2 * mesh%mu(ie) + mesh%lambda(ie)) * a_sigma_tmp(:,1))/ (4*mesh%mu(ie)**2 + 4 * mesh%lambda(ie)*mesh%mu(ie))
                            e(iv,2) = ((2 * mesh%mu(ie) + mesh%lambda(ie)) * a_sigma_tmp(:,2) - mesh%lambda(ie)*a_sigma_tmp(:,1))/(4*mesh%mu(ie)**2 + 4*mesh%lambda(ie)*mesh%mu(ie))
                            e(iv,3) = a_sigma_tmp(:,3)/(2 * mesh%mu(ie))
                        endif
                    endif
                end do !element loop
                call sync_mpi()
            end do ! RK-LOOP

            do ie=1,mesh%nelem
                iv=mesh%ibool(:,ie)
                do i=1,Np
                    if (par%poroelastic) then
                        energy_kin(it) = energy_kin(it) + 1./mesh%A(4,1,ie) * (q(iv(i),4)**2 + q(iv(i),5)**2)
                    else
                        energy_kin(it) = energy_kin(it) + mesh%rho(ie) * (q(iv(i),4)**2 + q(iv(i),5)**2)
                    endif
                    energy_pot(it) = energy_pot(it) + 0.5 * (q(iv(i),1)*e(iv(i),1) + q(iv(i),2)*e(iv(i),2) + q(iv(i),3)*e(iv(i),3))
                end do
            end do
            energy(it) = energy_kin(it) + energy_pot(it)
            call sum_real_all(energy(it),all_energy(it), CUSTOM_REAL)

            ! check if energy grows to much, its done with sta/lta trigger to switch off the pml if necessary
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

            if (mod(it,movie%frame) == 0 .or. it == par%nt) then
                call maxval_real(maxval(uplot),maxu,CUSTOM_REAL)
                call maxval_real(maxval(vplot),maxv,CUSTOM_REAL)
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
                    if (par%poroelastic .and. mesh%nfluids >= 1) then
                        if (movie%p1) then
                            call writePlotBinaries(q(:,6), "p1", myrank, it)
                        endif
                        if (movie%v1) then
                            call writePlotBinaries(v1plot, "normv1", myrank, it)
                            call writePlotBinaries(q(:,7), "v1x", myrank, it)
                            call writePlotBinaries(q(:,8), "v1z", myrank, it)
                        endif
                        if (mesh%nfluids == 2) then
                            if (movie%p2) then
                                call writePlotBinaries(q(:,9), "p2", myrank, it)
                            endif
                            if (movie%v2) then
                                call writePlotBinaries(v2plot, "normv2", myrank, it)
                                call writePlotBinaries(q(:,10), "v2x", myrank, it)
                                call writePlotBinaries(q(:,11), "v2z", myrank, it)
                            endif
                        end if
                    end if
                endif

                if (myrank == 0) then
                    write(filename,"('out/energytemp')")
                    call plotPoints2D(stf(:),all_energy,filename)
                    call outputstamp(par, localtime, it, timestamp, maxu, maxv)
                end if
            end if
        end do! timeloop

        ! save seismos
        if (mesh%has_rec) then
            do r=1,rec%nrec
                write(filename,"('out/seismo.x.',i7.7,'.sdu')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotux(r,:),filename)
                write(filename,"('out/seismo.z.',i7.7,'.sdu')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotuz(r,:),filename)
                write(filename,"('out/seismo.x.',i7.7,'.sdv')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotv(r,:),filename)
                write(filename,"('out/seismo.z.',i7.7,'.sdv')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotw(r,:),filename)
                write(filename,"('out/seismo.x.',i7.7,'.sda')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotax(r,:),filename)
                write(filename,"('out/seismo.z.',i7.7,'.sda')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotaz(r,:),filename)
                if (par%poroelastic .and. mesh%nfluids >= 1) then
                    write(filename,"('out/seismo.p.',i7.7,'.sdp1')") rec%recnr(r)
                    call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotp1(r,:),filename)
                    write(filename,"('out/seismo.x.',i7.7,'.sdv1')") rec%recnr(r)
                    call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotv1x(r,:),filename)
                    write(filename,"('out/seismo.z.',i7.7,'.sdv1')") rec%recnr(r)
                    call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotv1z(r,:),filename)
                    if (mesh%nfluids == 2) then
                        write(filename,"('out/seismo.p.',i7.7,'.sdp2')") rec%recnr(r)
                        call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotp2(r,:),filename)
                        write(filename,"('out/seismo.x.',i7.7,'.sdv2')") rec%recnr(r)
                        call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotv2x(r,:),filename)
                        write(filename,"('out/seismo.z.',i7.7,'.sdv2')") rec%recnr(r)
                        call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotv2z(r,:),filename)
                    endif
                endif
                if (par%div) then
                    write(filename,"('out/seismo.r.',i7.7,'.sdv')") rec%recnr(r)
                    call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plot_r(r,:),filename)
                end if
                if (par%curl) then
                    write(filename,"('out/seismo.t.',i7.7,'.sdv')") rec%recnr(r)
                    call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plot_t(r,:),filename)
                end if
            end do
        end if

        if (myrank==0) then
            write(filename,"('out/energy')")
            call plotPoints2D(stf(:)-par%plott0,all_energy,filename)
        end if

        if (par%log.and.myrank==0) then
            write(*,"(a80)") "|                              timeloop finished                               |"
            write(*,"(a80)") "|------------------------------------------------------------------------------|"
        end if
        close(27)

        deallocate(ux,uz,ax,az,uplot,vplot,v1plot,v2plot)
        deallocate(rQ,Q,Qn,Qm,e,resQ, resU, u)
        deallocate(energy, all_energy, energy_kin, energy_pot, a_sigma_tmp)

        deallocate(a_rQ,theta,thetan,thetam, resTheta, int_theta)
        deallocate(stf)
        deallocate(elasticfluxvar, APA, T, invT, VT, VTfree)
        if (par%attenuation) then
            deallocate(anelasticvar, aVT, aVTfree, aT, aAPA)
        endif

        deallocate(srcelemV)
        deallocate(recelemV)

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
            if (par%poroelastic .and. mesh%nfluids >= 1) then
                deallocate(plotp1)
                deallocate(plotv1x)
                deallocate(plotv1z)
                if (mesh%nfluids == 2) then
                    deallocate(plotp2)
                    deallocate(plotv2x)
                    deallocate(plotv2z)
                endif
            endif
        end if

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
            deallocate(np_zeros)
            deallocate(resFprime, resGprime)
        end if
        ! stop MPI
        call finalize_mpi()
    end subroutine timeloop2d
end module timeloopMod
