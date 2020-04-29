!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2020 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2020 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2020 Marc S. Boxberg (RWTH Aachen University, GER)
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
program solver
    ! Program to solve the 2D (poro)elastic wave equation with the discontinous galerkin methode in 2D on triangular mesh
    use parameterMod
    use meshMod
    use timeloopMod
    use sourceReceiverMod
    use errorMessage
    use slipInterfaceMod
    use adjointMod ! TBU
    use collectMovieMod ! TBU

    implicit none

    type(error_message) :: errmsg
    type (parameterVar) :: par
    type(lsi_parameter) :: lsipar
    type (meshVar) :: mesh_change
    type (srcVar) :: src
    type(Invvar) :: inv
    type(movie_parameter) :: movie
    integer :: world_size
    integer :: myrank, i, i_src, BFGS_m, BFGS_i, BFGS_j, l
    integer :: t1, t2, count_max, rate
    integer :: flag, nanalert, nanalert_save
    integer :: n_src_tmp, n_src
    real(kind=CUSTOM_REAL) :: freq_save, misfit, bestalpha, misfittemp, maxstep, minstep, chi_model
    real(kind=CUSTOM_REAL) :: deltat
    real(kind=CUSTOM_REAL) :: dt,tges,tsave
    real(kind=CUSTOM_REAL), dimension(3) :: misfit_total
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: BFGS_p_vp, BFGS_p_vs,BFGS_y_vp, BFGS_y_vs
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: cvs,cvp
    logical :: file_exists
    logical, dimension(:), allocatable :: act_src
    character(len=80) :: filename
    character(len=6) :: myname = "solver"
    character(len=7) :: inv_key='forward'

    !Create new errormessagetrace
    call new(errmsg,myname)

    ! start MPI
    call init_mpi()
    call comm_size(world_size)
    call comm_rank(myrank)

    if (myrank==0) call writeLogo()

    ! read Parfile
    call readParfile(par, lsipar, movie, myrank, errmsg)

    write(filename,"('/meshVar',i6.6)") myrank+1    ! read in mesh in separate variable to change mesh later
    filename=trim(outpath)//trim(filename)
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
        call readMeshVar(mesh_change,filename, errmsg)
    else
        write(*,*) "error in databases, files not existing"
        call add(errmsg, 2, "Error in databases, file "//trim(filename)//" does not exist!", myname, filename)
        call print(errmsg)
        call stop_mpi()
    end if

    if (.level.errmsg == 2) then
        call print(errmsg)
        stop
    endif

        ! get total number of sources
    n_src_tmp=0
    if (mesh_change%has_src) then
        write(filename,"('srcVar',i6.6)") myrank+1
        inquire(file=trim(outpath)//trim(filename), exist=file_exists)
        if (file_exists) then
            call readSrcVar(src,trim(outpath)//filename)
        else
            call add(errmsg, 2, "File does not exist!", myname, filename)
            call print(errmsg)
            call stop_mpi()
        end if
        n_src_tmp=src%nsrc
        call deallocSrcVar(src)
    endif
    call sum_int_all(n_src_tmp, n_src)

    if (par%nproc > 1) call sync_mpi()

    if (par%log .and. myrank == 0) then
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a80)') "|                             start NEXD2D mpi solver                          |"
    end if

    t1 = 0
    t2 = 0
    if (par%log.and.myrank==0) call system_clock(t1, rate, count_max)
    if (par%log.and.myrank==0) tges=0.
    allocate(act_src(n_src))
    act_src=.true.

    call sync_mpi()

    if (.not. par%inversion) then
        ! If no inversion just perform a single forward simulation
        call timeloop2d(par, lsipar, inv_key, par%lowfreq, zero, 0, 1, act_src, movie, myrank, errmsg)
    else
        if (myrank==0) then         ! prepare output streams for inversion
            filename="misfit"
            filename=trim(invpath)//trim(filename)
            open(unit=50,file=trim(filename),status='unknown')
            filename="alpha"
            filename=trim(invpath)//trim(filename)
            open(unit=51,file=trim(filename),status='unknown')
            filename="frequency"
            filename=trim(invpath)//trim(filename)
            open(unit=52,file=trim(filename),status='unknown')
            filename="time"
            filename=trim(invpath)//trim(filename)
            open(unit=54,file=trim(filename),status='unknown')
        endif

        BFGS_m=10               ! m of BFGS calculation defining the number of former iteration steps used for calculation of search direction
        BFGS_j=0                ! iteration variable for BFGS
        freq_save=0.            ! variable to save current used cutoff frequency

        do i=1, par%inversionSteps  ! loop over inversion steps
            call deallocMeshvar(mesh_change)                ! might be helpful to write a function checking if mesh variable is allocated in the first iteration
            write(filename,"('/meshVar',i6.6)") myrank+1    ! read in mesh in separate variable to change mesh later
            filename=trim(outpath)//trim(filename)
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readMeshVar(mesh_change,filename, errmsg)
            else
                write(*,*) "error in databases, files not existing"
                call add(errmsg, 2, "Error in databases, file "//trim(filename)//" does not exist!", myname, filename)
                call print(errmsg)
                call stop_mpi()
            end if

            nanalert=0          ! variable to check later for NaN in results
            tsave=tges          ! to measure computational time
            misfit_total=0.     ! to sum up misfit from different sources
            call readInvFile(inv,myrank,par, errmsg)    ! read file with information on inversion
            if (abs(freq_save-inv%upperfreq) < epsilon(freq_save)) then      ! if frequency is changed compared to former iteration reset BFGS
                BFGS_j=0
                freq_save=inv%upperfreq
            endif

            if (inv%inv_type == 3) then     ! Set BFGS iteration variable
                BFGS_j=BFGS_j+1
            else
                BFGS_j= 1
            endif

            if (inv%inv_type == 3) then     ! override step size for first iteration of BFGS
                if (BFGS_j==1) then
                    if (i==1) then
                        inv%stepsize(1)=0.
                        inv%stepsize(2)=0.03
                        inv%stepsize(3)=0.1
                    else
                        inv%stepsize(1)=0.
                        inv%stepsize(2)=0.003
                        inv%stepsize(3)=0.01
                    endif
                endif
            endif

            if (i==1) then  ! get time step locally in this function to calculate total simulation time
                if (par%autodt) then
                    deltat = mesh_change%dtfactor*par%cfl
                else
                    deltat = par%dt
                endif
                if (.not.par%autont) then
                    par%t_total = (par%nt - 1) * deltat
                    par%autont = .true.
                endif
            endif

            if (BFGS_j==1) then                 ! allocate BFGS vectors and search directions cvs and cvp
                if (.not. allocated(cvs)) then
                    allocate(cvs(mesh_change%nelem),cvp(mesh_change%nelem))
                    allocate(BFGS_p_vp(mesh_change%nelem,BFGS_m),BFGS_p_vs(mesh_change%nelem,BFGS_m))
                    allocate(BFGS_y_vp(mesh_change%nelem,BFGS_m),BFGS_y_vs(mesh_change%nelem,BFGS_m))
                endif
                BFGS_p_vp=0.
                BFGS_p_vs=0.
                BFGS_y_vp=0.
                BFGS_y_vs=0.
            endif
            if (BFGS_j>BFGS_m) then ! counter of BFGS iteration
                BFGS_i=BFGS_m
            else
                BFGS_i=BFGS_j
            endif

            if (par%log.and.myrank==0.and.i==1) then
                write(*,"(a80)") "|------------------------------------------------------------------------------|"
                write(*,"(a80)") "| Start first forward simulation                                               |"
            endif


            do i_src=1,n_src
                !                    if (myrank==0) call system('rm -f '//trim(temppath)//'wavefield*')       ! delete temp wavefield files
                call sync_mpi()
                act_src=.false.     ! reset active sources
                act_src(i_src)=.true.
                inv_key='forward'   ! define kind of forward simulation
                call timeloop2d(par, lsipar, inv_key, par%lowfreq, inv%upperfreq, i, 1, act_src, movie, myrank, errmsg) ! call first forward simulation

                call sync_mpi()

                if (par%log.and.myrank==0) then
                    call  system_clock(t2, rate, count_max) ! call timer to measure run time of the program
                    dt=(t2 - t1)/rate
                    tges=tges+dt
                    write(*,"(a80)") "|------------------------------------------------------------------------------|"
                    write(*,"(a41,i3,a15,i3,a18)") "| Finished forward simulation for source ", i_src, " of iteration: ", i, "                 |"
                    write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                    write(*,"(a80)") "|                                                                              |"
                    write(*,"(a80)") "| Start calculating misfit and make adjoint sources                            |"
                    t1=t2
                endif

                call getMisfit(misfit,mesh_change,myrank,inv%upperfreq,i_src)   ! get misfit for the seismograms of the current model for this rank
                call sync_mpi()
                call sum_real_all(misfit,misfittemp,CUSTOM_REAL)                    ! sum up the misfit for all ranks

                misfit_total(1)=misfit_total(1) + misfittemp                        ! sum up misfit for all sources

                if (i_src==1) write(52,*) i,inv%upperfreq                           ! write cut off frequency to file
                if (par%log.and.myrank==0) then
                    write(*,"(a10,es10.3,a12,i3,a4,i3,a38)") '| Misfit= ', misfit_total(1), ' for source ', i_src, ' of ', n_src, '                                     |'
                endif

                call makeAdjSourceElastic(par,myrank,inv,i,i_src)                   ! make adjoint source time functions
                call sync_mpi()

                if (par%log.and.myrank==0) then
                    call  system_clock(t2, rate, count_max)
                    dt=(t2 - t1)/rate
                    tges=tges+dt
                    write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                    write(*,"(a80)") "|                                                                              |"
                    write(*,"(a80)") "| Start adjoint forward simulation                                             |"
                    t1=t2
                endif

                inv_key='adjoint'   ! change flag for kind of simulation to adjoint
                call timeloop2d(par, lsipar, inv_key, par%lowfreq, inv%upperfreq, i, 1, act_src, movie, myrank, errmsg) ! cal adjoint simulation
                call sync_mpi()

                if (par%log.and.myrank==0) then
                    call  system_clock(t2, rate, count_max)
                    dt=(t2 - t1)/rate
                    tges=tges+dt
                    write(*,"(a80)") "|------------------------------------------------------------------------------|"
                    write(*,"(a41,i3,a14,i3,a19)") "| Finished adjoint simulation for source ", i_src, " of iteration ", i, "                  |"
                    write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                    if (i_src==n_src) then
                        write(*,"(a80)") "|                                                                              |"
                        write(*,"(a80)") "| Start calculation of the search direction and change model                   |"
                    else
                        write(*,"(a80)") "|                                                                              |"
                        write(*,"(a80)") "| Start forward simulation for next source                                     |"
                    endif
                    t1=t2
                endif
            enddo   ! loop over sources
            call collectKernels(myrank,i,mesh_change,n_src)  ! collect kernels from all sources
            !call collectKernels(myrank,i,mesh_change,n_src, inv,chi_model)  ! collect kernels from all sources
            if (myrank==0) write(50,*) i-1, misfit_total(1), chi_model

            if (inv%inv_type == 1) call getSteepestDescent(cvs,cvp,mesh_change,i,myrank)    ! depending on the inversion type get the corresponding search direction
            if (inv%inv_type == 2) call getConjugateGradient(cvs,cvp,mesh_change,i,myrank)
            if (inv%inv_type == 3) then
                if (BFGS_j>1) call getBFGS_y(BFGS_y_vp, BFGS_y_vs,mesh_change,i,BFGS_j,myrank,BFGS_m)   ! get new y vector for BFGS
                call getBFGS(cvs,cvp,BFGS_p_vp, BFGS_p_vs,BFGS_y_vp, BFGS_y_vs, mesh_change,i,BFGS_j,myrank,BFGS_m)
            endif

            !            if (inv%smooth_grad) call smoothArray(cvs, cvp, myrank, mesh_change, inv%depth_grad, inv%weight_grad)

            if (par%nproc > 1) call sync_mpi()
            ! change model with search direction and test step width
            call changeModel(inv%stepsize(2),mesh_change,par,BFGS_p_vp(:,BFGS_i), BFGS_p_vs(:,BFGS_i), cvs, cvp, myrank, inv,  i, BFGS_j,0, errmsg)

            if (par%log.and.myrank==0) then
                call  system_clock(t2, rate, count_max)
                dt=(t2 - t1)/rate
                tges=tges+dt
                write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                write(*,"(a80)") "|                                                                              |"
                write(*,"(a80)") "| Start first test simulation                                                  |"
                t1=t2
            endif

            do i_src=1,n_src
                act_src=.false.         ! reset active sources
                act_src(i_src)=.true.
                inv_key='forward'
                call timeloop2d(par, lsipar, inv_key, par%lowfreq, inv%upperfreq, i, 2, act_src, movie, myrank, errmsg) ! first test forward simulation

                call getMisfit(misfit,mesh_change,myrank,inv%upperfreq,i_src)       ! get misfit for first test simulation

                if (par%nproc > 1) call sync_mpi()                      ! sum up misfit over ranks and sources
                if (par%nproc > 1) then
                    call sum_real_all(misfit,misfittemp,CUSTOM_REAL)
                else
                    misfittemp=misfit
                endif
                misfit_total(2) = misfit_total(2) + misfittemp

                if (par%log.and.myrank==0) then
                    call  system_clock(t2, rate, count_max)
                    dt=(t2 - t1)/rate
                    tges=tges+dt
                    write(*,"(a80)") "|------------------------------------------------------------------------------|"
                    write(*,"(a52,i3,a25)") "| Finished first test forward simulation for source ", i_src,"                        |"
                    write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                    if (i_src==n_src) then
                        write(*,"(a80)") "|                                                                              |"
                        write(*,"(a80)") "| Change model and start second test forward simulation                        |"
                    else
                        write(*,"(a80)") "|                                                                              |"
                        write(*,"(a80)") "| Start first test simulation for next source                                  |"
                    endif
                    t1=t2
                endif
            enddo
            ! change model for second test simulation
            call changeModel(inv%stepsize(3),mesh_change,par,BFGS_p_vp(:,BFGS_i), BFGS_p_vs(:,BFGS_i), cvs, cvp, myrank, inv,  i, BFGS_j,1, errmsg)

            do i_src=1,n_src    ! do the same as for the first test simulation
                act_src=.false.
                act_src(i_src)=.true.
                call timeloop2d(par, lsipar, inv_key, par%lowfreq, inv%upperfreq, i, 3, act_src, movie, myrank, errmsg)

                call getMisfit(misfit,mesh_change,myrank,inv%upperfreq,i_src)

                if (par%nproc > 1) call sync_mpi()
                if (par%nproc > 1) then
                    call sum_real_all(misfit,misfittemp,CUSTOM_REAL)
                else
                    misfittemp=misfit
                endif
                misfit_total(3)=misfit_total(3) + misfittemp

                if (par%log.and.myrank==0) then
                    call system_clock(t2, rate, count_max)
                    dt=(t2 - t1)/rate
                    tges=tges+dt
                    write(*,"(a80)") "|------------------------------------------------------------------------------|"
                    write(*,"(a53,i3,a24)") "| Finished second test forward simulation for source ", i_src,"                       |"
                    write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                    if (i_src==n_src) then
                        write(*,"(a80)") "|                                                                              |"
                        write(*,"(a80)") "| Start search for best step size                                              |"
                    else
                        write(*,"(a80)") "|                                                                              |"
                        write(*,"(a80)") "| Start second test simulation for next source                                 |"
                    endif
                    t1=t2
                endif
            enddo

            if (par%log.and.myrank==0) then
                write(*,"(a16,f10.7,a2,f10.7,a2,f10.7,a30)") '| step sizes:   ', inv%stepsize(1), ', ', inv%stepsize(2), ', ', inv%stepsize(3),'                             |'
                write(*,"(a16,es10.3,a2,es10.3,a2,es10.3,a30)") '| misfit values:', misfit_total(1), ', ', misfit_total(2), ', ', misfit_total(3),'                             |'
                if (BFGS_j==1) then
                    write(*,"(a80)") "| Step sizes are percentages of the maximum velocity change                    |"
                else
                    write(*,"(a80)") "| Step sizes are normalized BFGS steps                                         |"
                endif
            endif
            call findBestAlpha(inv%stepsize,misfit_total,bestalpha,nanalert,flag)    ! find the best step length
            nanalert_save=nanalert                                              ! save if one misfit calculation produced NaN

            if (BFGS_j==1) then                 ! get maximum and minimum values for the test step lengths to prevent infinite loops
                maxstep=inv%max_step
                minstep=inv%min_step
            else
                maxstep=inv%max_step_BFGS
                minstep=inv%min_step_BFGS
            endif

            do while((flag==-1 .or. flag==-2) .and. (inv%stepsize(3)<=maxstep .and. inv%stepsize(3)>=minstep))
                if (flag==-1) then     ! halve test step size if no optimal step is found but can be expected at smaller values
                    inv%stepsize(1)=inv%stepsize(1)
                    inv%stepsize(3)=inv%stepsize(2)
                    inv%stepsize(2)=0.5*inv%stepsize(2)
                    misfit_total(1)=misfit_total(1)
                    misfit_total(3)=misfit_total(2)
                    misfit_total(2)=0.
                    if (par%log.and.myrank==0) then
                        write(*,"(a80)") "|------------------------------------------------------------------------------|"
                        write(*,"(a80)") "| Used step lengths are too large. Try smaller ones.                           |"
                        write(*,"(a24,f10.5,a2,f10.5,a2,f10.5,a22)") "| New step lengths are: ", inv%stepsize(1), ', ', inv%stepsize(2), ', ', inv%stepsize(3),'                     |'
                        write(*,"(a80)") "| Change model for new step length                                             |"
                        t1=t2
                    endif
                    call changeModel(inv%stepsize(2),mesh_change,par,BFGS_p_vp(:,BFGS_i), BFGS_p_vs(:,BFGS_i), cvs, cvp, myrank, inv,  i, BFGS_j,1, errmsg)
                else        ! double the test step lengths if the minimum of the quadratic function is at step lengths larger than the current values
                    inv%stepsize(1)=inv%stepsize(2)
                    inv%stepsize(2)=inv%stepsize(3)
                    inv%stepsize(3)=2*inv%stepsize(3)
                    misfit_total(1)=misfit_total(2)
                    misfit_total(2)=misfit_total(3)
                    misfit_total(3)=0.
                    if (par%log.and.myrank==0) then
                        write(*,"(a80)") "|------------------------------------------------------------------------------|"
                        write(*,"(a80)") "| Used step lengths are too small. Try larger ones.                            |"
                        write(*,"(a24,f10.5,a2,f10.5,a2,f10.5,a22)") "| New step lengths are: ", inv%stepsize(1), ', ', inv%stepsize(2), ', ', inv%stepsize(3),'                     |'
                        write(*,"(a80)") "| Change model for new step length                                             |"
                        t1=t2
                    endif
                    call changeModel(inv%stepsize(3),mesh_change,par,BFGS_p_vp(:,BFGS_i), BFGS_p_vs(:,BFGS_i), cvs, cvp, myrank, inv,  i, BFGS_j,1, errmsg)
                endif

                do i_src=1,n_src    ! run an additional simulation for all sources
                    act_src=.false.
                    act_src(i_src)=.true.
                    if (flag == -1) then ! distinguish between new first or second test simulation
                        call timeloop2d(par, lsipar, inv_key, par%lowfreq, inv%upperfreq, i, 2, act_src, movie, myrank, errmsg)
                    else
                        call timeloop2d(par, lsipar, inv_key, par%lowfreq, inv%upperfreq, i, 3, act_src, movie, myrank, errmsg)
                    endif
                    call getMisfit(misfit,mesh_change,myrank,inv%upperfreq, i_src)

                    if (par%nproc > 1) call sync_mpi()
                    if (par%nproc > 1) then
                        call sum_real_all(misfit,misfittemp,CUSTOM_REAL)
                    else
                        misfittemp=misfit
                    endif
                    if (flag == -1) then ! change the misfit value depending on the kind of the new test simulation
                        misfit_total(2)=misfit_total(2) + misfittemp
                    else
                        misfit_total(3)=misfit_total(3) + misfittemp
                    endif
                    if (par%log.and.myrank==0) then
                        call  system_clock(t2, rate, count_max)
                        dt=(t2 - t1)/rate
                        tges=tges+dt
                        write(*,"(a80)") "|------------------------------------------------------------------------------|"
                        write(*,"(a51,f10.5,a11,i3,a5)") "| Finished simulation with new largest step length ",inv%stepsize(3)," for source", i_src, '    |'
                        write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                        if (i_src==n_src) then
                            write(*,"(a80)") "|                                                                              |"
                            write(*,"(a80)") "| Find best step size                                                          |"
                        else
                            write(*,"(a80)") "|                                                                              |"
                            write(*,"(a80)") "| Start second test simulation for next source                                 |"
                        endif
                        write(*,"(a80)") "|------------------------------------------------------------------------------|"
                        t1=t2
                    endif
                enddo
                if (par%log.and.myrank==0) then
                    write(*,"(a16,f10.7,a2,f10.7,a2,f10.7,a30)") '| step sizes:   ', inv%stepsize(1), ', ', inv%stepsize(2), ', ', inv%stepsize(3),'                             |'
                    write(*,"(a16,es10.3,a2,es10.3,a2,es10.3,a30)") '| misfit values:', misfit_total(1), ', ', misfit_total(2), ', ', misfit_total(3),'                             |'
                    if (BFGS_j==1) then
                        write(*,"(a80)") "| Step sizes are percentages of the maximum velocity change                    |"
                    else
                        write(*,"(a80)") "| Step sizes are normalized BFGS steps                                         |"
                    endif
                endif
                call findBestAlpha(inv%stepsize,misfit_total,bestalpha,nanalert,flag)    ! find again an optimal step length

                if (nanalert_save==1 .and. nanalert==0) then
                    write(*,"(a80)") "|------------------------------------------------------------------------------|"
                    write(*,"(a80)") "| EXIT unstable simulation!                                                    |"
                    write(*,"(a80)") "|------------------------------------------------------------------------------|"
                    exit
                endif
                nanalert_save=nanalert
            enddo
            if (flag == -2 .or. flag == -1) then
                bestalpha=inv%stepsize(3)
                bestalpha=inv%stepsize(3)
                write(*,"(a80)") "|------------------------------------------------------------------------------|"
                write(*,"(a7,i3,a2)") "| Cannot find proper step length. Just do step with last tried step length ", inv%stepsize(3)," |"
                write(*,"(a80)") "|------------------------------------------------------------------------------|"
            endif

            if (myrank==0) write(51,*) i, bestalpha
            if (par%log .and. myrank==0) then
                call  system_clock(t2, rate, count_max)
                dt=(t2 - t1)/rate
                tges=tges+dt
                write(*,"(a80)") "|------------------------------------------------------------------------------|"
                write(*,"(a24,f10.5,a46)") '| Found best step size: ', bestalpha, '                                             |'
                write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                write(*,"(a80)") "|                                                                              |"
                write(*,"(a34, i3,a43 )") "| Change model for inversion step ", i, ' and create movie files                   |'
                t1=t2
            endif

            ! change model with the final optimal step length
            call changeModel(bestalpha,mesh_change,par,BFGS_p_vp(:,BFGS_i), BFGS_p_vs(:,BFGS_i), cvs, cvp, myrank, inv,  i, BFGS_j,3, errmsg)
            call deallocMeshvar(mesh_change)    ! delete current mesh
            call sync_mpi()
            if (myrank==0) call collectMovie(par, movie,  i, errmsg)    ! make vtk files for the current iteration step
            call sync_mpi()
            if (par%log .and. myrank==0) then
                call  system_clock(t2, rate, count_max)
                dt=(t2 - t1)/rate
                tges=tges+dt
                write(*,"(a80)") "| Finished movie files                                                         |"
                write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                write(*,"(a80)") "|                                                                              |"
                write(*,"(a80)") "| Iteration step finished                                                      |"
                write(*,"(a80)") "|------------------------------------------------------------------------------|"
                t1=t2
            endif

            if (myrank==0) write(54,*) i, tges-tsave, tges
            if (BFGS_j>BFGS_m) then                 ! reset the entries of the BFGS vectors
                do l=1,BFGS_m-1
                    BFGS_p_vs(:,l)=BFGS_p_vs(:,l+1)
                    BFGS_p_vp(:,l)=BFGS_p_vp(:,l+1)
                    BFGS_y_vs(:,l)=BFGS_y_vs(:,l+1)
                    BFGS_y_vp(:,l)=BFGS_y_vp(:,l+1)
                enddo
                BFGS_p_vs(:,BFGS_m)=0.
                BFGS_p_vp(:,BFGS_m)=0.
                BFGS_y_vs(:,BFGS_m)=0.
                BFGS_y_vp(:,BFGS_m)=0.
            endif
        enddo

        if (par%log .and. myrank==0) then
            write(*,"(a80)") "| Perform an additional forward simulation to obtain the final misfit          |"
            t1=t2
        endif

        if (par%nproc > 1) call sync_mpi()
        write(filename,"('/meshVar',i6.6)") myrank+1    ! read mesh again for final simulation
        filename=trim(outpath)//trim(filename)
        inquire(file=trim(filename), exist=file_exists)
        if (file_exists) then
            call readMeshVar(mesh_change,filename, errmsg)
        else
            write(*,*) "error in databases, files not existing"
            call add(errmsg, 2, "Error in databases, file "//trim(filename)//" does not exist!", myname, filename)
            call print(errmsg)
            call stop_mpi()
        end if
        misfit_total(1)=0.

        do i_src=1,n_src    ! do a final forward simulation to obtain misfit for final model
            act_src=.false.
            act_src(i_src)=.true.
            inv_key='forward'
            call timeloop2d(par, lsipar, inv_key, par%lowfreq, inv%upperfreq, par%inversionSteps+1, 1, act_src, movie, myrank, errmsg)

            call getMisfit(misfit,mesh_change,myrank,inv%upperfreq,i_src)

            if (par%nproc > 1) call sync_mpi()
            if (par%nproc > 1) then
                call sum_real_all(misfit,misfittemp,CUSTOM_REAL)
            else
                misfittemp=misfit
            endif
            misfit_total(1) = misfit_total(1) + misfittemp

            if (par%log.and.myrank==0) then
                call  system_clock(t2, rate, count_max)
                dt=(t2 - t1)/rate
                tges=tges+dt
                write(*,"(a80)")                 "|------------------------------------------------------------------------------|"
                write(*,"(a37,i3,a40)") "| Finished last simulation for source", i_src, '                                       |'
                write(*,"(a19,f7.1,a34,f10.1,a10)") "| Calculation time ", dt , " seconds, total computation time: ", tges/3600., " hours   |"
                if (i_src==n_src) then
                    write(*,"(a80)") "| Finalize inversion by deleting old wave field files                          |"
                else
                    write(*,"(a80)") "|                                                                              |"
                    write(*,"(a80)") "| Start last simulation for next source                                        |"
                endif
                t1=t2
            endif
        enddo
        if (myrank==0) then
            write(50,*) par%inversionSteps, misfit_total(1), chi_model
            close(50)
            close(51)
            close(52)
        endif
        if (par%nproc > 1) call sync_mpi()
        call deallocMeshvar(mesh_change)
        deallocate(cvs, cvp,BFGS_p_vp, BFGS_p_vs,BFGS_y_vp,BFGS_y_vs)

        if (myrank==0) then
            call system('rm -f '//trim(temppath)//'/wavefield*')       ! delete temp wavefield files
        endif
        if (par%log .and. myrank == 0)  write(*,'(a80)') "|------------------------------------------------------------------------------|"
    endif
    deallocate(act_src)

    if (par%log .and. myrank == 0) then
        write(*,'(a80)') "|                            end NEXD2D mpi solver                             |"
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
    end if

    ! stop MPI
    call finalize_mpi()
end program solver
