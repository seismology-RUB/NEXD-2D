!-----------------------------------------------------------------------
!   Copyright 2014-2020 Andre Lamert (Ruhr-Universität Bochum, GER)
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
module adjointMod
    !use constantsMod
    use meshMod
    use parameterMod
    use mpiMod
    use timeSeries
    use dtMod
    use collectMovieMod
    implicit none

    interface makeTimeSeries
        module procedure makeTimeSeriesDP
        module procedure makeTimeSeriesSP
    end interface
    interface cubicInterpolation
        module procedure cubicInterpolationArray
        module procedure cubicInterpolationSnglVal
    end interface


    type :: invVar
        real(kind=CUSTOM_REAL) :: upperfreq                                         ! cutoff frequency for low-pass filtering
        integer :: inv_type                                                         ! flag for the search type of the search direction (1: Steepest descent, 2: Conjugate gradient, 3: BFGS)
        real(kind=CUSTOM_REAL), dimension(3) :: stepsize                            ! Length of test step sizes to be used
        real(kind=CUSTOM_REAL) :: min_step, max_step, min_step_BFGS, max_step_BFGS  ! minimum and maximum values for test step sizes during search for best step size
        real(kind=CUSTOM_REAL) :: min_vp, max_vp, min_vs, max_vs                    ! Minimum and maximum valocity values for the reconstructed models
        !logical :: smooth_grad                                                      ! Smooth the gradient by averaging over neighbor elements (EXPERIMENTAL!)
        !real(kind=CUSTOM_REAL) :: weight_grad                                       ! Weighting of smoothing (EXPERIMENTAL!)
        !integer :: depth_grad                                                       ! Number of iterations used for smoothing of gradient (EXPERIMENTAL!)
        !logical :: smooth_model                                                     ! Smooth model (similar to gradient) (EXPERIMENTAL!)
        !real(kind=CUSTOM_REAL) :: weight_model                                      ! Weighting for model smoothing (EXPERIMENTAL!)
        !integer :: depth_model                                                      ! Number of iterations for model smoothing (EXPERIMENTAL!)
        !logical :: use_smoothing                                                    ! Use Smoothing at all
        !real(kind=CUSTOM_REAL) :: smooth_alpha                                      ! Weighting between model misfit and misfit from measurements

    end type invVar
contains

    subroutine readInvfile(this,myrank,par, errmsg)
    ! Read the file "data/invpar" with information regarding the inversion
        implicit none
        type(error_message) :: errmsg
        type (InvVar) :: this
        type(parameterVar):: par
        integer ::myrank

        !local variables
        character(len=80) :: filename
        integer :: ier

        filename=trim('data/invpar')
        open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
        if (myrank==0) then
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a26, a14, a12, a28)") "|                          ","Begin reading ", filename, "...                        |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        endif
        ! Cycle to read the parameters form the parameter file. The order of appearence of the parameters is not important

        call readFloatPar(this%upperfreq, "upperfreq", filename, 0, errmsg)
        call readIntPar(this%inv_type, "inv_type", filename, 0, errmsg)
        call readFloatPar(this%stepsize(1), "step_1", filename, 0, errmsg)
        call readFloatPar(this%stepsize(2), "step_2", filename, 0, errmsg)
        call readFloatPar(this%stepsize(3), "step_3", filename, 0, errmsg)
        call readFloatPar(this%min_step, "min_step", filename, 0, errmsg)
        call readFloatPar(this%max_step, "max_step", filename, 0, errmsg)
        call readFloatPar(this%min_step_BFGS, "min_step_BFGS", filename, 0, errmsg)
        call readFloatPar(this%max_step_BFGS, "max_step_BFGS", filename, 0, errmsg)
        call readFloatPar(this%min_vp, "min_vp", filename, 0, errmsg)
        call readFloatPar(this%max_vp, "max_vp", filename, 0, errmsg)
        call readFloatPar(this%min_vs, "min_vs", filename, 0, errmsg)
        call readFloatPar(this%max_vs, "max_vs", filename, 0, errmsg)
!        call readLogicalPar(this%smooth_grad, "smooth_grad", filename, 0, errmsg)
!        call readFloatPar(this%weight_grad, "weight_grad", filename, 0, errmsg)
!        call readIntPar(this%depth_grad, "depth_grad", filename, 0, errmsg)
!        call readLogicalPar(this%smooth_model, "smooth_model", filename, 0, errmsg)
!        call readFloatPar(this%weight_model, "weight_model", filename, 0, errmsg)
!        call readIntPar(this%depth_model, "depth_model", filename, 0, errmsg)
!        call readLogicalPar(this%use_smoothing, "use_smoothing", filename, 0, errmsg)
!        call readFloatPar(this%smooth_alpha, "smooth_alpha", filename, 0, errmsg)


        ! print information read from file to screen
        if (par%log.and.myrank==0) then
            write (*,"(a40, f10.1, a30)") "|                    Highest frequency: ", this%upperfreq, "                             |"
            write (*,"(a40, i10, a30)")   "|                     Inversion method: ", this%inv_type,  "                             |"
            write (*,"(a40, f10.3, a30)") "|                            Stepsizes: ", this%stepsize(1),"                             |"
            write (*,"(a40, f10.3, a30)") "|                                       ", this%stepsize(2),"                             |"
            write (*,"(a40, f10.3, a30)") "|                                       ", this%stepsize(3),"                             |"
            write (*,"(a40, f10.1, a30)") "|               Minimum test step size: ", this%min_step, "                             |"
            write (*,"(a40, f10.1, a30)") "|               Maximum test step size: ", this%max_step, "                             |"
            write (*,"(a40, f10.1, a30)") "|               Minimum BFGS step size: ", this%min_step_BFGS, "                             |"
            write (*,"(a40, f10.1, a30)") "|               Maximum BFGS step size: ", this%max_step_BFGS, "                             |"
            write (*,"(a40, f10.1, a30)") "|                  Minimal vp velocity: ", this%min_vp,  "                             |"
            write (*,"(a40, f10.1, a30)") "|                  Maximal vp velocity: ", this%max_vp,  "                             |"
            write (*,"(a40, f10.1, a30)") "|                  Minimal vs velocity: ", this%min_vs,  "                             |"
            write (*,"(a40, f10.1, a30)") "|                  Maximal vs velocity: ", this%max_vs,  "                             |"
!            write (*,"(a40, l10, a30)")   "|                      Smooth gradient: ", this%smooth_grad, "                             |"
!            write (*,"(a40, f10.1, a30)") "|           Smoothing weights gradient: ", this%weight_grad, "                             |"
!            write (*,"(a40, i10, a30)")   "|             Smoothing depth gradient: ", this%depth_grad,  "                             |"
!            write (*,"(a40, l10, a30)")   "|                         Smooth model: ", this%smooth_model, "                             |"
!            write (*,"(a40, f10.1, a30)") "|              Smoothing weights model: ", this%weight_model, "                             |"
!            write (*,"(a40, i10, a30)")   "|                Smoothing depth model: ", this%depth_model,  "                             |"
!            write (*,"(a40, l10, a30)")   "|                        Use smoothing: ", this%use_smoothing,  "                             |"
!            write (*,"(a40, f10.1, a30)") "|                      Smoothing alpha: ", this%smooth_alpha,  "                             |"
        endif
        close(19)
    end subroutine readInvfile

    subroutine LowPassZeroPhaseFilter(stf,nt,deltat,upperfreq, tapern, frac)
    ! Lowpass filter a timeseries to a given cutoff frequency
        type(time_Series) :: timeseries, temp                       ! Use time_series type from timeSeriesMod
        real(kind=CUSTOM_REAL), dimension(:) :: stf                 ! source time function to be filtered
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: stf_temp   ! temp array for stf
        real(kind=CUSTOM_REAL) :: deltat, upperfreq                 ! deltat: time sampling, upperfreq: cutoff frequency
        integer :: nt,j, frac                                       ! nt: length of stf, frac: fraction of stf to be tapered on both ends
        logical :: tapern                                           ! flag for tappering

        allocate(stf_temp(size(stf)))

        call makeTimeSeries(timeseries,nt,deltat,stf)               ! Create TimeSeries to apply following functions
        if (tapern) timeseries=hanningTaperTimeSeries(timeseries,real(nt/frac*deltat))  ! Tapering
        timeseries=lowPassButterworthRecursiveTimeSeries(timeseries,real(upperfreq),4)  ! Low-pass filter of degree 4
        temp=ReverseTimeSeries(timeseries)                                              ! Reverse time series to get zero-time-shift filter
        temp=lowPassButterworthRecursiveTimeSeries(temp,real(upperfreq),4)              ! Filter again
        timeseries=ReverseTimeSeries(temp)                                              ! Reverse time series to old order
        do j=1, nt                                                                      ! Extract information from timeSeries type to old stf array
            stf_temp(j)=getSampleTimeSeries(timeseries,j)
        enddo
        stf=stf_temp
        deallocate(stf_temp)
    end subroutine LowPassZeroPhaseFilter

    subroutine cubicInterpolationArray(t_old, y_old, t_new, y_new)
        ! do an interpolation based on four points and a cubic polynomial. Finds the interpolated value y_new at times t_new.
        ! t_old and t_new need to be ordered at equally spaced!
        real(kind=8), dimension(:) :: t_old, t_new              ! time samples
        real(kind=8), dimension(:) :: y_old, y_new              ! values at time samples
        real(kind=8) ::a,b,c,d,y12,t23,t34,t24,y23,t12,t13,y34,t1t1,t1t2,t2t3,t3t3,t2t2,t3t4,t4t4, t1pt2    ! interpolation variables
        integer :: i, ind, length_old, length_new               ! counter
        length_old=size(t_old)                                  ! just get the size of old and new array
        length_new=size(t_new)
        do i=1,length_new                                       ! for all values of the new array
            if (i <= length_old) then                           ! Need this if for next if
                if (abs(t_new(i)-t_old(i)) < 1.E-6) then        ! If both time values are extremly close to each other, just copy the y value
                    y_new(i)=y_old(i)
                else
                    ! Find the position of the t_old value closest but smaller as currrently used t_new value
                    ! Found index is always the second of four points for the cubic interpolation
                    call findindex(real(t_old,kind=CUSTOM_REAL),real(t_new(i),kind=CUSTOM_REAL),ind)
                    if (ind==1) then    ! Reset index if the interpolation is close to the ends of the old array
                        ind=2
                    else if (length_old-ind==1) then
                        ind=length_old-2
                    else if (length_old-ind==0) then
                        ind=length_old-2
                    endif
                    ! Get parameters for interpolation
                    y12=y_old(ind-1)-y_old(ind)
                    t23=t_old(ind)-t_old(ind+1)
                    t34=t_old(ind+1)-t_old(ind+2)
                    t24=t_old(ind)-t_old(ind+2)
                    y23=y_old(ind)-y_old(ind+1)
                    t12=t_old(ind-1)-t_old(ind)
                    t13=t_old(ind-1)-t_old(ind+1)
                    y34=y_old(ind+1)-y_old(ind+2)
                    t1t1=t_old(ind-1)*t_old(ind-1)
                    t1t2=t_old(ind-1)*t_old(ind)
                    t2t3=t_old(ind)*t_old(ind+1)
                    t3t3=t_old(ind+1)*t_old(ind+1)
                    t2t2=t_old(ind)*t_old(ind)
                    t3t4=t_old(ind+1)*t_old(ind+2)
                    t4t4=t_old(ind+2)*t_old(ind+2)
                    t1pt2=t_old(ind-1)+t_old(ind)
                    ! get parameters of cubic polynomial
                    a=(y12*t23*t34*t24-y23*t12*t34*(t24+t13)+y34*t23*t13*t12)/(t12*t23*t34)/&
                        ((t1t1+t1t2-t2t3-t3t3)*t24-(t2t2+t2t3-t3t4-t4t4)*t13)
                    b=(y23*t34-y34*t23)/(t23*t34*t24)-a*(t2t2+t2t3-t3t4-t4t4)/t24
                    c=y12/t12-a*(t1t1+t1t2+t2t2)-b*t1pt2
                    d=y_old(ind)-a*t_old(ind)**3-b*t_old(ind)**2-c*t_old(ind)
                    ! Get new value
                    y_new(i)=a*t_new(i)**3+b*t_new(i)**2+c*t_new(i)+d
                endif
            else
                ! Similar as above
                call findindex(real(t_old,kind=CUSTOM_REAL),real(t_new(i),kind=CUSTOM_REAL),ind)
                if (ind==1) then
                    ind=2
                else if (length_old-ind==1) then
                    ind=length_old-2
                else if (length_old-ind==0) then
                    ind=length_old-2
                endif
                y12=y_old(ind-1)-y_old(ind)
                t23=t_old(ind)-t_old(ind+1)
                t34=t_old(ind+1)-t_old(ind+2)
                t24=t_old(ind)-t_old(ind+2)
                y23=y_old(ind)-y_old(ind+1)
                t12=t_old(ind-1)-t_old(ind)
                t13=t_old(ind-1)-t_old(ind+1)
                y34=y_old(ind+1)-y_old(ind+2)
                t1t1=t_old(ind-1)*t_old(ind-1)
                t1t2=t_old(ind-1)*t_old(ind)
                t2t3=t_old(ind)*t_old(ind+1)
                t3t3=t_old(ind+1)*t_old(ind+1)
                t2t2=t_old(ind)*t_old(ind)
                t3t4=t_old(ind+1)*t_old(ind+2)
                t4t4=t_old(ind+2)*t_old(ind+2)
                t1pt2=t_old(ind-1)+t_old(ind)
                a=(y12*t23*t34*t24-y23*t12*t34*(t24+t13)+y34*t23*t13*t12)/(t12*t23*t34)/&
                    ((t1t1+t1t2-t2t3-t3t3)*t24-(t2t2+t2t3-t3t4-t4t4)*t13)
                b=(y23*t34-y34*t23)/(t23*t34*t24)-a*(t2t2+t2t3-t3t4-t4t4)/t24
                c=y12/t12-a*(t1t1+t1t2+t2t2)-b*t1pt2
                d=y_old(ind)-a*t_old(ind)**3-b*t_old(ind)**2-c*t_old(ind)
                y_new(i)=a*t_new(i)**3+b*t_new(i)**2+c*t_new(i)+d
            endif
        enddo
    end subroutine cubicInterpolationArray

    subroutine cubicInterpolationSnglVal(t_old, y_old, t_new, y_new)
    ! Does the same as cubicInterpolationArray for just one value of t_new
        real(kind=8), dimension(:) :: t_old, y_old
        real(kind=8) :: t_new
        real(kind=CUSTOM_REAL) :: y_new
        real(kind=8) ::a,b,c,d,y12,t23,t34,t24,y23,t12,t13,y34,t1t1,t1t2,t2t3,t3t3,t2t2,t3t4,t4t4, t1pt2
        integer :: ind, length_old
        length_old=size(t_old)

        call findindex(real(t_old,kind=CUSTOM_REAL),real(t_new,kind=CUSTOM_REAL),ind)
        if (ind==1) then
            ind=2
        else if (length_old-ind==1) then
            ind=length_old-2
        else if (length_old-ind==0) then
            ind=length_old-2
        endif
        y12=y_old(ind-1)-y_old(ind)
        t23=t_old(ind)-t_old(ind+1)
        t34=t_old(ind+1)-t_old(ind+2)
        t24=t_old(ind)-t_old(ind+2)
        y23=y_old(ind)-y_old(ind+1)
        t12=t_old(ind-1)-t_old(ind)
        t13=t_old(ind-1)-t_old(ind+1)
        y34=y_old(ind+1)-y_old(ind+2)
        t1t1=t_old(ind-1)*t_old(ind-1)
        t1t2=t_old(ind-1)*t_old(ind)
        t2t3=t_old(ind)*t_old(ind+1)
        t3t3=t_old(ind+1)*t_old(ind+1)
        t2t2=t_old(ind)*t_old(ind)
        t3t4=t_old(ind+1)*t_old(ind+2)
        t4t4=t_old(ind+2)*t_old(ind+2)
        t1pt2=t_old(ind-1)+t_old(ind)
        a=(y12*t23*t34*t24-y23*t12*t34*(t24+t13)+y34*t23*t13*t12)/(t12*t23*t34)/&
            ((t1t1+t1t2-t2t3-t3t3)*t24-(t2t2+t2t3-t3t4-t4t4)*t13)
        b=(y23*t34-y34*t23)/(t23*t34*t24)-a*(t2t2+t2t3-t3t4-t4t4)/t24
        c=y12/t12-a*(t1t1+t1t2+t2t2)-b*t1pt2
        d=y_old(ind)-a*t_old(ind)**3-b*t_old(ind)**2-c*t_old(ind)
        y_new=real(a*t_new**3+b*t_new**2+c*t_new+d, kind=CUSTOM_REAL)
    end subroutine cubicInterpolationSnglVal

    subroutine findindex(t_old,t_new,ind)
    ! finds the index (ind) of a value t_new within the array t_old. The found index corresponds to the closest but smaller value in t_old compared to t_new
        real(kind=CUSTOM_REAL), dimension(:) :: t_old
        real(kind=CUSTOM_REAL) :: t_new
        integer :: ind, left, right, middle, size_t

        ind=-1
        left=1
        size_t=size(t_old)
        right=size_t
        middle=(right+left)/2
        if (t_new >= t_old(right)) then     ! check if t_new is larger than the largest t_old value
            left=right
            ind=right   ! set index to largest value of t_old
            ! if t_new is only one time step larger as the largest t_old value it can still be used
            if (t_new > 1.0001*(2*t_old(right)-t_old(right-1))) then ! factor 1.0001 to handle bugs caused by numeric accuracy of numbers
                write(*,*) 't_new value outside of t_old range', t_new, t_old(right), 2*t_old(right)-t_old(right-1)
                call stop_mpi()
            endif
        endif
        if (t_new<t_old(1)) then    ! check if t_new value is too small
            write(*,*) 't_new value smaller than smallest value of t_old array', t_new, t_old(1)
            call stop_mpi()
        endif
        do while (right-left > 1)   ! divide search array in halfs until value is found
            middle=(right+left)/2

            if (t_old(middle) > t_new) then
                right=middle
            else if (t_old(middle) < t_new) then
                left=middle
            else
                ind=middle
                exit
            endif
        enddo
        if (ind/=middle) ind=left   ! set ind to the found index value
    end subroutine findindex

    subroutine getMisfit(misfit,mesh,myrank,upperfreq,i_src)
    ! caluclates the misfit between two functions based on the squared difference of these functions
        type(meshVar) :: mesh
        type(recVar) :: rec
        real(kind=CUSTOM_REAL) :: misfit, deltat, upperfreq
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: seismo_model, seismo_meas,t_model,t_meas
        real(kind=SIZE_DOUBLE), dimension(:), allocatable :: seismo_interp
        logical :: file_exists
        integer :: i,j,l, myrank, ntmeas, ntmodel, ier, i_src, length
        character(len=1), dimension(:), allocatable :: dimension
        character(len=1) :: letter
        character(len=256) ::filename, line
        integer :: frac

        frac=10     ! fraction of seismogram which is tappered
        length=2    ! number of seismogram components
        letter='u'  ! last letter of the seismograms marking displacement, velocity or acceleration seismogram. Here: displacement

        allocate(dimension(length))
        misfit=0.0
        if (mesh%has_rec) then
            write(filename,"('/recVar',i6.6)") myrank+1 ! Read receiver information
            filename=trim(outpath)//trim(filename)
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readRecVar(rec,filename)
            else
                write(*,*) "error in recVar, files not existing"
                call stop_mpi()
            end if
            dimension(1)='x'        ! define the two used dimensions. The program is written in a way that it can easily be adapted to other seismogram types
            dimension(2)='z'
            do l=1,length
                do i=1,rec%nrec
                    !obtain number of samples in the simulated seismograms
                    write(filename,"('/seismo.',A1,'.',i7.7,'.sd',A1)") dimension(l),rec%recnr(i),letter
                    filename=trim(outpath)//trim(filename)
                    open(unit=31,file=trim(filename),status='old')
                    ntmodel=0
                    ier=0
                    do while (.not. is_iostat_end(ier))
                        read(31,*,iostat=ier) line
                        if (is_iostat_end(ier)) exit
                        ntmodel = ntmodel + 1
                    enddo
                    !obtain number of samples in the measured seismograms
                    write(filename,"('seismo.',A1,'.',i7.7,'.sd',A1,'.',i3.3,'.meas')") dimension(l),rec%recnr(i),letter, i_src
                    filename=trim(invpath)//trim(filename)
                    open(unit=32,file=trim(filename),status='old')
                    ntmeas=0
                    ier=0
                    do while (.not. is_iostat_end(ier))
                        read(32,*,iostat=ier) line
                        if (is_iostat_end(ier)) exit
                        ntmeas = ntmeas + 1
                    enddo
                    !reset pointer in both seismogram files
                    rewind(31)
                    rewind(32)

                    allocate(seismo_model(ntmodel),seismo_meas(ntmeas),t_model(ntmodel),t_meas(ntmeas))
                    ! read both seismograms
                    do j=1,ntmodel
                        read(31,*) t_model(j), seismo_model(j)
                    enddo
                    do j=1,ntmeas
                        read(32,*) t_meas(j), seismo_meas(j)
                    enddo
                    ! obtain an interpolation array for the shorter seismogram
                    if (t_model(ntmodel) <= t_meas(ntmeas)) then
                        allocate(seismo_interp(ntmodel))
                        deltat=t_model(2)-t_model(1)
                    else
                        allocate(seismo_interp(ntmeas))
                        deltat=t_meas(2)-t_model(1)
                    endif
                    ! Since measured and modeled seismograms can have different frequency content, low-pass filter the measurements
                    call lowPassZeroPhaseFilter(seismo_meas,ntmeas,t_meas(2)-t_meas(1),upperfreq,.true.,frac)
                    ! Depending on which seismogram contains longer times, interpolate one seismogram onto the time sampling of the other
                    if (t_model(ntmodel) <= t_meas(ntmeas)) then
                        call cubicInterpolation(dble(t_meas),dble(seismo_meas),dble(t_model),seismo_interp)
                        do j=1,int(ntmodel*((frac-1.)/frac))    ! calculate misfit but without values in the tapered zone
                            misfit=misfit + 0.5 * deltat * (seismo_model(j)-real(seismo_interp(j), kind=CUSTOM_REAL))**2
                        enddo
                    else
                        call cubicInterpolation(dble(t_model),dble(seismo_model),dble(t_meas),seismo_interp)
                        do j=1,int(ntmeas*((frac-1.)/frac))
                            misfit=misfit + 0.5 * deltat * (seismo_meas(j)-real(seismo_interp(j), kind=CUSTOM_REAL))**2
                        enddo
                    endif

                    close(31)
                    close(32)

                    deallocate(seismo_model,seismo_meas,t_model,t_meas,seismo_interp)
                enddo
            enddo
            deallocate(dimension)
            call deallocRecVar(rec)
        end if
    end subroutine getMisfit

    subroutine makeAdjSourceElastic(par,myrank,inv,inv_step,i_src)
    ! Calculates the source time function for adjoint sources based on the difference between measured and modeled seismograms
        type(parameterVar) :: par
        type(meshVar) :: mesh
        type(recVar) :: rec
        type(srcVar) :: src
        type(time_Series), dimension(:), allocatable ::timeseries_x, timeseries_z
        type(invVar) :: inv
        type(error_message) :: errmsg
        integer :: myrank, inv_step, ier, nt, nt_synth, ntmin, i_src
        character(len=256) ::filename
        logical :: file_exists
        integer :: i,j, ind
        real(kind=CUSTOM_REAL) :: deltat, dummy
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: measured_seis, synth_seis, adj_source, adj_temp
        real(kind=SIZE_DOUBLE), dimension(:,:), allocatable :: dble_prec_seis
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: synth_seis_temp, meas_seis_temp
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: time, time_synth, timemin, time_synth_new

        write(filename,"('/meshVar',i6.6)") myrank+1            ! read mesh information
        filename=trim(outpath)//trim(filename)
        inquire(file=trim(filename), exist=file_exists)
        if (file_exists) then
            call readMeshVar(mesh,filename, errmsg)
        else
            write(*,*) "error in databases, files not existing"
            call stop_mpi()
        end if

        deltat=mesh%dtfactor*par%cfl            ! get current time step

        if (mesh%has_rec) then                  ! adjoint sources can only appear if stations are present
            write(filename,"('/recVar',i6.6)") myrank+1     ! Get information on stations
            filename=trim(outpath)//trim(filename)
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readRecVar(rec,filename)
            else
                write(*,*) "error in recVar, files not existing"
                call stop_mpi()
            end if

            allocate(timeseries_x(rec%nrec), timeseries_z(rec%nrec))    ! allocate two arrays for every stations
            do i=1,rec%nrec
                write(filename,"('/seismo.x.',i7.7,'.sdu')") rec%recnr(i)   ! get sample number of modeled seismograms
                filename=trim(outpath)//trim(filename)
                open(unit=27,file=trim(filename),status='unknown')
                j=0
                ier=0
                do while (.not. is_iostat_end(ier))
                    read(27,*,iostat=ier) dummy, dummy
                    if (is_iostat_end(ier)) exit
                    j = j + 1
                enddo
                close(27)
                nt_synth=j
                write(filename,"('/seismo.x.',i7.7,'.sdu')") rec%recnr(i)   ! read modeled seismograms
                filename=trim(outpath)//trim(filename)
                open(unit=27,file=trim(filename),status='unknown')
                write(filename,"('/seismo.z.',i7.7,'.sdu')") rec%recnr(i)
                filename=trim(outpath)//trim(filename)
                open(unit=28,file=trim(filename),status='unknown')

                allocate(time_synth(nt_synth), synth_seis_temp(nt_synth,2))
                do j=1,nt_synth
                    read(27,*) time_synth(j),synth_seis_temp(j,1)
                    read(28,*) time_synth(j),synth_seis_temp(j,2)
                enddo

                close(27)
                close(28)

                write(filename,"('seismo.x.',i7.7,'.sdu.',i3.3,'.meas')") rec%recnr(i), i_src   ! get number of samples for measured seismograms
                filename=trim(invpath)//trim(filename)
                open(unit=27,file=trim(filename),status='unknown')
                j=0
                ier=0
                do while (.not. is_iostat_end(ier))
                    read(27,*,iostat=ier) dummy, dummy
                    if (is_iostat_end(ier)) exit
                    j = j + 1
                enddo
                close(27)
                nt=j
                allocate(time(nt), meas_seis_temp(nt,2))

                write(filename,"('seismo.x.',i7.7,'.sdu.',i3.3,'.meas')") rec%recnr(i), i_src   ! read measured seismograms
                filename=trim(invpath)//trim(filename)
                open(unit=27,file=trim(filename),status='unknown')
                write(filename,"('seismo.z.',i7.7,'.sdu.',i3.3,'.meas')") rec%recnr(i), i_src
                filename=trim(invpath)//trim(filename)
                open(unit=28,file=trim(filename),status='unknown')

                do j=1,nt
                    read(27,*) time(j),meas_seis_temp(j,1)
                    read(28,*) time(j),meas_seis_temp(j,2)
                enddo
                close(27)
                close(28)

                ntmin=nt_synth                          ! use sample number of modeled seismograms for arrays regarding adjoint sources
                allocate(adj_source(ntmin,rec%nrec,2), adj_temp(ntmin,rec%nrec,2), timemin(ntmin))
                timemin=time_synth

                if (time(nt) < time_synth(nt_synth)) then   ! if measured seismograms are shorter, shorten the modeled seismograms
                    call findindex(time_synth, time(nt), ind)
                    allocate(synth_seis(ind,rec%nrec,2),time_synth_new(ind),measured_seis(ind,rec%nrec,2),dble_prec_seis(ind,2))
                    synth_seis(:,i,:)=synth_seis_temp(1:ind,:)
                    time_synth_new=time_synth(1:ind)
                else
                    ind=ntmin
                    allocate(synth_seis(ntmin,rec%nrec,2),time_synth_new(ntmin),measured_seis(ntmin,rec%nrec,2), dble_prec_seis(ntmin,2))
                    synth_seis(:,i,:)=synth_seis_temp(:,:)
                    time_synth_new=time_synth
                endif

                ! interpolate the measured seismogram onto the time samples of the modeled seismogram
                call cubicInterpolation(dble(time),dble(meas_seis_temp(:,1)),dble(time_synth_new),dble_prec_seis(:,1))
                call cubicInterpolation(dble(time),dble(meas_seis_temp(:,2)),dble(time_synth_new),dble_prec_seis(:,2))
                measured_seis(:,i,1)=real(dble_prec_seis(:,1), kind=CUSTOM_REAL)
                measured_seis(:,i,2)=real(dble_prec_seis(:,2), kind=CUSTOM_REAL)

                deallocate (synth_seis_temp, meas_seis_temp)

                ! Since the measured seismograms can contain higher frequencies as the modeled seismograms, apply a low-pass filter
                call lowPassZeroPhaseFilter(measured_seis(:,i,1),ind,timemin(2)-timemin(1),inv%upperfreq,.true.,10)
                call lowPassZeroPhaseFilter(measured_seis(:,i,2),ind,timemin(2)-timemin(1),inv%upperfreq,.true.,10)

                ! write out filtered measured seismograms to compare them later to the modeled seismograms
                write(filename,"('seismo.x.',i7.7,'.sdu.',i3.3,'.meas.filter',i4.4)") rec%recnr(i),i_src,int(inv%upperfreq)
                filename=trim(invpath)//trim(filename)
                open(unit=27,file=trim(filename),status='unknown')
                write(filename,"('seismo.z.',i7.7,'.sdu.',i3.3,'.meas.filter',i4.4)") rec%recnr(i),i_src,int(inv%upperfreq)
                filename=trim(invpath)//trim(filename)
                open(unit=28,file=trim(filename),status='unknown')
                do j=1,ind
                    write(27,*) timemin(j), measured_seis(j,i,1)
                    write(28,*) timemin(j), measured_seis(j,i,2)
                end do
                close(27)
                close(28)

                ! obtain adjoint sources
                do j=1,ind
                    adj_source(j,i,:)=measured_seis(j,i,:)-synth_seis(j,i,:)
                enddo
                ! if the modeled seismograms are longer than the measured ones, add zeros to the adjoint source
                if (time(nt) <= time_synth(nt_synth)) then
                    adj_source(ind+1:,i,:)=0.
                endif
                ! Reverse adjoint sources as they need to be modeled backwards in time
                do j=1,ntmin
                    adj_temp(ntmin-j+1,i,1)=adj_source(j,i,1)
                    adj_temp(ntmin-j+1,i,2)=adj_source(j,i,2)
                enddo

                adj_source=adj_temp

                ! Finally low-pass filter and taper the adjoint sources again
                call lowPassZeroPhaseFilter(adj_source(:,i,1),ntmin,timemin(2)-timemin(1),inv%upperfreq,.true.,10)
                call lowPassZeroPhaseFilter(adj_source(:,i,2),ntmin,timemin(2)-timemin(1),inv%upperfreq,.true.,10)

                ! write out the adjoint sources
                write(filename,"('adj_stf',i7.7,'iter',i3.3,'src',i3.3)") rec%recnr(i), inv_step, i_src
                filename=trim(adjpath)//trim(filename)
                open(unit=27,file=trim(filename),status='unknown')
                do j=1,ntmin
                    write(27,*) timemin(j), adj_source(j,i,1), adj_source(j,i,2)
                end do
                close(27)
                deallocate(time_synth, time, measured_seis, synth_seis, adj_source, adj_temp, timemin, dble_prec_seis, time_synth_new)
            enddo
            deallocate(timeseries_x, timeseries_z)
        endif
        if (mesh%has_src) call deallocSrcVar(src)
        if (mesh%has_rec) call deallocRecVar(rec)
    end subroutine makeAdjSourceElastic

    subroutine collectKernels(myrank, inv_step,mesh, n_src)
    !subroutine collectKernels(myrank, inv_step,mesh, n_src, inv, chi_model)
    ! Collect the kernels from multiple sources to a global kernel
        integer :: myrank, inv_step, n_src, i_src
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: Kvs, Kvp, model_cp, model_cs
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: Kvs_src, Kvp_src
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rec_vp, rec_vs
        !real(kind=CUSTOM_REAL) ::  chi_model
        type(MeshVar) :: mesh
        !type(InvVar) :: inv
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: maskdummy

        ! The smoothing part is in experimental status and not included in the program yet
        !if (inv%use_smoothing) allocate(rec_vp(mesh%mpi_ne,mesh%mpi_nn), rec_vs(mesh%mpi_ne,mesh%mpi_nn))
        !if (inv%use_smoothing) allocate(model_cp(mesh%nelem), model_cs(mesh%nelem))

        allocate(Kvs(mesh%nelem),Kvp(mesh%nelem), maskdummy(mesh%nelem))
        allocate(Kvs_src(mesh%nelem,n_src),Kvp_src(mesh%nelem,n_src))
        Kvs=0.
        Kvp=0.
        maskdummy=1.
        do i_src=1, n_src   ! just read the different kernels from disk and sum them up
            call readKernel(Kvs_src(:,i_src), Kvp_src(:,i_src),myrank,inv_step,i_src,maskdummy)
            Kvs=Kvs + Kvs_src(:,i_src)
            Kvp=Kvp + Kvp_src(:,i_src)
        enddo
        !        if (inv%use_smoothing) then
        !            call MPI_Send_Elem(mesh%vp, rec_vp, mesh, myrank)
        !            call MPI_Send_Elem(mesh%vs, rec_vs, mesh, myrank)
        !            call get_modelmisfit(myrank, mesh, inv, rec_vp, rec_vs, chi_model)
        !            call sync_mpi()
        !            if (myrank==0) write(*,*) 'finished misfit'
        !            call get_modelsmoothing(myrank, mesh, inv, rec_vp, rec_vs, model_cp, model_cs)
        !            if (myrank==0) write(*,*) 'finished smoothing'
        !            ! add Modelmisfit
        !        endif

        ! Write out the summed kernels
        call writePlotBinariesInv(Kvp, adjpath, 'Kvp', myrank, inv_step)
        call writePlotBinariesInv(Kvs, adjpath, 'Kvs', myrank, inv_step)

        deallocate(Kvs, Kvp, Kvs_src, Kvp_src, maskdummy)
        if (allocated(rec_vp)) deallocate(rec_vp, rec_vs, model_cp, model_cs)
    end subroutine collectKernels

    subroutine readKernel(Kvs,Kvp,myrank,iter_step,i_src,srcrecmask)
    ! Read kernels from file
        real(kind=CUSTOM_REAL), dimension(:) :: Kvs,Kvp,srcrecmask
        integer :: myrank, iter_step, i_src
        character(len=256) ::filename
        type(error_message) :: errmsg

        if (i_src==0) then  ! i_src is used as the number of source and at the same time i_src==0 is used as flag to read an already summed kernel
            call readPlotBinariesInv(1, size(Kvp), Kvp, adjpath, 'Kvp', myrank+1, iter_step, errmsg)
            call readPlotBinariesInv(1, size(Kvs), Kvs, adjpath, 'Kvs', myrank+1, iter_step, errmsg)
        else
            write(filename,"('Kvp','_src',i3.3,'_',i6.6,'_iter',i3.3,'.bin')") i_src,myrank+1,iter_step
            filename=trim(adjpath)//trim(filename)
            open(unit=30,file=trim(filename),form = "UNFORMATTED",status='unknown')
            read(30) Kvp
            close(30)
            write(filename,"('Kvs','_src',i3.3,'_',i6.6,'_iter',i3.3,'.bin')") i_src,myrank+1,iter_step
            filename=trim(adjpath)//trim(filename)
            open(unit=30,file=trim(filename),form = "UNFORMATTED",status='unknown')
            read(30) Kvs
            close(30)
        endif
        ! Apply the mask to reduce the kernel values around sources, receivers and free boundaries
        Kvs=Kvs*srcrecmask
        Kvp=Kvp*srcrecmask
    end subroutine readKernel

    subroutine getSteepestDescent(cvs,cvp,mesh,inv_step, myrank)
    ! Get search direction of steepest descent just be reading the kernel and directly using it
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: Kvs,Kvp,cvs,cvp
        integer :: inv_step, myrank
        type(meshVar) :: mesh

        allocate(Kvs(mesh%nelem),Kvp(mesh%nelem))
        call readKernel(Kvs, Kvp,myrank,inv_step,0,mesh%srcrecmask)
        cvs=Kvs
        cvp=Kvp
        deallocate(Kvs,Kvp)
    end subroutine getSteepestDescent

    subroutine getConjugateGradient(cvs,cvp,mesh,inv_step,myrank)
    ! Get conjugate gradient as search direction
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: Kvs,Kvp,cvs,cvp
        real(kind=SIZE_DOUBLE) :: beta_top_temp, beta_bottom_temp, beta_top, beta_bottom, beta
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: Kvs_old,Kvp_old
        integer :: i, inv_step,myrank
        type(meshVar) :: mesh

        if (inv_step == 1) then ! Use steepest descent in first iteration
            allocate(Kvs(mesh%nelem),Kvp(mesh%nelem))
            call readKernel(Kvs, Kvp,myrank,inv_step,0,mesh%srcrecmask)
            cvs=Kvs
            cvp=Kvp
            deallocate(Kvs,Kvp)
        else
            allocate(Kvs(mesh%nelem),Kvp(mesh%nelem),Kvs_old(mesh%nelem),Kvp_old(mesh%nelem))
            call readKernel(Kvs, Kvp,myrank,inv_step,0,mesh%srcrecmask) ! read current gradient
            call readKernel(Kvs_old, Kvp_old,myrank,inv_step-1,0,mesh%srcrecmask)   ! read gradient from former iteration step
            beta_top_temp=0.
            beta_bottom_temp=0.
            do i=1,mesh%nelem   ! Get beta value to scale conjugate gradient for every rank (Eq. 2.127 Lamert (2019))
                beta_top_temp= beta_top_temp - Kvp(i)*(Kvp_old(i)-Kvp(i))
                beta_top_temp= beta_top_temp - Kvs(i)*(Kvs_old(i)-Kvs(i))
                beta_bottom_temp=beta_bottom_temp + Kvp(i)*Kvp(i)
                beta_bottom_temp=beta_bottom_temp + Kvs(i)*Kvs(i)
            enddo
            call sum_real_allDP(beta_top_temp,beta_top) ! Sum up the different values for beta for every rank
            call sum_real_allDP(beta_bottom_temp,beta_bottom)
            beta=beta_top/beta_bottom
            if (beta < 0.) beta=0.  ! reset to steepest descent if beta is negative
            cvs=Kvs+real(beta, kind=CUSTOM_REAL)*cvs        ! get search direction (Eq. 2.126 Lamert(2019))
            cvp=Kvp+real(beta, kind=CUSTOM_REAL)*cvp
            deallocate(Kvs,Kvp,Kvs_old,Kvp_old)
        endif
    end subroutine getConjugateGradient

    subroutine getBFGS(cvs,cvp,BFGS_p_vp, BFGS_p_vs,BFGS_y_vp, BFGS_y_vs, mesh,inv_step,BFGS_j,myrank,m)
    ! Get search direction from BFGS method based on the last m iterations
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: cvs, cvp, q_vp, q_vs
        real(kind=SIZE_DOUBLE), dimension(:), allocatable :: alpha, rho
        real(kind=CUSTOM_REAL), dimension(:,:) :: BFGS_p_vp, BFGS_p_vs, BFGS_y_vp, BFGS_y_vs
        real(kind=SIZE_DOUBLE) :: sum_sgl, sum_all, temp, gamma
        integer :: inv_step,myrank,i,m, start, j, BFGS_j
        type(meshVar) :: mesh

        allocate (q_vp(mesh%nelem), q_vs(mesh%nelem))
        call readKernel(cvs, cvp,myrank,inv_step,0,mesh%srcrecmask)

        q_vp=-cvp
        q_vs=-cvs
        start=min(BFGS_j,m)-1   ! BFGS uses the last m iterations. However, if there are only less iterations all of them are used

        allocate(alpha(start),rho(start))
        do i=start,1,-1     ! Do for the last iterations but backwards

            sum_sgl=0.
            do j=1,mesh%nelem       ! get rho of Eq. 2.130 Lamert (2019)
                sum_sgl = sum_sgl + BFGS_y_vp(j,i)*BFGS_p_vp(j,i) + BFGS_y_vs(j,i)*BFGS_p_vs(j,i)
            enddo
            call sum_real_allDP(sum_sgl,sum_all)
            rho(i)=1./sum_all

            sum_sgl=0.
            do j=1,mesh%nelem       ! get gamma of Eq. 2.131 Lamert (2019), here it is named alpha!
                sum_sgl = sum_sgl +q_vp(j)*BFGS_p_vp(j,i) + q_vs(j)*BFGS_p_vs(j,i)
            enddo
            call sum_real_allDP(sum_sgl,sum_all)
            alpha(i)=rho(i)*sum_all

            do j=1,mesh%nelem       ! Obtain vector q of Eq. 2.131 Lamert (2019)
                q_vp(j) = q_vp(j)-real(alpha(i), kind=CUSTOM_REAL)*BFGS_y_vp(j,i)
                q_vs(j) = q_vs(j)-real(alpha(i), kind=CUSTOM_REAL)*BFGS_y_vs(j,i)
            enddo
        enddo

        if (start>0) then
            sum_sgl=0.
            do j=1,mesh%nelem   ! Obtain numerator and denominator of Eq. 2.132 Lamert (2019)
                sum_sgl = sum_sgl + BFGS_y_vp(j,start)*BFGS_p_vp(j,start) + BFGS_y_vs(j,start)*BFGS_p_vs(j,start)
            enddo
            call sum_real_allDP(sum_sgl,sum_all)
            temp=sum_all
            sum_sgl=0.
            do j=1,mesh%nelem
                sum_sgl = sum_sgl + BFGS_y_vp(j,start)*BFGS_y_vp(j,start) + BFGS_y_vs(j,start)*BFGS_y_vs(j,start)
            enddo
            call sum_real_allDP(sum_sgl,sum_all)
            gamma=temp/sum_all  ! Gamma are the diagonal entries of H_j^0 Eq. 2.132 Lamert (2019)

            q_vp = real(gamma, kind=CUSTOM_REAL)*q_vp   ! Eq. 2.133 Lamert (2019)
            q_vs = real(gamma, kind=CUSTOM_REAL)*q_vs
        endif

        do i=1,start
            sum_sgl=0.
            do j=1,mesh%nelem   ! Obtain beta of Eq. 2.134 Lamert (2019), here named gamma
                sum_sgl = sum_sgl + q_vp(j)*BFGS_y_vp(j,i) + q_vs(j)*BFGS_y_vs(j,i)
            enddo
            call sum_real_allDP(sum_sgl,sum_all)
            gamma=rho(i)*sum_all

            q_vp(:) = q_vp(:) + BFGS_p_vp(:,i) * real(alpha(i) - gamma, kind=CUSTOM_REAL) ! Apply the second line of Eq. 2.134 Lamert (2019)
            q_vs(:) = q_vs(:) + BFGS_p_vs(:,i) * real(alpha(i) - gamma, kind=CUSTOM_REAL)
        enddo
        cvp=-q_vp   ! Use negative q as search direction
        cvs=-q_vs
        deallocate(q_vp, q_vs,alpha,rho)

    end subroutine getBFGS

    subroutine getBFGS_y(BFGS_y_vp, BFGS_y_vs,mesh,inv_step,BFGS_j,myrank,m)
    ! Obtain the new entry for the y vector of the BFGS method (Eq. 2.130 Lamert (2019))
        real(kind=CUSTOM_REAL), dimension(:,:) :: BFGS_y_vp, BFGS_y_vs
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: Kvs, Kvp, Kvs_old, Kvp_old
        integer :: inv_step,myrank, j, m, ind, BFGS_j
        type(meshVar) :: mesh

        allocate(Kvs(mesh%nelem),Kvp(mesh%nelem),Kvs_old(mesh%nelem),Kvp_old(mesh%nelem))
        ! Read kernel from current and former iteration step
        call readKernel(Kvs, Kvp,myrank,inv_step,0,mesh%srcrecmask)
        call readKernel(Kvs_old, Kvp_old,myrank,inv_step-1,0,mesh%srcrecmask)

        ind=min(BFGS_j,m)-1
        do j=1,mesh%nelem
            BFGS_y_vp(j,ind)=Kvp_old(j)-Kvp(j)  ! Just use the difference as in the equation
            if (mesh%vs(j)>0.) then             ! If the material of an element is a fluid, no search directions needs to be caluclated for this element for vs
                BFGS_y_vs(j,ind)=Kvs_old(j)-Kvs(j)
            else
                BFGS_y_vs(j,ind)=0.
            endif
        enddo

        deallocate(Kvs, Kvp, Kvs_old, Kvp_old)

    end subroutine getBFGS_y

    subroutine changeModel(alpha, mesh, par, cur_BFGS_p_vp, cur_BFGS_p_vs, cvs, cvp, myrank, inv, inv_step, BFGS_j, run_number, errmsg)
    ! Apply the search direction to the current model and obtain the model for the next iteration step
        real(kind=CUSTOM_REAL) :: alpha
        type(MeshVar) :: mesh, meshtemp
        type(InvVar) :: inv
        type(parameterVar) :: par
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:) :: cvs,cvp, cur_BFGS_p_vp, cur_BFGS_p_vs
        integer :: i, myrank, inv_step, run_number, BFGS_j
        character(len=256) ::filename
        logical :: maxreached
        real(kind=CUSTOM_REAL) :: temp, rtemp, temp2
        real(kind=CUSTOM_REAL) :: maxvs, maxvp, maxvstemp, maxvptemp, vpmaxtemp, vpmax, maxkernel
        real(kind=CUSTOM_REAL) :: min_vp_vs_temp, min_vp_vs

        call copyMesh(mesh, meshtemp)   ! get a temporal second mesh variable
        maxkernel = 0.

        if (BFGS_j==1) then             ! In the first step of the BFGS method or in case of the two other methods no scaling of the seach
                                        ! direction to the model values can be obtained. Thus the maximum value of the search direction is
                                        ! scaled to the current model
            maxvptemp=0
            maxvstemp=0

            do i=1,meshtemp%nelem       ! Find the absolute maximum of the search direction
                if (abs(cvs(i)) > maxvstemp) maxvstemp=abs(cvs(i))
                if (abs(cvp(i)) > maxvptemp) maxvptemp=abs(cvp(i))
            enddo

            if (par%nproc > 1) then ! get the maximum of all ranks
                call maxval_real_all(maxvptemp,maxvp,CUSTOM_REAL)
                call maxval_real_all(maxvstemp,maxvs,CUSTOM_REAL)
            else
                maxvp=maxvptemp
                maxvs=maxvstemp
            endif
            maxkernel=max(maxvs,maxvp)
        endif

        min_vp_vs_temp=100.
        do i=1,meshtemp%nelem
            maxreached=.false.
            if (BFGS_j==1) then     ! do the relative change using the maximum of the search direction for the first step of BFGS and for CG and SD
                if (run_number==3) then !get the new entry for the BFGS vector q (Eq 2.130 Lamert (2019))
                    cur_BFGS_p_vp(i)=meshtemp%vp(i) * alpha * cvp(i) / maxkernel
                    cur_BFGS_p_vs(i)=meshtemp%vs(i) * alpha * cvs(i) / maxkernel
                endif
                meshtemp%vp(i) = meshtemp%vp(i) * (1 + alpha * cvp(i) / maxkernel)  ! Update vp and vs
                if (meshtemp%vs(i)>0.) meshtemp%vs(i) = meshtemp%vs(i)  * (1 + alpha * cvs(i) / maxkernel)
            else                    ! For later BFGS steps use the scaling from the BFGS method
                if (run_number==3) then
                    cur_BFGS_p_vp(i)=alpha * cvp(i)
                    cur_BFGS_p_vs(i)=alpha * cvs(i)
                endif
                meshtemp%vp(i) = meshtemp%vp(i) + alpha * cvp(i)
                if (meshtemp%vs(i)>0.) meshtemp%vs(i) = meshtemp%vs(i)  + alpha * cvs(i)
            endif

            if (meshtemp%vp(i) > inv%max_vp) then   ! Check if the maximum or minimum values for vp and vs are reached
                meshtemp%vp(i) = inv%max_vp
                maxreached=.true.
            endif
            if (meshtemp%vp(i) < inv%min_vp) then
                meshtemp%vp(i) = inv%min_vp
                maxreached=.true.
            endif
            if (meshtemp%vs(i) > inv%max_vs) then
                meshtemp%vs(i) = inv%max_vs
                maxreached=.true.
            endif
            if (meshtemp%vs(i) < inv%min_vs .and. meshtemp%vs(i) < epsilon(meshtemp%vs(i))) then
                meshtemp%vs(i) = inv%min_vs
                maxreached=.true.
            endif
            if (meshtemp%vs(i)>0.) then             ! Check the vo/vs ratio
                if (meshtemp%vp(i)/meshtemp%vs(i) < 1.3) then
                    meshtemp%vs(i)=meshtemp%vp(i)/1.3
                    maxreached=.true.
                endif
            endif
            if (min_vp_vs_temp>meshtemp%vp(i)/meshtemp%vs(i)) then      ! if smoothing is included later this needs to be moved behind the smoothing
                min_vp_vs_temp=meshtemp%vp(i)/meshtemp%vs(i)
            endif
        enddo

        !        if (inv%smooth_model) call smoothArray(meshtemp%vp, meshtemp%vs, myrank, meshtemp, inv%depth_model, inv%weight_model)

        if (par%nproc > 1) then                     ! get the overall minimum vp/vs ratio
            call minval_real_all(min_vp_vs_temp,min_vp_vs,CUSTOM_REAL)
            if (myrank==0) write(*,"(a24,f4.2,a52)") '| Smallest vp/vs ratio: ', min_vp_vs,'                                                   |'
        else
            write(*,"(a24,f4.2,a52)") '| Smallest vp/vs ratio: ', min_vp_vs_temp,'                                                   |'
        endif

        do i=1,meshtemp%nelem
            if (par%attenuation) then   ! recalculate the anelastic coefficients for every element. Needs to be done for every element here since every element can have different material values
                call calcAnelasticCoefficientsElement(par,meshtemp%qp(i),meshtemp%qs(i),meshtemp%vp(i),meshtemp%vs(i),meshtemp%rho(i),meshtemp%ylambda(:,i),meshtemp%ymu(:,i),meshtemp%lambda(i),meshtemp%mu(i), errmsg)
                meshtemp%vpu(i) = sqrt((meshtemp%lambda(i)+2*meshtemp%mu(i))/meshtemp%rho(i))
                meshtemp%vsu(i) = sqrt(meshtemp%mu(i)/meshtemp%rho(i))
            else
                meshtemp%vpu(i) = meshtemp%vp(i)
                meshtemp%vsu(i) = meshtemp%vs(i)
                meshtemp%mu(i)  = meshtemp%vs(i)**2*meshtemp%rho(i)
                meshtemp%lambda(i) = meshtemp%vp(i)**2*meshtemp%rho(i)-2*meshtemp%mu(i)
            endif
            ! check for NaN
            if (isnan(meshtemp%vp(i)) .or. isnan(meshtemp%vs(i)) .or. isnan(meshtemp%rho(i)) .or. isnan(meshtemp%mu(i)) .or. isnan(meshtemp%lambda(i))) then
                write(*,"(a80)") "|  WARNING!!!                                                                  |"
                write(*,"(a23,i3,a9,i3,a12)") '| Parameter in element ', i, ' of rank ', myrank,' is infinity                             |'
                write(*,"(a80)") '|      rho,       vp,       vs                                                 |'
                write(*,"(a2,f8.3,a2,f8.3,a2,f8.3,a2,f8.3,a2,f8.3,a30)") '| ',meshtemp%rho(i),', ', meshtemp%vp(i),', ', meshtemp%vs(i),', ',meshtemp%mu(i), ', ',meshtemp%lambda(i), '                             |'
                call stop_mpi()
            endif
        enddo

        ! recalculate time step since the maximum velocity might have changed
        rtemp=1e7
        vpmaxtemp=0.
        do i=1,meshtemp%nelem
            temp=2.0/3.0 * meshtemp%minGLL*(meshtemp%sDT(i)/meshtemp%vpu(i))
            if (temp<rtemp) rtemp=temp
            if (vpmaxtemp<meshtemp%vpu(i)) vpmaxtemp=meshtemp%vpu(i)
        end do

        meshtemp%dtfactor = rtemp
        if (par%nproc > 1) then
            call minval_real_all(meshtemp%dtfactor, temp2,CUSTOM_REAL)
            if (par%autodt) meshtemp%dtfactor=temp2
            call maxval_real_all(vpmaxtemp, vpmax,CUSTOM_REAL)
        endif

        if (par%autodt) then
            if (myrank==0) then
                write(*,"(a17,es12.5,a25,f8.3,a18)") "| New time step: ", meshtemp%dtfactor*par%cfl, ', with highest velocity: ', vpmax, '                 |'
            endif
        else
            if (myrank==0) then
                write(*,"(a18,es12.5,a24,f8.3,a22)") "| Used time step: ", par%dt, ', with highest velocity:', vpmax, '                     |'
                write(*,"(a80)") "| Please check your chosen time step for stability with the new velocity!      |"
            endif
        endif
        ! write out new mesh variables
        write(filename,"('/meshVar',i6.6)") myrank+1
        filename=trim(outpath)//trim(filename)
        call writeMeshVar(meshtemp,filename)

        if (run_number==3) then ! write out the new model as vtk files for the final model change for every iteration step
            call writePlotBinariesInv(meshtemp%vp, invpath, 'model_vp', myrank, inv_step)
            call writePlotBinariesInv(meshtemp%vs, invpath, 'model_vs', myrank, inv_step)
            call writePlotBinariesInv(cvp, invpath, 'cvp', myrank, inv_step)
            call writePlotBinariesInv(cvs, invpath, 'cvs', myrank, inv_step)
        endif
        call deallocMeshvar(meshtemp)
    end subroutine changeModel

    subroutine findBestAlpha(alpha,misfit,bestalpha,nanalert,flag)
    ! finds the minimum of a quadratic function based on three data points and makes sure that the minimum is a real minimum and that it is located within the data range
        real(kind=CUSTOM_REAL), dimension(:) :: alpha, misfit
        real(kind=CUSTOM_REAL), dimension(2) :: a
        real(kind=CUSTOM_REAL) :: bestalpha
        integer :: nanalert, flag

        ! get first two coefficients of the quadratic function (the constant of the function is not important here)
        a(1)=((misfit(1)-misfit(2))/(alpha(1)-alpha(2))-(misfit(2)-misfit(3))/(alpha(2)-alpha(3)))/(alpha(1)-alpha(3))
        a(2)=(misfit(1)-misfit(2))/(alpha(1)-alpha(2))-a(1)*(alpha(1)+alpha(2))

        ! get position of minimum
        bestalpha=-a(2)/(2*a(1))
        nanalert=0
        flag = 0
        if (a(1) > 0.0) then    ! if quadratic function has a minimum
            if (bestalpha > alpha(1) .and. bestalpha < alpha(3)) then   ! minimum is in data range
                if (misfit(3)<=misfit(2)) then  ! not really a quadratic function
                    flag = -2             ! Flag to double data range
                endif
            else if (bestalpha <= alpha(1)) then
                flag = -1                 ! Flag to halve data range
            else
                flag = -2
            endif
            if (isnan(misfit(1)) .or. isnan(misfit(2)) .or. isnan(misfit(3))) then  ! test for nan misfit
                flag = -1
                nanalert = 1
            endif
        else    ! function has a maximum
            if (bestalpha > alpha(1) .and. bestalpha < alpha(3)) then   ! outside data range
                if (misfit(1) <= misfit(3)) then
                    flag = -1
                else
                    bestalpha=alpha(3)      ! if misfit smaller than for the model before, use this large step width
                endif
            else if (bestalpha <= alpha(1)) then
                flag = -2
            else
                flag = -1
            endif
            if (isnan(misfit(1)) .or. isnan(misfit(2)) .or. isnan(misfit(3))) then  ! test for nan misfit
                flag = -1
                nanalert = 1
            endif
        endif

    end subroutine findBestAlpha

    subroutine makeTimeSeriesDP(timeseries, nt, deltat, stf)
    ! Interface function to create a time series from timeSeriesMod for double precision
        type(Time_Series) :: timeseries
        real(kind=8) :: deltat
        real(kind=8), dimension(:) :: stf
        integer :: nt

        call createDPFromDataTimeSeries(timeseries,nt,dble(0.),real(deltat),real(stf))
    end subroutine

    subroutine makeTimeSeriesSP(timeseries, nt, deltat, stf)
    ! Interface function to create a time series from timeSeriesMod for single precision
        type(Time_Series) :: timeseries
        real(kind=4) :: deltat
        real(kind=4), dimension(:) :: stf
        integer :: nt

        call createSPFromDataTimeSeries(timeseries,nt,0.,deltat,stf)
    end subroutine

!   Some functions in experimental state for smoothing
!
!    subroutine MPI_Send_Elem(arr_send, arr_rec, mesh, myrank)
!        real(kind=CUSTOM_REAL), dimension(:) :: arr_send
!        real(kind=CUSTOM_REAL), dimension(:,:) :: arr_rec
!        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: send
!        type(MeshVar) :: mesh
!        integer :: i,j, dest, tag, req, reqrec, myrank
!
!        allocate(send(mesh%mpi_ne,mesh%mpi_nn))
!        do i=1,mesh%mpi_nn
!            do j=1,mesh%mpi_ne ! loop over interface elements
!                if ( mesh%mpi_connection(i,j,1) >0) then
!                    send(j,i)= arr_send(mesh%mpi_connection(i,j,1))   ! sende Werte
!                end if
!            end do
!        end do
!        tag=0
!        do i=1,mesh%mpi_nn
!            dest=mesh%mpi_neighbor(i)-1
!            call isendV_real(send(:,i),mesh%mpi_ne,dest,tag,req,CUSTOM_REAL)
!            call irecV_real(arr_rec(:,i),mesh%mpi_ne,dest,tag,reqrec,CUSTOM_REAL)
!            call wait_req(req)
!            call wait_req(reqrec)
!        end do
!        deallocate(send)
!    end subroutine MPI_Send_Elem
!
!    subroutine get_modelmisfit(myrank, mesh, inv, rec_vp, rec_vs, chi_model)
!        integer :: myrank
!        type(meshVar) :: mesh
!        type(InvVar) :: inv
!        integer :: i, j
!        real(kind=CUSTOM_REAL), dimension(:,:) :: rec_vp, rec_vs
!        real(kind=CUSTOM_REAL) :: chi_model_sngl, chi_model, u
!
!        chi_model_sngl=0.
!        do i=1,mesh%nelem
!            u=mesh%vp(i)*mesh%smooth_A(i,1)
!            if (myrank==0) write(*,*) i, 1, u, mesh%vp(i), mesh%smooth_A(i,:)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) then
!                    u = u + mesh%vp(mesh%smooth_A(i,j))
!                    if (myrank==0) write(*,*) i, j, u, mesh%vp(mesh%smooth_A(i,j)), mesh%smooth_A(i,j), 'normal'
!                endif
!                if (mesh%smooth_A(i,j)<0) then
!                    u = u + rec_vp(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!                    if (myrank==0) write(*,*) i, j, u, rec_vp(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i)), mesh%smooth_A(i,j), 'MPI'
!                endif
!            enddo
!            chi_model_sngl=chi_model_sngl + u**2
!            u=mesh%vs(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) u = u + mesh%vs(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) u = u + rec_vs(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!            chi_model_sngl=chi_model_sngl + u**2
!        enddo
!        call sum_real_all(chi_model_sngl,chi_model,CUSTOM_REAL)
!    end subroutine get_modelmisfit
!
!
!    subroutine get_modelsmoothing(myrank, mesh, inv, rec_vp, rec_vs, model_cp, model_cs)
!        integer :: myrank
!        type(meshVar) :: mesh
!        type(InvVar) :: inv
!        real(kind=CUSTOM_REAL), dimension(:,:) :: rec_vp, rec_vs
!        real(kind=CUSTOM_REAL), dimension(:) :: model_cp, model_cs
!        real(kind=CUSTOM_REAL), dimension(:), allocatable :: u
!        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rec_u
!        integer :: i,j
!
!        allocate(u(mesh%nelem))
!        allocate(rec_u(mesh%mpi_ne,mesh%mpi_nn))
!
!        do i=1,mesh%nelem   ! u=Am
!            u(i)=mesh%vp(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) u(i) = u(i) + mesh%vp(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) u(i) = u(i) + rec_vp(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!        enddo
!        call MPI_Send_Elem(u, rec_u, mesh, myrank)  ! MPI for u
!        do i=1,mesh%nelem   ! cp=Au
!            model_cp(i)=u(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) model_cp(i) = model_cp(i) + u(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) model_cp(i) = model_cp(i) + rec_u(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!        enddo
!        !same for vs
!        do i=1,mesh%nelem   ! u=Am
!            u(i)=mesh%vs(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) u(i) = u(i) + mesh%vs(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) u(i) = u(i) + rec_vs(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!        enddo
!        call MPI_Send_Elem(u, rec_u, mesh, myrank)  ! MPI for u
!        do i=1,mesh%nelem   ! cs=Au
!            model_cs(i)=u(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) model_cs(i) = model_cs(i) + u(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) model_cs(i) = model_cs(i) + rec_u(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!        enddo
!        deallocate(u, rec_u)
!    end subroutine get_modelsmoothing

end module adjointMod
