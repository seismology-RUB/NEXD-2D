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

!-------------------------------------------------------------
!> \brief Routines that operate with time series
!-------------------------------------------------------------
module timeSeries
    use fourierTransform
    use filterCoefficientsDecimate
    use recursiveFilterCoefficients
    use dateTime
    use realloc
    implicit none
    interface dealloc; module procedure deallocTimeSeries; end interface
    interface multiplyTimeSeries
        module procedure scalarMultiplyTimeSeries
        module procedure arrayMultiplyTimeSeries
    end interface
    interface energyTimeSeries
        module procedure energyNormalizedTimeSeries
        module procedure energyUnnormalizedTimeSeries
    end interface
    interface createZeroTimeSeries
        module procedure createSPZeroTimeSeries
        module procedure createDPZeroTimeSeries
        module procedure createWithDateZeroTimeSeries
    end interface
    interface createLinkTimeSeries
        module procedure createSPLinkTimeSeries
        module procedure createDPLinkTimeSeries
        module procedure createWithDateLinkTimeSeries
    end interface
    interface createFromDataTimeSeries
        module procedure createSPFromDataTimeSeries
        module procedure createDPFromDataTimeSeries
        module procedure createWithDateFromDataTimeSeries
    end interface
    interface createEmptyTimeSeries
        module procedure createSPEmptyTimeSeries
        module procedure createDPEmptyTimeSeries
        module procedure createWithDateEmptyTimeSeries
    end interface
    interface adjustTimeSeries
        module procedure adjustSPTimeSeries
        module procedure adjustDPTimeSeries
    end interface
    interface operator (.nsamp.); module procedure nsampTimeSeries; end interface
    interface operator (.dt.); module procedure dtTimeSeries; end interface
    interface operator (.sample.); module procedure getSampleTimeSeries; end interface
    interface operator (.stime.); module procedure getStimeTimeSeries; end interface
    interface operator (.tanf.); module procedure startTimeSeries; end interface
    interface operator (.tanfdp.); module procedure startDPTimeSeries; end interface
    interface operator (.tend.); module procedure endTimeSeries; end interface
    interface operator (.tenddp.); module procedure endDPTimeSeries; end interface
    interface operator (.trace.); module procedure traceTimeSeries; end interface
    interface operator (.tfs.); module procedure startFullSecondsTimeSeries; end interface
    interface operator (.tns.); module procedure startNanoSecondsTimeSeries; end interface
    type time_series
        private
        integer :: nsamp                             !< Number of samples
        real :: dt                                   !< Sampling interval
        double precision :: tanf                     !< Start time in sec after midnight in double precision
        type (date_time) :: stime                    !< Date time object with date and time of first sample
        logical :: exdate                            !< True if stime member was explicitly set
        logical :: link                              !< If true data is just a pointer
        real, dimension(:), pointer :: y => null()   !< Pointer to data
    end type
!
    contains
    !--------------------------------------------------------
    !> \brief Create zero time series
    !> \param this Time series object
    !> \param nsamp Number of samples (space for data will be allocated)
    !> \param tanf Starting time after midnight in seconds
    !> \param dt Sampling interval
    !
    subroutine createDPZeroTimeSeries(this,nsamp,tanf,dt)
        type (time_series) :: this
        integer :: nsamp
        real :: dt
        double precision :: tanf
        allocate(this%y(nsamp))
        this%nsamp = nsamp; this%tanf = tanf; this%dt = dt
        this%link = .false.
        this%exdate = .false.
        this%y = 0.
        call createDateTime(this%stime,1800,1,0,0,0,0)           ! some out of range default value if not specified
        this%stime = (this%stime).plus.(this%tanf)
    end subroutine createDPZeroTimeSeries
    !
    ! for compatibility
    !
    subroutine createSPZeroTimeSeries(this,nsamp,tanf,dt)
        type (time_series) :: this
        integer :: nsamp
        real :: dt,tanf
        call createDPZeroTimeSeries(this,nsamp,dble(tanf),dt)
    end subroutine createSPZeroTimeSeries
    !
    ! with date time object for first sample
    !
    subroutine createWithDateZeroTimeSeries(this,nsamp,stime,dt)
        type (time_series) :: this
        type (date_time) :: stime
        integer :: nsamp
        real :: dt
        call createDPZeroTimeSeries(this,nsamp,.tanfdp.stime,dt)
        this%stime = stime
        this%exdate = .true.
    end subroutine createWithDateZeroTimeSeries
    !--------------------------------------------------------
    !> \brief  create time_series object with just a link to the data
    !
    subroutine createDPLinkTimeSeries(this,nsamp,tanf,dt,s)
        type (time_series) :: this
        integer :: nsamp
        real :: dt
        double precision :: tanf
        real, dimension(:), target :: s
        !
        this%link = .true.
        this%nsamp = nsamp; this%dt = dt; this%tanf = tanf
        this%y => s
        call createDateTime(this%stime,1800,1,0,0,0,0)           ! some out of range default value if not specified
        this%stime = (this%stime).plus.(this%tanf)
        this%exdate = .false.
    end subroutine createDPLinkTimeSeries
    !
    ! for compatibility
    !
    subroutine createSPLinkTimeSeries(this,nsamp,tanf,dt,s)
        type (time_series) :: this
        integer :: nsamp
        real :: dt,tanf
        real, dimension(:), target :: s
        call createDPLinkTimeSeries(this,nsamp,dble(tanf),dt,s)
    end subroutine createSPLinkTimeSeries
    !
    ! with date time object for first sample
    !
    subroutine createWithDateLinkTimeSeries(this,nsamp,stime,dt,s)
        type (time_series) :: this
        type (date_time) :: stime
        real, dimension(:), target :: s
        integer :: nsamp
        real :: dt
        call createDPLinkTimeSeries(this,nsamp,.tanfdp.stime,dt,s)
        this%stime = stime
        this%exdate = .true.
    end subroutine createWithDateLinkTimeSeries
    !--------------------------------------------------------
    !> \brief  create time_series object with a true copy of the data
    !
    subroutine createDPFromDataTimeSeries(this,nsamp,tanf,dt,s)
        type (time_series) :: this
        integer :: nsamp
        real :: dt
        double precision :: tanf
        real, dimension(:) :: s
        !
        this%link = .false.
        this%nsamp = nsamp; this%dt = dt; this%tanf = tanf
        allocate(this%y(nsamp))
        this%y = s(1:nsamp)
        call createDateTime(this%stime,1800,1,0,0,0,0)           ! some out of range default value if not specified
        this%stime = (this%stime).plus.(this%tanf)
        this%exdate = .false.
    end subroutine createDPFromDataTimeSeries
    !
    ! for compatibility
    !
    subroutine createSPFromDataTimeSeries(this,nsamp,tanf,dt,s)
        type (time_series) :: this
        integer :: nsamp
        real :: dt,tanf
        real, dimension(:) :: s
        call createDPFromDataTimeSeries(this,nsamp,dble(tanf),dt,s)
    end subroutine createSPFromDataTimeSeries
    !
    ! with date time object for first sample
    !
    subroutine createWithDateFromDataTimeSeries(this,nsamp,stime,dt,s)
        type (time_series) :: this
        type (date_time) :: stime
        real, dimension(:), target :: s
        integer :: nsamp
        real :: dt
        call createDPFromDataTimeSeries(this,nsamp,.tanfdp.stime,dt,s)
        this%stime = stime
        this%exdate = .true.
    end subroutine createWithDateFromDataTimeSeries
    !-------------------------------------------------------------
    !> \brief Create a time series with empty (unassopciated) data vector
    !> \param this Time series object
    !> \param tanf Start time
    !> \param dt Sampling interval
    !
    subroutine createDPEmptyTimeSeries(this,tanf,dt)
        type (time_series) :: this
        real :: dt
        double precision :: tanf
        this%nsamp = 0; this%tanf = tanf; this%dt = dt; this%link = .false.
        this%y => null()
        call createDateTime(this%stime,1800,1,0,0,0,0)           ! some out of range default value if not specified
        this%stime = (this%stime).plus.(this%tanf)
        this%exdate = .false.
    end subroutine createDPEmptyTimeSeries
    !
    ! for compatibility
    !
    subroutine createSPEmptyTimeSeries(this,tanf,dt)
        type (time_series) :: this
        real :: dt,tanf
        call createDPEmptyTimeSeries(this,dble(tanf),dt)
    end subroutine createSPEmptyTimeSeries
    !
    ! with date time object for first sample
    !
    subroutine createWithDateEmptyTimeSeries(this,nsamp,stime,dt)
        type (time_series) :: this
        type (date_time) :: stime
        integer :: nsamp
        real :: dt
        call createDPZeroTimeSeries(this,nsamp,.tanfdp.stime,dt)
        this%stime = stime
        this%exdate = .true.
    end subroutine createWithDateEmptyTimeSeries
    !-------------------------------------------------------------------------
    !> \brief  absolute maximum of a time series
    !> \param this Time series object
    !
    real function absmaxTimeSeries(this)
        type (time_series) :: this
        absmaxTimeSeries = maxval(abs(this%y(1:this%nsamp)))
    end function absmaxTimeSeries
    !-------------------------------------------------------------------------
    !> \brief  Sample index of absolute maximum of a time series
    !> \param this Time series object
    !
    integer function absmaxLocationTimeSeries(this)
        type (time_series) :: this
        absmaxLocationTimeSeries = maxloc(abs(this%y(1:this%nsamp)),1)
    end function absmaxLocationTimeSeries
    !--------------------------------------------------------
    !> \brief  calculate average of time series
    !
    real function averageTimeSeries(this)
        type (time_series) :: this
        averageTimeSeries = sum(this%y)/this%nsamp
    end function averageTimeSeries
    !--------------------------------------------------------
    !> \brief  Add a scaled time series to a second one
    !> Just add sample by sample without taking care of tanf
    !< check for equal dt and equal length
    !
    subroutine addScaledTimeSeries(this,scale,other)
        type (time_series) :: this,other
        real :: scale
        if (abs((other%dt-this%dt)/this%dt) > tiny(this%dt)) then
            print *,'dt of time series differ!'
            stop
        endif
        if (this%nsamp /= other%nsamp) then
            print *,'nsamp of time series differ!'
            stop
        endif
        other%y = other%y+scale*this%y
    end subroutine addScaledTimeSeries
    !--------------------------------------------------------------------------
    !> \brief Add a constant to a time series sample by sample
    !
    subroutine addConstantTimeSeries(this,add)
        type (time_series) :: this
        real :: add
        this%y = this%y+add
    end subroutine addConstantTimeSeries
    !-----------------------------------------------------------------------------
    !> \brief  adjust time series such that that therd is on a sample
    !>  and that desired begin and end of time series is honoured
    !!  time series can be interpolated to a new higher sample rate
    !> \param  therd source time
    !> \param  tbh	desired time span before therd
    !> \param  tah	desired time span after therd
    !> \param  dtnew new sampling interval
    !> \param  ats adjusted time series
    !-----------------------------------------------------------------------------
    function adjustDPTimeSeries(this,therd,tbh,tah,dtnew,ierr) result(ats)
        type (time_series) :: this,ats
        type (date_time) :: stime
        integer :: ns,j,n
        real :: frac,dtnew
        double precision :: therd,tbh,tah,tbeg,t,tend
        integer, optional :: ierr
        !
        !  check feasibility
        !
        if ((dtnew-this%dt)/this%dt > 1.e-6) then
            if (present(ierr)) then
                ierr = 1; return
            else
                print *,'new sampling interval is larger than old one!'
                print *,'dtnew = ',dtnew,', current dt = ',this%dt
                stop
            endif
        endif
        if (this%tanf.gt.therd-tbh) then
            if (present(ierr)) then
                print *,'Gap: ',this%tanf-therd+tbh
                ierr = 2; return
            else
                print *,'Error: Data start too late.'
                print *,'Gap: ',this%tanf-therd+tbh
                stop
            endif
        endif
        tend = this%tanf+(this%nsamp-1)*this%dt
        if ((tend/(therd+tah)-1.) < -this%dt ) then
            if (present(ierr)) then
                print *,'Missing time range: ',therd+tah-tend
                print *,therd,tah,therd+tah,tend
                ierr = 3; return
            else
                print *,'Error: Data too short'
                print *,'Missing time range: ',therd+tah-tend
                print *,therd,tah,therd+tah,tend
                stop
            endif
        endif
        if (present(ierr)) ierr = 0
        !
        !  beginning of time series,
        !  an integer number of samples before therd and close to therd-tbh
        !
        tbeg=therd-int(tbh/dtnew)*dtnew
        ns = int((therd+tah-tbeg)/dtnew)+1
        !
        !  new time series object
        !
        if (this%exdate) then
            if (tbeg-this%tanf >= 0.d0) then
                stime = (this%stime).plus.(abs(tbeg-this%tanf))
            else
                stime = (this%stime).minus.(abs(tbeg-this%tanf))
            endif
            call createZeroTimeSeries(ats,ns,stime,dtnew)
        else
            call createZeroTimeSeries(ats,ns,tbeg,dtnew)
        endif
        !
        !  interpolate time series at new sample positions starting with tbeg
        !  and ending just before therd+tah
        !
        do j=1,ns
            t=tbeg+(j-1)*dtnew
            n=min( int((t-this%tanf)/this%dt)+1, this%nsamp-1)
            frac = real((t-this%tanf)/this%dt-(n-1),4)
            ats%y(j)=(1.-frac)*this%y(n)+frac*this%y(n+1)
        enddo
    end function adjustDPTimeSeries
    !-----------------------------------------------------------------------
    !  for compatibility
    !
    function adjustSPTimeSeries(this,therd,tbh,tah,dtnew) result(ats)
        type (time_series) :: this,ats
        real :: dtnew
        real :: therd,tbh,tah
        ats = adjustDPTimeSeries(this,dble(therd),dble(tbh),dble(tah),dtnew)
    end function adjustSPTimeSeries
    !------------------------------------------------------------------------------------
    !> \brief Akaike Information Criterium of a time series
    !
    function aicTimeSeries(this) result(res)
        type (time_series) :: this,res,cs
        integer :: k
        cs = cumulativeSumTimeSeries(this)
        where (cs%y < tiny(1.0)) cs%y = tiny(1.0)
        if (this%exdate) then
            call createZeroTimeSeries(res,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(res,this%nsamp,this%tanf,this%dt)
        endif
        do k = 2,this%nsamp
            res%y(k) = (k-1)*log(cs%y(k)/k)+(this%nsamp-k+1)*log((cs%y(this%nsamp)-cs%y(k-1))/(this%nsamp-k+1))
        enddo
        res%y(1) = res%y(2)
        call dealloc(cs)
    end function aicTimeSeries
    !------------------------------------------------------------------------------------
    !> \brief Append data to data vector of time series
    !> Adjust sample number. Do this only for time series with link = .false..
    !
    !> \param this Time series object
    !> \param d vector with additional data
    !
    subroutine appendDataTimeSeries(this,d)
        type (time_series) :: this
        real, dimension(:) :: d
        integer :: n
        !
        if (this%link) then
            print *,'Data not owned by time series object. Do nothing'
            return
        endif
        !
        n = size(d)
        this%y => reallocate(this%y,this%nsamp+n)
        this%y(this%nsamp+1:this%nsamp+n) = d(1:n)
        this%nsamp = this%nsamp+n
    end subroutine appendDataTimeSeries
    !--------------------------------------------------------------------------
    !> \brief Multiply a time series with an array sample by sample
    !
    subroutine arrayMultiplyTimeSeries(this,factor_array)
        type (time_series) :: this
        real, dimension(:) :: factor_array
        integer :: nf,nm
        nf = size(factor_array)
        nm = min(this%nsamp,nf)
        this%y(1:nm) = this%y(1:nm)*factor_array(1:nm)
    end subroutine arrayMultiplyTimeSeries
    !--------------------------------------------------------------------------
    !> \brief Divide a time series by an array sample by sample
    !
    subroutine arrayDivideTimeSeries(this,denom_array)
        type (time_series) :: this
        real, dimension(:) :: denom_array
        integer :: nf,nm
        nf = size(denom_array)
        nm = min(this%nsamp,nf)
        this%y(1:nm) = this%y(1:nm)/denom_array(1:nm)
    end subroutine arrayDivideTimeSeries
    !--------------------------------------------------------
    !> \brief  create a true copy of given time series
    !
    subroutine copyTimeSeries(this,copy,n)
        type (time_series) :: this,copy
        integer, optional :: n
        integer :: nsamp
        if (present(n)) then; nsamp = n; else; nsamp = this%nsamp; endif
        copy%nsamp = nsamp
        copy%dt = this%dt
        copy%tanf = this%tanf
        copy%stime = this%stime
        copy%exdate = this%exdate
        allocate(copy%y(nsamp))
        copy%y = this%y(1:nsamp)                  ! true copy
        copy%link = .false.
    end subroutine copyTimeSeries
    !--------------------------------------------------------------------------------
    !> \brief  causally convolve time series with some filter response
    !>
    !> \param  h filter response
    !> \param  cts convolved time series
    !
    function convolveTimeSeries(this,h) result(cts)
        type (time_series) :: this,cts
        real, dimension(:) :: h
        integer :: nh,j,k,nd,nc
        !
        nd = this%nsamp
        nh = size(h)
        nc = max(nd,nh)
        !
        if (this%exdate) then
            call createZeroTimeSeries(cts,nc,this%stime,this%dt)
        else
            call createZeroTimeSeries(cts,nc,this%tanf,this%dt)
        endif
        !
        !  perform convolution, explanation see causfalt.txt
        !
        if (nd >= nh) then
            do k=1,nc
                cts%y(k) = 0.
                do j=1,min(k,nh)
                    cts%y(k) = cts%y(k) + h(j)*this%y(k-j+1)
                enddo
            enddo
        else
            do k=1,nc
                cts%y(k) = 0.
                do j=1,min(k,nd)
                    cts%y(k) = cts%y(k) + this%y(j)*h(k-j+1)
                enddo
            enddo
        endif
        cts%y = cts%y*this%dt
    end function convolveTimeSeries
    !--------------------------------------------------------------
    !> \brief  Do a cross correlation in the time domain
    !
    !>  If k >= 0: c(k) = sum_{i=1}^{min(ns,nd-k)} d(k+i)*s(i)*dt.
    !>  If k <  0: c(k) = sum_{i=|k|+1}^{min(nd+|k|,ns)} d(k+i)*s(i)*dt
    !>
    !> \param  ncc dimension of c in calling program (input)
    !> \param  nd number of d samples (input)
    !> \param  ns number of s samples (input)
    !> \param  nc number of cross correlation samples (output)
    !> \param  d d-array (input)
    !> \param  s s-array (input)
    !> \param  c cross correlation function (ouptut)
    !> \param  kmin true index of leftmost sample (output)
    !> \param  kmax true index of rightmost sample (output)
    !> \param  kl desired index of leftmost sample of c (minimum = -ns+1) (optional)
    !> \param  kr desired index of rightmost sample of c (maximum = nd-1) (optional)
    !>
    !>  Note: index k=0 corresponds to zero lag
    !>        Maximum number of samples of cross-correlation is nd+ns-1.
    !>        Index k=kmin is mapped to first element of c in calling program,
    !>        Index k=0 is mapped to (first element-kmin) of c in calling program if kmin < 0.
    !>        Else k=0 is not considered.
    !<        Index k=kmax is mapped to (first element+|kmax-kmin|) of c in calling program.
    !-------------------------------------------------------------------
    function crossCorrelateTimeSeries(this,s,kmin,kmax,kl,kr)  result(croco)
        type (time_series) :: this,croco
        integer, optional :: kl, kr
        real, dimension(:) :: s
        real :: help
        integer ::i,k,kmin,kmax,nd,ns,nc
        !
        nd = this%nsamp
        ns = size(s)
        if(.not.present(kl)) then; kmin = -ns+1; else; kmin = max(kl,-ns+1); endif
        if(.not.present(kr)) then; kmax = nd-1;  else; kmax = min(kr, nd-1); endif
        nc = kmax-kmin+1
        call createZeroTimeSeries(croco,nc,dble(kmin*this%dt),this%dt)
        !
        do k=kmin,-1
            help=0.
            do i=iabs(k)+1,min(nd+iabs(k),ns)
                help=help+this%y(k+i)*s(i)
            enddo
            croco%y(k-kmin+1)=help*this%dt
        enddo
        do k=max(kmin,0),kmax
            help=0.
            do i=1,min(nd-k,ns)
                help=help+this%y(k+i)*s(i)
            enddo
            croco%y(k-kmin+1)=help*this%dt
        enddo
    end function crossCorrelateTimeSeries
    !---------------------------------------------------------------------------
    !> \brief Cumulative sum of a time series
    !
    function cumulativeSumTimeSeries(this) result(ts)
        type (time_series) :: this,ts
        integer :: i
        !
        if (this%exdate) then
            call createZeroTimeSeries(ts,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(ts,this%nsamp,this%tanf,this%dt)
        endif
        ts%y(1) = this%y(1)
        do i = 2,this%nsamp
            ts%y(i) = ts%y(i-1)+this%y(i)
        enddo
    end function cumulativeSumTimeSeries
    !--------------------------------------------------------
    !> \brief  free time_series object
    !
    subroutine deallocTimeSeries(this)
        type (time_series) :: this
        if (associated(this%y)) then
            if (this%link) then; nullify(this%y); else; deallocate(this%y); endif
        endif
    end subroutine deallocTimeSeries
    !-------------------------------------------------------------------------
    !> \brief decimate a time series using a given decimation sequence
    !
    function decimateTimeSeries(this,decfac) result(rts)
        type (time_series) :: this,rts,temp
        integer, dimension(:) :: decfac
        integer :: i
        !
        call copyTimeSeries(this,rts)
        do i=1,size(decfac)
            if(decfac(i).eq.2) then
                temp = resampleTimeSeries(rts,filCoeffDec2,decfac(i))
            else if(decfac(i).eq.3) then
                temp = resampleTimeSeries(rts,filCoeffDec3,decfac(i))
            else if(decfac(i).eq.5) then
                temp = resampleTimeSeries(rts,filCoeffDec5,decfac(i))
            else
                print *,'Decimation factor not implemented, use 2, 3 or 5!'
                stop
            endif
            call dealloc(rts)      ! free intermediate storage
            rts = temp
        enddo
    end function decimateTimeSeries
    !------------------------------------------------------------------------
    !> \brief  detrend a time series
    !
    function detrendTimeSeries(this) result(fts)
        type (time_series) :: this,fts
        double precision :: gn,alpha,beta,det,sx,sjx,a,b
        integer :: j
        !
        if (this%exdate) then
            call createZeroTimeSeries(fts,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(fts,this%nsamp,this%tanf,this%dt)
        endif
        gn = this%nsamp
        alpha = 0.5d0*gn*(gn+1.d0)
        beta = (2.d0*gn+1.d0)*(gn+1.d0)*gn/6.d0
        det = gn*beta-alpha*alpha
        sx = 0.d0
        sjx = 0.d0
        do j=1,this%nsamp
            sx = sx+this%y(j)
            sjx = sjx+this%y(j)*j
        enddo
        a = (sx*beta-sjx*alpha)/det
        b = (sjx*gn-sx*alpha)/det
        do j=1,this%nsamp
            fts%y(j) = real(this%y(j)-a-b*j,4)
        enddo
    end function detrendTimeSeries
    !--------------------------------------------------------------------------------
    !> \brief  differentiate time series
    !>
    !> \param  dts differentiated time series (output)
    !>
    !!  dp(i) = (d(i+1)-d(i-1))/(2*dt)
    !!  dp(1) = dp(2)
    !!  dp(n) = dp(n-1)
    !
    function diffTimeSeries(this) result(dts)
        type (time_series) :: this,dts
        integer :: n,i
        !
        n = this%nsamp
        if (this%exdate) then
            call createZeroTimeSeries(dts,n,this%stime,this%dt)
        else
            call createZeroTimeSeries(dts,n,this%tanf,this%dt)
        endif
        do i=2,n-1
            dts%y(i) = (this%y(i+1)-this%y(i-1))/(2.*this%dt)
        enddo
        dts%y(1) = dts%y(2)
        dts%y(n) = dts%y(n-1)
    !
    end function diffTimeSeries
    !-------------------------------------------------------
    !> \brief  get dt of time series
    !
    real function dtTimeSeries(this)
        type (time_series), intent(in) :: this
        dtTimeSeries = this%dt
    end function dtTimeSeries
    !-------------------------------------------------------------------------
    !> \brief end time of a time series (double)
    !
    function endDPTimeSeries(this) result(tend)
        type (time_series), intent(in) :: this
        double precision :: tend
        tend = this%tanf+(this%nsamp-1)*this%dt
    end function endDPTimeSeries
    !-------------------------------------------------------------------------
    !> \brief end time of a time series (single)
    !
    function endTimeSeries(this) result(tend)
        type (time_series), intent(in) :: this
        real :: tend
        tend = sngl(this%tanf)+(this%nsamp-1)*this%dt
    end function endTimeSeries
    !-------------------------------------------------------------------------
    !> \brief  energy of a time series
    !
    real function energyNormalizedTimeSeries(this,sigma)
        type (time_series) :: this
        real, dimension(:) :: sigma
        energyNormalizedTimeSeries = sum((this%y(1:this%nsamp)/sigma(1:this%nsamp))**2)*this%dt
    end function energyNormalizedTimeSeries
    !-------------------------------------------------------------------------
    !> \brief  energy of a time series without normalization
    !
    real function energyUnnormalizedTimeSeries(this)
        type (time_series) :: this
        energyUnnormalizedTimeSeries = sum(this%y(1:this%nsamp)**2)
    end function energyUnnormalizedTimeSeries
    !---------------------------------------------------------------------------
    !> \brief Extend a pointer array of time_series-objects
    !
    function extendArrayTimeSeries(array,n) result(newarray)
        type (time_series), dimension(:), pointer :: array
        type (time_series), dimension(:), pointer :: newarray
        integer :: n,nold,i
        !
        allocate(newarray(n))
        if (.not. associated(array)) return
        nold = min(size(array),n)
        newarray(1:nold) = array(1:nold)
        do i = 1,nold
            call unlinkDataTimeSeries(array(i))
        enddo
        do i = nold+1,size(array)
            call dealloc(array(i))
        enddo
        deallocate(array)
    end function extendArrayTimeSeries
    !-------------------------------------------------------
    !> \brief get i-th sample of time series
    !
    real function getSampleTimeSeries(this,i)
        type (time_series), intent(in) :: this
        integer, intent(in) :: i
        getSampleTimeSeries = this%y(i)
    end function getSampleTimeSeries
    !-------------------------------------------------------------------
    !> \brief Get start time of time series in date_time object
    !
    function getStimeTimeSeries(this) result(stime)
        type (time_series), intent(in) :: this
        type (date_time) :: stime
        stime = this%stime
    end function getStimeTimeSeries
    !-----------------------------------------------------------------------------
    !> \brief Hanning taper of time series at both ends. Use 0.5*(1-cos(pi*t/w))
    !
    function hanningTaperTimeSeries(this,w) result(res)
        type (time_series) :: this
        type (time_series) :: res
        real :: w,tend,t
        integer :: i,nt
        !
        if (this%exdate) then
            call createZeroTimeSeries(res,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(res,this%nsamp,this%tanf,this%dt)
        endif
        nt = int(w/this%dt)
        do i = 1,nt
            t = (i-1)*this%dt
            res%y(i) = this%y(i)*0.5*(1.-cos(mc_pi*t/w))
        enddo
        do i = nt+1,this%nsamp-nt
            res%y(i) = this%y(i)
        enddo
        tend = (this%nsamp-1)*this%dt
        do i = this%nsamp-nt+1,this%nsamp
            t = (i-1)*this%dt
            res%y(i) = this%y(i)*0.5*(1.-cos(mc_pi*(tend-t)/w))
        enddo
    end function hanningTaperTimeSeries
    !-----------------------------------------------------------------------------
    !> \brief Hanning taper of time series at tail only. Use 0.5*(1-cos(pi*(tend-t)/w))
    !
    function hanningTaperTailTimeSeries(this,w) result(res)
        type (time_series) :: this
        type (time_series) :: res
        real :: w,tend,t
        integer :: i,nt
        !
        if (this%exdate) then
            call createZeroTimeSeries(res,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(res,this%nsamp,this%tanf,this%dt)
        endif
        nt = int(w/this%dt)
        do i = 1,this%nsamp-nt
            res%y(i) = this%y(i)
        enddo
        tend = (this%nsamp-1)*this%dt
        do i = this%nsamp-nt+1,this%nsamp
            t = (i-1)*this%dt
            res%y(i) = this%y(i)*0.5*(1.-cos(mc_pi*(tend-t)/w))
        enddo
    end function hanningTaperTailTimeSeries
    !------------------------------------------------------------------------
    !> \brief  high pass Butterworth filtering of time series using recursive filters
    !> \param  fc corner frequency (Hz)
    !> \param  m order
    !>  use a hp1 and a hp2 filter to do the job
    !
    function highPassButterworthRecursiveTimeSeries(this,fc,m) result(fts)
        type (time_series) :: this,fts,temp
        type (recursive_filter_coeff) :: rfchp1,rfchp2
        real :: fc,tc,h
        double precision :: pih
        integer :: m,mm,j
        !
        tc = 1./fc
        call copyTimeSeries(this,fts)
        mm = m/2
        if (m > 2*mm) then                                                  ! odd order filter
            call highPassO1RecursiveFilterCoefficients(rfchp1,tc,this%dt)
            temp = recursiveFilterTimeSeries(fts,rfchp1)
            call dealloc(fts)
            fts = temp
        endif
        !
        pih=0.5d0*pi
        do j=1,mm
            h = real(sin(pih*(2*j-1)/m),4)
            call highPassO2RecursiveFilterCoefficients(rfchp2,tc,h,this%dt)
            temp = recursiveFilterTimeSeries(fts,rfchp2)
            call dealloc(fts)
            fts = temp
        enddo
    end function highPassButterworthRecursiveTimeSeries
    !-----------------------------------------------------------------------------
    !> \brief  integrate time series using trapezoidal rule
    !
    function integrateTimeSeries(this) result(its)
        type (time_series) :: this,its
        integer :: i
        !
        if (this%exdate) then
            call createZeroTimeSeries(its,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(its,this%nsamp,this%tanf,this%dt)
        endif
        !
        its%y(1) = 0.
        do i=2,this%nsamp
            its%y(i) = its%y(i-1)+0.5*(this%y(i)+this%y(i-1))*this%dt
        enddo
    end function integrateTimeSeries
    !--------------------------------------------------------
    !> \brief Length of a time series
    !
    real function lengthTimeSeries(this)
        type (time_series) :: this
        lengthTimeSeries = (this%nsamp-1)*this%dt
    end function lengthTimeSeries
    !--------------------------------------------------------
    !> \brief  linearly combine two time series
    !
    function linearCombineTimeSeries(this,other,a,b) result(lcts)
        type (time_series) :: this,other,lcts
        real :: a,b
        if (abs((other%dt-this%dt)/this%dt) > tiny(this%dt)) then
            print *,'dt of time series differ!'
            stop
        endif
        if (this%nsamp /= other%nsamp) then
            print *,'nsamp of time series differ!'
            stop
        endif
        if (this%exdate) then
            call createZeroTimeSeries(lcts,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(lcts,this%nsamp,this%tanf,this%dt)
        endif
        lcts%y = a*this%y + b*other%y
    end function linearCombineTimeSeries
    !------------------------------------------------------------------------
    !> \brief Low pass Butterworth filtering of time series using recursive filters
    !> \param  fc corner frequency (Hz)
    !> \param  m order
    !>  use a lp1 and a lp2 filter to do the job
    !
    function lowPassButterworthRecursiveTimeSeries(this,fc,m) result(fts)
        type (time_series) :: this,fts,temp
        type (recursive_filter_coeff) :: rfclp1,rfclp2
        real :: fc,tc,h
        double precision :: pih
        integer :: m,mm,j
        !
        tc = 1./fc
        call copyTimeSeries(this,fts)
        mm = m/2
        if (m > 2*mm) then                                                  ! odd order filter
            call lowPassO1RecursiveFilterCoefficients(rfclp1,tc,this%dt)
            temp = recursiveFilterTimeSeries(fts,rfclp1)
            call dealloc(fts)
            fts = temp
        endif
        !
        pih=0.5d0*pi
        do j=1,mm
            h = real(sin(pih*(2*j-1)/m),4)
            call lowPassO2RecursiveFilterCoefficients(rfclp2,tc,h,this%dt)
            temp = recursiveFilterTimeSeries(fts,rfclp2)
            call dealloc(fts)
            fts = temp
        enddo
    end function lowPassButterworthRecursiveTimeSeries
    !-------------------------------------------------------------------------
    !> \brief  maximum of a time series
    !> \param this Time series object
    !
    real function maxTimeSeries(this)
        type (time_series) :: this
        maxTimeSeries = maxval(this%y(1:this%nsamp))
    end function maxTimeSeries
    !-------------------------------------------------------------------------
    !> \brief  minimum of a time series
    !> \param this Time series object
    !
    real function minTimeSeries(this)
        type (time_series) :: this
        minTimeSeries = minval(this%y(1:this%nsamp))
    end function minTimeSeries
    !-------------------------------------------------------------------------
    !> \brief  location of maximum of a time series
    !> \param this Time series object
    !
    integer function maxLocationTimeSeries(this)
        type (time_series) :: this
        maxLocationTimeSeries = maxloc(this%y(1:this%nsamp),1)
    end function maxLocationTimeSeries
    !-------------------------------------------------------------------------
    !> \brief misift between corresponding parts of time series
    !>  it is assumed that sigma fits to this
    !!  and that dt is the same and common sample times are identical
    !
    real function misfitTimeSeries(this,s,sigma)
        type (time_series), intent(in) :: this,s
        real, dimension(:) :: sigma
        double precision :: tend,tanf
        integer :: na1,na2,ne1,ne2
        !
        !  find common part of time series
        !
        tend = min(.tend.this,.tend.s)
        tanf = max(this%tanf,s%tanf)
        na1 = nint((tanf-this%tanf)/this%dt)+1
        ne1 = nint((tend-this%tanf)/this%dt)+1
        na2 = nint((tanf-s%tanf)/this%dt)+1
        ne2 = nint((tend-s%tanf)/this%dt)+1

        misfitTimeSeries = sum(((this%y(na1:ne1)-s%y(na2:ne2))/sigma(na1:ne1))**2)
    end function misfitTimeSeries
    !--------------------------------------------------------------------------------
    !> \brief Calclulate n-th order moving moment of a time series
    !! To avoid differences of powers of data values (because of loss of accuracy)
    !! Use b^4-a^4 = (b-a)*(b+a)*(b^2+a^2)
    !! Use b^3-a^3 = (b-a)*(b^2+ab+a^2)
    !
    function movingMomentTimeSeries(this,nord,nwin) result(res)
        type (time_series) :: this
        type (time_series) :: res
        integer :: nord,nwin,ns,j
        double precision, dimension(:), allocatable :: var
        double precision, dimension(:), allocatable :: s
        double precision, dimension(:), allocatable :: y
        !
        if (this%exdate) then
            call createZeroTimeSeries(res,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(res,this%nsamp,this%tanf,this%dt)
        endif
        !
        ns = this%nsamp
        allocate(var(ns),y(ns),s(ns))
        y = dble(this%y)
        if (nord > 1) then
            do j = 1,nwin-1                                                ! sum up j available samples
                var(j) = sum(y(1:j)**2)/j
            enddo
            var(nwin) = sum(y(1:nwin)**2)
            do j = nwin+1,ns                                               ! now recursively with nwin samples
                var(j) = var(j-1)+(y(j)+y(j-nwin))*(y(j)-y(j-nwin))
            enddo
            var(nwin:ns) = var(nwin:ns)/nwin
        endif
        if (nord == 2) then
            res%y = sngl(var)
        else if (nord == 3 .or. nord == 4) then
            do j = 1,nwin-1
                s(j) = sum(y(1:j)**nord)/j/var(j)**(nord/2.)
            enddo
            s(nwin) = sum(y(1:nwin)**nord)
            if (nord == 3) then
                do j = nwin+1,ns
                    s(j) = s(j-1)+(y(j)-y(j-nwin))*(y(j)**2+y(j)*y(j-nwin)+y(j-nwin)**2)
                enddo
            else if (nord == 4) then
                do j = nwin+1,ns
                    s(j) = s(j-1)+(y(j)-y(j-nwin))*(y(j)**2+y(j-nwin)**2)*(y(j)+y(j-nwin))
                enddo
            endif
            s(nwin:ns) = s(nwin:ns)/nwin/var(nwin:ns)**(nord/2.)
            res%y = sngl(s)
        else if (nord == 1) then                                                      ! moving average of nwin samples
            do j = 1,nwin-1
                s(j) = sum(y(1:j))/j
            enddo
            s(nwin) = sum(y(1:nwin))
            do j = nwin+1,ns
                s(j) = s(j-1)-y(j-nwin)+y(j)
            enddo
            s(nwin:ns) = s(nwin:ns)/nwin
            res%y = real(s,4)
        endif
        deallocate(var,s,y)
    end function movingMomentTimeSeries
    !-------------------------------------------------------
    !> \brief  get nsamp of time series
    !
    integer function nsampTimeSeries(this)
        type (time_series), intent(in) :: this
        nsampTimeSeries = this%nsamp
    end function nsampTimeSeries
    !--------------------------------------------------------
    !> \brief  Normalize time series to absolute maximum
    !
    subroutine normalizeTimeSeries(this)
        type (time_series) :: this
        real :: rmax
        rmax = maxval(abs(this%y))
        this%y = this%y/rmax
    end subroutine normalizeTimeSeries
    !-------------------------------------------------------------------------
    !> \brief plane misift between time series, no weighting, divide by nsamp
    !
    real function plainMisfitTimeSeries(this,s)
        type (time_series), intent(in) :: this,s
        double precision :: tend,tanf
        integer :: na1,na2,ne1,ne2
        !
        !  find common part of time series
        !
        tend = min(.tend.this,.tend.s)
        tanf = max(this%tanf,s%tanf)
        na1 = nint((tanf-this%tanf)/this%dt)+1
        ne1 = nint((tend-this%tanf)/this%dt)+1
        na2 = nint((tanf-s%tanf)/this%dt)+1
        ne2 = nint((tend-s%tanf)/this%dt)+1
        plainMisfitTimeSeries = sum(((this%y(na1:ne1)-s%y(na2:ne2)))**2)/(ne1-na1+1)
    end function plainMisfitTimeSeries
    !--------------------------------------------------------
    !> \brief  multiply two time series and optionally weight by sigma
    !
    function productTimeSeries(this,other,sig) result(prts)
        type (time_series) :: this,other,prts
        real, dimension(:), optional :: sig
        real, dimension(:), allocatable :: sigma
        !
        if (abs((other%dt-this%dt)/this%dt) > tiny(this%dt)) then
            print *,'dt of time series differ!'
            stop
        endif
        if (this%nsamp /= other%nsamp) then
            print *,'nsamp of time series differ!'
            stop
        endif
        allocate(sigma(this%nsamp))
        if (present(sig)) then; sigma = sig; else; sigma = 1.; endif
        if (this%exdate) then
            call createZeroTimeSeries(prts,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(prts,this%nsamp,this%tanf,this%dt)
        endif
        prts%y = this%y/sigma*other%y/sigma
        deallocate(sigma)
    end function productTimeSeries
    !--------------------------------------------------------------------
    !> \brief apply set of recursive filter coefficients to time series
    !
    function recursiveFilterTimeSeries(this,rfc) result(fts)
        type (time_series) :: this,fts
        type (recursive_filter_coeff) :: rfc
        integer :: n,j
        double precision :: xn,xa,xaa,ya,yaa,yn,xoff
        !
        n = this%nsamp
        if (this%exdate) then
            call createZeroTimeSeries(fts,n,this%stime,this%dt)
        else
            call createZeroTimeSeries(fts,n,this%tanf,this%dt)
        endif
        !
        if (rfc%typ == 'hp1' .or. rfc%typ == 'hp2') then
            xoff = this%y(1)                             ! first sample as baseline
        else
            xoff = 0.d0
        endif
        xa=0.d0; xaa=0.d0; ya=0.d0; yaa=0.d0
        do j=1,n
            xn = this%y(j)-xoff
            yn = rfc%f0*xn+rfc%f1*xa+rfc%f2*xaa+rfc%g1*ya+rfc%g2*yaa
            xaa=xa; xa=xn; yaa=ya; ya=yn;
            fts%y(j) = real(yn,4)
        enddo
    end function recursiveFilterTimeSeries
    !--------------------------------------------------------------------------------
    !> \brief  resample time series and do anti alias filtering
    !>
    !> \param  h filter coefficients
    !> \param  ndec	decimation factor
    !> \param  rts resampled time series
    !
    function resampleTimeSeries(this,h,ndec) result(rts)
        type (time_series) :: this,rts
        integer :: ndec
        real, dimension(:) :: h
        integer :: nh,i,j,k,skip,nf,nd
        !
        !  allocate f
        !
        nd = this%nsamp
        nf = (nd-1)/ndec+1
        if (this%exdate) then
            call createZeroTimeSeries(rts,nf,this%stime,this%dt*ndec)
        else
            call createZeroTimeSeries(rts,nf,this%tanf,this%dt*ndec)
        endif
        !
        !  perform the convolution, sum from k+1-nh to min(k,nd) to ensure
        !  that a) k-j+1 >= 1 and b) k-j+1 <= nh and c) j <= nd
        !  throw away the first (nh-1)/2 samples to compensate for phase shift of filter
        !  also throw away the last (nh-1)/2 samples at the end
        !  only calculate every ndec'th sample
        !
        nh=size(h)
        skip = (nh-1)/2
        i=0
        do k=skip+1,nd+skip,ndec
            i=i+1
            rts%y(i)=0.
            do j=max(k+1-nh,1),min(k,nd)
                rts%y(i) = rts%y(i)+this%y(j)*h(k-j+1)
            enddo
        enddo
    end function resampleTimeSeries
    !------------------------------------------------------------------------
    !> \brief  remove average of time series
    !
    function removeAverageTimeSeries(this) result(fts)
        type (time_series) :: this,fts
        real :: avg
        !
        if (this%exdate) then
            call createZeroTimeSeries(fts,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(fts,this%nsamp,this%tanf,this%dt)
        endif
        avg = averageTimeSeries(this)
        fts%y = this%y-avg
    end function removeAverageTimeSeries
    !-------------------------------------------------------------------------
    !> \brief  rms of a time series
    !
    real function rmsTimeSeries(this)
        type (time_series) :: this
        rmsTimeSeries = sqrt(sum(this%y(1:this%nsamp)**2)/this%nsamp)
    end function rmsTimeSeries
    !-------------------------------------------------------------------------
    !> \brief Normalize time series by setting the mean of the n largest amplitudes to ampnorm (usually 1)
    !
    function robustNormalizeTimeSeries(this,n,ampnorm) result(fts)
        type (time_series) :: this,fts
        integer :: n
        real, dimension(:),allocatable :: ampmax
        integer, dimension(:), allocatable :: imax
        real :: ampnorm,val,avgmax
        integer :: i,k
        !
        if (this%exdate) then
            call createZeroTimeSeries(fts,this%nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(fts,this%nsamp,this%tanf,this%dt)
        endif
        allocate(ampmax(n),imax(n))
        ampmax(1:n) = 0.
        do i = 1,this%nsamp
            val = abs(this%y(i))
            if (val > ampmax(1)) then; ampmax(1) = val; imax(1) = i; cycle; endif
            do k=2,n
                if (val > ampmax(k) .and. val <= ampmax(k-1)) then; ampmax(k) = val; imax(k) = i; exit; endif
            enddo
        enddo
        !	do k = 1,n
        !		print *,k,imax(k),ampmax(k)
        !	enddo
        avgmax = sum(ampmax)/n
        fts%y = this%y/avgmax*ampnorm
        deallocate(ampmax,imax)
    end function robustNormalizeTimeSeries
    !--------------------------------------------------------
    !> \brief  multiply time series by a factor
    !
    subroutine scalarMultiplyTimeSeries(this,factor)
        type (time_series) :: this
        real :: factor
        this%y = factor*this%y
    end subroutine scalarMultiplyTimeSeries
    !-------------------------------------------------------
    !> \brief  get tanf of time series in single precision (compatibility)
    !
    real function startTimeSeries(this)
        type (time_series), intent(in) :: this
        startTimeSeries = sngl(this%tanf)
    end function startTimeSeries
    !-------------------------------------------------------
    !> \brief  get tanf of time series in single precision
    !
    double precision function startDPTimeSeries(this)
        type (time_series), intent(in) :: this
        startDPTimeSeries = this%tanf
    end function startDPTimeSeries
    !--------------------------------------------------------
    !> \brief  starting time in full seconds after midnight
    !
    function startFullSecondsTimeSeries(this) result(tfs)
        type (time_series), intent(in) :: this
        integer :: tfs
        tfs = int(this%tanf)
    end function startFullSecondsTimeSeries
    !--------------------------------------------------------
    !> \brief  nanosecond part of starting time after midnight
    !
    function startNanoSecondsTimeSeries(this) result(tns)
        type (time_series), intent(in) :: this
        integer :: tns
        tns = nint((this%tanf-int(this%tanf))*1.d9)
    end function startNanoSecondsTimeSeries
    !--------------------------------------------------------
    !> \brief  set link flag to true
    !
    subroutine setLinkTrueTimeSeries(this)
        type (time_series) :: this
        this%link = .true.
    end subroutine setLinkTrueTimeSeries
    !--------------------------------------------------------
    !> \brief  set link flag to false
    !
    subroutine setLinkFalseTimeSeries(this)
        type (time_series) :: this
        this%link = .false.
    end subroutine setLinkFalseTimeSeries
    !--------------------------------------------------------
    !> \brief  shallow copy of time series
    !
    subroutine shallowCopyTimeSeries(this,copy)
        type (time_series), intent(in) :: this
        type (time_series), intent(out) :: copy
        copy%nsamp = this%nsamp
        copy%dt = this%dt
        copy%tanf = this%tanf
        copy%y => this%y
        copy%link = .true.
        copy%stime = this%stime
    end subroutine shallowCopyTimeSeries
    !-------------------------------------------------------------------------
    !> \brief Shift starting time of time series by given number of seconds
    !
    subroutine shiftTimeSeries(this,tshift)
        type (time_series) :: this
        double precision :: tshift
        !
        this%tanf = this%tanf+tshift
        if (tshift >= 0.d0) then
            this%stime = (this%stime).plus.tshift
        else
            this%stime = (this%stime).minus.(abs(tshift))
        endif
    end subroutine shiftTimeSeries
    !------------------------------------------------------------------------
    !> \brief Calculate signal to noise ratio of time series
    !> \param tn start of noise window in sec after start of time series
    !> \param wn length of noise window in seconds
    !> \param ts start of signal window in sec after start of time series
    !> \param ws length of signal window in seconds
    !
    function snrTimeSeries(this,tn,wn,ts,ws,ierr) result(res)
        type (time_series) :: this
        real :: tn,wn,ts,ws,res,rms_s,rms_n
        integer :: ierr,jn1,jn2,js1,js2
        !
        res = -1.
        jn1 = max(ceiling(tn/this%dt),1)
        jn2 = min(ceiling((tn+wn)/this%dt),this%nsamp)
        if (jn2 <= jn1) then; ierr = 1; return; endif
        !
        js1 = max(ceiling(ts/this%dt),1)
        js2 = min(ceiling((ts+ws)/this%dt),this%nsamp)
        if (js2 <= js1) then; ierr = 2; return; endif
        !
        rms_s = sqrt(sum(this%y(js1:js2)**2)/(js2-js1+1))
        rms_n = sqrt(sum(this%y(jn1:jn2)**2)/(jn2-jn1+1))
        res = rms_s/rms_n
        ierr = 0
    end function snrTimeSeries
    !--------------------------------------------------------
    !> \brief  get a pointer to the time series data
    !
    function traceTimeSeries(this) result(p)
        type (time_series), intent(in) :: this
        real, dimension(:), pointer :: p
        p => this%y
    end function traceTimeSeries
    !-----------------------------------------------------------------------------
    !> \brief  truncate time series to end time tend (up to closest sample after tend)
    !
    function truncateTimeSeries(this,tend) result(tts)
        type (time_series) :: this,tts
        double precision :: tend
        integer :: nsamp
        !
        nsamp = int((tend-this%tanf)/this%dt)+2     ! tend = tanf+(nsamp-1)*dt
        if (nsamp < 1) then
            print *,'truncateTimeSeries: tend is before tanf'
            stop
        endif
        if (this%exdate) then
            call createZeroTimeSeries(tts,nsamp,this%stime,this%dt)
        else
            call createZeroTimeSeries(tts,nsamp,this%tanf,this%dt)
        endif
        tts%y = this%y(1:nsamp)
    end function truncateTimeSeries
    !--------------------------------------------------------------------------------------
    !> \brief Nullify pointers to  time_series data
    !
    subroutine unlinkDataTimeSeries(this)
        type (time_series) :: this
        if (associated(this%y)) nullify(this%y)
    end subroutine unlinkDataTimeSeries
    !--------------------------------------------------------------------------------------
    !> \brief Reverse time axis of time_series
    !
    function ReverseTimeSeries(this) result(temp)
        type (time_series) :: this, temp
        integer :: i
        call createDPFromDataTimeSeries(temp,this%nsamp,this%tanf,this%dt,this%y)
        do i=1,this%nsamp
            temp%y(this%nsamp-i+1)=this%y(i)
        enddo
    end function ReverseTimeSeries
!
end module timeSeries
