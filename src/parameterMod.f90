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
module parameterMod

    ! module to read a parameter file to set up the simulation
    use constantsMod
    use fileParameterMod
    use errorMessage
    use slipInterfaceMod

    implicit none

    type :: parameterVar
        !private
        sequence
        character(len=80) :: title              !Title of the simulation
        logical :: log = .true.                 !Write information onto screen
        logical :: debug                        !Parameter to enable certain output for debugging purposes

        !external model (obsolete?)
        !logical :: external                    !read external model?
        character(len=255) :: externalfilename  !File containing the external model
        logical :: extvel                       !read external model?

        integer :: nproc                        !Number of processors for the calculation

        !Materials
        integer :: matn                         !Number of different materials

        !attenuation/viscoelasticity
        logical :: attenuation                  !Viscoelasticity used?
        real(kind=CUSTOM_REAL) :: f0_att        !frequency where the elastic parameters are equal to the anelastic ones
        real(kind=CUSTOM_REAL) :: f_max_att     !Maximum of anelastic frequency band
        real(kind=CUSTOM_REAL) :: fac_att       !Factor to describe lower end frequency band

        !Timeintegration
        integer :: timeint                      !Type of timeintegration (1:euler 2:rk2 3:rk3)
        integer :: nt                           !Number of timesteps
        logical :: autodt                       !Calculate timestep automatically?
        real(kind=CUSTOM_REAL) :: dt            !timestep (needs to be set if autodt = .false.)
        real(kind=CUSTOM_REAL) :: cfl           !Courantnumber (cfl value for dt)

        !Absorbing boundaries / PML parameters
        logical :: set_pml                      !pml
        real(kind=CUSTOM_REAL) :: pml_delta !pml
        real(kind=CUSTOM_REAL) :: pml_rc, pml_kmax, pml_afac !pml
        logical :: use_trigger                  !use sta_lta trigger for energy monitoring
        integer :: avg_window1, avg_window2     !lta and sta windows
        real(kind=CUSTOM_REAL) :: sta_lta_trigger !threshold for trigger

        !receiver
        real(kind=CUSTOM_REAL) :: rec_angle     !rotate receivers

        !Seismograms
        logical :: div                          ! Enables the seperate calculation of the radial component of the seismogram
        logical :: curl                         ! Enables the seperate calculation of the tangential component of the seismogram
        logical :: autoshift                    ! Enables the automatic shift of the seismogram by 1.2/f0
        real(kind=CUSTOM_REAL) :: plott0        ! Offset for the timeaxis for the seismogram

        !Flux
        integer :: fluxtype         !Specifies in which way the flux is calculated.
    end type parameterVar


    type movie_parameter
        !A series of switches to select which files are created for plotting
        integer :: frame        !Number of timesteps for each frame of the movie
        logical :: movie        !Select if any moviefiles are created
        logical :: trimesh      !Create files with average in each element
        logical :: points       !Create files with data for each point
        logical :: displacement !Plot displacement field
        logical :: velocity     !Plot velocity field
        logical :: stress       !Plot stress field
    end type

    contains

    subroutine readParfile(this, lsipar, movie, myrank, errmsg)
        implicit none
        type(error_message) :: errmsg
        type (parameterVar) :: this
        type(lsi_parameter) :: lsipar
        type(movie_parameter) :: movie
        integer, intent(in) :: myrank

        !local variables
        character(len=80) :: filename
        character(len=11) :: myname = "readParfile"

        integer :: ier
        logical :: log = .true.

        call addTrace(errmsg,myname)

        filename=trim('data/parfile')
        open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
        if (ier /= 0) then
            call add(errmsg,2,'Could not open file.' ,myname, filename)
            call print(errmsg)
            stop
        endif

        ! Cycle to read the parameters form the parameter file. The order of appearence of the parameters is not important
        call readStringPar(this%title, "title", filename, 0, errmsg)
        call readStringPar(this%externalfilename, "externalfilename", filename, 0, errmsg)
        call readIntPar(this%nproc, "nproc", filename, 0, errmsg)
        call readIntPar(this%timeint, "timeint", filename, 0, errmsg)
        call readIntPar(movie%frame, "frame", filename, 0, errmsg)
        call readIntPar(this%nt, "nt", filename, 0, errmsg)
        call readIntPar(this%avg_window1, "avg_window1", filename, 0, errmsg)
        call readIntPar(this%avg_window2, "avg_window2", filename, 0, errmsg)
        call readIntPar(this%fluxtype, "fluxtype", filename, 0, errmsg)
        call readFloatPar(this%cfl, "cfl", filename, 0, errmsg)
        call readFloatPar(this%pml_delta, "pml_delta", filename, 0, errmsg)
        call readFloatPar(this%pml_rc,"pml_rc", filename, 0, errmsg)
        call readFloatPar(this%pml_kmax, "pml_kmax", filename, 0, errmsg)
        call readFloatPar(this%pml_afac,"pml_afac", filename, 0, errmsg)
        call readFloatPar(this%sta_lta_trigger,"sta_lta_trigger", filename, 0, errmsg)
        call readFloatPar(this%f_max_att,"f_max_att", filename, 0, errmsg)
        call readFloatPar(this%f0_att,"f0_att", filename, 0, errmsg)
        call readFloatPar(this%fac_att,"att_factor", filename, 0, errmsg)
        call readFloatPar(this%rec_angle,"rec_angle", filename, 0, errmsg)
        call readFloatPar(this%dt,"dt", filename, 0, errmsg)
        call readLogicalPar(this%log, "log", filename, 0, errmsg)
        call readLogicalPar(this%extvel, "extvel", filename, 0, errmsg)
        call readLogicalPar(this%set_pml, "set_pml", filename, 0, errmsg)
        call readLogicalPar(this%attenuation, "attenuation", filename, 0, errmsg)
        call readLogicalPar(this%autodt, "autodt", filename, 0, errmsg)
        call readLogicalPar(this%debug, "debug", filename, 0, errmsg)
        call readLogicalPar(this%use_trigger, "use_trigger", filename, 0, errmsg)
        call readLogicalPar(lsipar%lsi, "lsi", filename, 0, errmsg)
        call readLogicalPar(lsipar%normal, "normal", filename, 0, errmsg)
        call readLogicalPar(lsipar%tangential, "tangential", filename, 0, errmsg)
        call readLogicalPar(movie%movie, "movie", filename, 0, errmsg)
        call readLogicalPar(movie%displacement, "save_movie_displacement", filename, 0, errmsg)
        call readLogicalPar(movie%velocity, "save_movie_velocity", filename, 0, errmsg)
        call readLogicalPar(movie%stress, "save_movie_stress", filename, 0, errmsg)
        call readLogicalPar(movie%trimesh, "save_movie_trimesh", filename, 0, errmsg)
        call readLogicalPar(movie%points, "save_movie_points", filename, 0, errmsg)
        call readLogicalPar(this%div, "div", filename, 0, errmsg)
        call readLogicalPar(this%curl, "curl", filename, 0, errmsg)
        call readLogicalPar(this%autoshift, "autoshift", filename, 0, errmsg)
        if (.not. this%autoshift) then
            call readFloatPar(this%plott0,"plott0", filename, 0, errmsg)
        end if
        close(19)

        !test if certain parameters have been read/entered correctly
        call ErrorMessages(this, lsipar, movie, myname, errmsg, filename)


        log = this%log

        if (log .and. myrank == 0) then
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, a10, a30)")   "|              Title of the simulation: ", this%title, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                               Basic parameters                               |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                           Create log: ", this%log, "                             |"
            write (*,"(a40, l10, a30)")   "|                           Debug mode: ", this%debug, "                             |"
            write (*,"(a40, i10, a30)")   "|                 Number of processors: ", this%nproc, "                             |"
            write (*,"(a40, l10, a30)")   "|        Load external velocitiy model: ", this%extvel, "                             |"
            if (this%extvel) write (*,*) trim(this%externalfilename)
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                Movie parameters                              |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                           Save movie: ", movie%movie, "                             |"
            if (movie%movie) then
                write (*,"(a40, i10, a30)")   "|   Number of timesteps for movieframe: ", movie%frame, "                             |"
                write (*,"(a40, l10, a30)")   "|           Create TriMesh-movie-files: ", movie%trimesh, "                             |"
                write (*,"(a40, l10, a30)")   "|            Create Points-movie-files: ", movie%points, "                             |"
                write (*,"(a40, l10, a30)")   "|                  Plot velocity field: ", movie%velocity, "                             |"
                write (*,"(a40, l10, a30)")   "|              Plot displacement field: ", movie%displacement, "                             |"
                write (*,"(a40, l10, a30)")   "|                    Plot stress field: ", movie%stress, "                             |"
            end if
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                         Parameters for timeintegration                       |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, i10, a30)")   "|                      Timeintegartion: ", this%timeint, "                             |"
            write (*,"(a40, i10, a30)")   "|                  Number of timesteps: ", this%nt, "                             |"
            write (*,"(a40, l10, a30)")   "|                               autodt: ", this%autodt, "                             |"
            if (.not. this%autodt) then
                write (*,"(a40, e10.3, a30)") "|                                   dt: ", this%dt, "                             |"
            endif
            write (*,"(a40, f10.1, a30)") "|                            cfl value: ", this%cfl, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                 PML parameters                               |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                          pml enabled: ", this%set_pml, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                        pml thickness: ", this%pml_delta, "                             |"
            write (*,"(a40, f10.3, a30)") "|           pml reflection coefficient: ", this%pml_rc, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                             pml kmax: ", this%pml_kmax, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                      factor for amax: ", this%pml_afac, "                             |"
            write (*,"(a40, i10, a30)")   "|                           lta window: ", this%avg_window1, "                             |"
            write (*,"(a40, i10, a30)")   "|                           sta window: ", this%avg_window2, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                    sta_lta threshold: ", this%sta_lta_trigger, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                          Parameters for attenuation                          |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   &
              "|                          Attenuation: ", this%attenuation, "                             |"
            if (this%attenuation) then
                write (*,"(a40, f10.1, a30)") "|                               f0_att: ", this%f0_att, "                             |"
                write (*,"(a40, f10.1, a30)") "|                            f_max_att: ", this%f_max_att, "                             |"
                write (*,"(a40, f10.1, a30)") "|                              fac_att: ", this%fac_att, "                             |"
            end if
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                          Parameters for receivers                            |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, f10.1, a30)") "|                            rec_angle: ", this%rec_angle, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                       Parameters regarding seismograms                       |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|  Auto-shift the seismogram by 1.2/f0: ", this%autoshift, "                             |"
            if (.not. this%autoshift) then
                write (*,"(a40, f10.7, a30)") "|                               Offset: ", this%plott0, "                             |"
            end if
            write (*,"(a40, l10, a30)")   "|                 Calculate divergence: ", this%div, "                             |"
            write (*,"(a40, l10, a30)")   "|                       Calculate curl: ", this%curl, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        end if
    end subroutine readParfile

    subroutine ErrorMessages(par, lsipar, movie, myname, errmsg, filename)
        ! This subroutine produces error messages if certain parameters are set in the wrong way.

        !input:
        type(lsi_parameter) :: lsipar
        type(parameterVar) :: par
        type(error_message) :: errmsg
        type(movie_parameter) :: movie
        character(len=*) :: myname, filename


        if (movie%movie) then
            if (movie%frame < 0) call add(errmsg, 2, "Number of time steps per movie frame needs to be positive.", myname, filename)
            if (movie%frame > par%nt) call add(errmsg, 2, "Number of time steps per movie frame less than the total number of time steps.", myname, filename)
            if (movie%frame == 0) call add(errmsg, 1, "If the number of time steps per movie frame is 0 no movie-files are created.", myname, filename)
        end if
        if (par%nproc <= 0) call add(errmsg, 2, "Number of cores needs to be positive.", myname, filename)
        if (par%nproc == 1) call add(errmsg, 1, "Single core calculations are not advised.", myname, filename)
        if (par%fac_att <= 0) call add(errmsg, 2, "Parameter 'fac_att' must be positive. Check the documentation for recommended values", myname, filename)
        if (par%f0_att <= 0) call add(errmsg, 2, "Parameter 'f0_att' must be positive. Check the documentation for recommended values", myname, filename)
        if (par%cfl <= 0.0) call add(errmsg, 2, "Parameter 'cfl' has to be positive. Recommended range is 0.0 < cfl <= 1.0", myname, filename)
        if (par%cfl > 1.0) call add(errmsg, 1, "Parameter 'cfl' might be too high. Recommended range is 0.0 < cfl <=0 1.0", myname, filename)
        if (.not. par%autodt .and. par%dt < epsilon(par%dt)) call add(errmsg, 2, "Choose dt /= 0, or autodt = .true.!", myname, filename )
        if( par%timeint < 1 .or. par%timeint > 4) call add(errmsg, 2, "Parameter to select time stepping variant is out of range. Select 1, 2, 3 or 4.", myname, filename)
        if (par%fluxtype == 0 .and. par%attenuation) call add(errmsg, 2, "Currenly the precalculated Riemannfluxes do not support attenuation!", myname, filename )
        if (lsipar%lsi) then
            if (par%fluxtype /= 0) then
                call add(errmsg, 2, "Parameter 'fluxtype' must be '0' for calculations including fractures.", myname, filename)
                return
            end if
            if (par%timeint == 3) then
                call add(errmsg, 2, "Calculation for the fracture influence is not available for tcd-rk3 integration!", myname, filename)
                return
            end if
            if (par%attenuation) then
                call add(errmsg, 2, "Calculation for the fracture influence is not compatible with inversion!", myname, filename)
                return
            end if
        end if
    end subroutine
end module parameterMod
