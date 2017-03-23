!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Marc S. Boxberg (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Thomas Möller (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2015-2017 Andre Lamert (Ruhr-Universitaet Bochum, Germany)
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
module parameterMod

    ! module to read a parameter file to set up the simulation
    use constantsMod
    use fileParameterMod
    use errorMessage

    implicit none

    type :: parameterVar
        !private
        sequence
        character(len=80) :: title              !Title of the simulation
        logical :: log = .true.                 !Write information onto screen
        logical :: debug                        !Parameter to enable certain output for debugging purposes
        
        !adjoint
        logical :: save_forward                 !save the forward wave field
        logical :: adjoint                      !compute adjoint kernels
        logical :: inv_model                    !use updated model?
        integer :: adj_nstep                    !nsteps for writing and reading the wavefield in adjoint simulations
        
        !external model (obsolete?)
        !logical :: external                    !read external model?
        character(len=255) :: externalfilename   !File containing the external model
        logical :: extvel                       !read external model?
        
        integer :: nproc                        !Number of processors for the calculation
        logical :: strongform                   !strong or weak form of DG
        
        !Materials
        integer :: matn                         !Number of different materials

        !attenuation/viscoelasticity
        logical :: attenuation                  !Viscoelasticity used?
        real(kind=CUSTOM_REAL) :: f0_att        !frequency for anelastic modulus def:500
        real(kind=CUSTOM_REAL) :: fr_att        !lower freq bound f0/3.4? def:144
        
        !Timeintegration
        integer :: timeint                      !Type of timeintegration (1:euler 2:rk2 3:rk3)
        integer :: nt                           !Number of timesteps
        logical :: autodt                       !Calculate timestep automatically?
        real(kind=CUSTOM_REAL) :: dt            !timestep (needs to be set if autodt = .false.)
        real(kind=CUSTOM_REAL) :: cfl           !Courantnumber (cfl value for dt)
        
        !Absorbing boundaries / PML parameters
        logical :: absorb_top, absorb_bottom, absorb_left, absorb_right !absorbing boundaries
        logical :: set_pml                      !pml
        real(kind=CUSTOM_REAL) :: xmin, xmax,zmin,zmax,pml_delta !pml
        real(kind=CUSTOM_REAL) :: pml_rc, pml_kmax, pml_afac !pml
        logical :: use_trigger                  !use sta_lta trigger for energy monitoring
        integer :: avg_window1, avg_window2     !lta and sta windows
        real(kind=CUSTOM_REAL) :: sta_lta_trigger !threshold for trigger
        
        !receiver
        real(kind=CUSTOM_REAL) :: rec_angle     !rotate receivers
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
        logical :: p1           !Plot pressure of the first fluid
        logical :: q1           !Plot velocity of the first fluid
        logical :: p2           !Plot pressure of the second fluid
        logical :: q2           !Plot velocity of the second fluid
    end type

    contains

    subroutine readParfile(this, movie, myrank, errmsg)
        implicit none
        type(error_message) :: errmsg
        type (parameterVar) :: this
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
            call add(errmsg,2,'Could not open file: '//trim(filename),myname)
            return
        endif

        ! Cycle to read the parameters form the parameter file. The order of appearence of the parameters is not important
        call readStringPar(this%title, "title", filename, 0)
        call readStringPar(this%externalfilename, "externalfilename", filename, 0)
        call readIntPar(this%nproc, "nproc", filename, 0)
        call readIntPar(this%timeint, "timeint", filename, 0)
        call readIntPar(movie%frame, "frame", filename, 0)
        call readIntPar(this%nt, "nt", filename, 0)
        call readIntPar(this%avg_window1, "avg_window1", filename, 0)
        call readIntPar(this%avg_window2, "avg_window2", filename, 0)
        call readIntPar(this%adj_nstep, "adj_nstep", filename, 0)
        call readFloatPar(this%cfl, "cfl", filename, 0)
        call readFloatPar(this%xmin, "xmin", filename, 0)
        call readFloatPar(this%xmax, "xmax", filename, 0)
        call readFloatPar(this%zmin, "zmin", filename, 0)
        call readFloatPar(this%zmax, "zmax", filename, 0)
        call readFloatPar(this%pml_delta, "pml_delta", filename, 0)
        call readFloatPar(this%pml_rc,"pml_rc", filename, 0)
        call readFloatPar(this%pml_kmax, "pml_kmax", filename, 0)
        call readFloatPar(this%pml_afac,"pml_afac", filename, 0)
        call readFloatPar(this%sta_lta_trigger,"sta_lta_trigger", filename, 0)
        call readFloatPar(this%f0_att,"f0_att", filename, 0)
        call readFloatPar(this%fr_att,"fr_att", filename, 0)
        call readFloatPar(this%rec_angle,"rec_angle", filename, 0)
        call readFloatPar(this%dt,"dt", filename, 0)
        call readLogicalPar(this%log, "log", filename, 0)
        call readLogicalPar(this%extvel, "extvel", filename, 0)
        call readLogicalPar(this%set_pml, "set_pml", filename, 0)
        call readLogicalPar(this%attenuation, "attenuation", filename, 0)
        call readLogicalPar(this%autodt, "autodt", filename, 0)
        call readLogicalPar(this%strongform, "strongform", filename, 0)
        call readLogicalPar(this%debug, "debug", filename, 0)
        call readLogicalPar(this%save_forward, "save_forward", filename, 0)
        call readLogicalPar(this%adjoint, "adjoint", filename, 0)
        call readLogicalPar(this%inv_model, "inv_model", filename, 0)
        call readLogicalPar(this%use_trigger, "use_trigger", filename, 0)
        call readLogicalPar(movie%movie, "movie", filename, 0)
        call readLogicalPar(movie%displacement, "save_movie_displacement", filename, 0)
        call readLogicalPar(movie%velocity, "save_movie_velocity", filename, 0)
        call readLogicalPar(movie%stress, "save_movie_stress", filename, 0)
        call readLogicalPar(movie%p1, "save_movie_p1", filename, 0)
        call readLogicalPar(movie%q1, "save_movie_q1", filename, 0)
        call readLogicalPar(movie%p2, "save_movie_p2", filename, 0)
        call readLogicalPar(movie%q2, "save_movie_q2", filename, 0)
        call readLogicalPar(movie%trimesh, "save_movie_trimesh", filename, 0)
        call readLogicalPar(movie%points, "save_movie_points", filename, 0)
        close(19)
      
        !test if certain parameters have been read/entered correctly
        if (.not. this%autodt .and. this%dt == 0.) then
            call add(errmsg, 2, "Choose dt /= 0, or autodt = .true.! Abort...", myname )
        endif

        log = this%log

        if (log .and. myrank == 0) then
            write (*, "(a80)") "--------------------------------------------------------------------------------"
            write (*,"(a40, a10, a30)")   "|              Title of the simulation: ", this%title, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                               Basic parameters                               |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                          Strong form: ", this%strongform, "                             |"
            write (*,"(a40, l10, a30)")   "|                           Create log: ", this%log, "                             |"
            write (*,"(a40, l10, a30)")   "|                           Debug mode: ", this%debug, "                             |"
            write (*,"(a40, i10, a30)")   "|                 Number of processors: ", this%nproc, "                             |"
            write (*,"(a40, l10, a30)")   "|        Load external velocitiy model: ", this%extvel, "                             |"
            if (this%extvel) write (*,*) trim(this%externalfilename)
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                              Adjoint Simulations                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                              adjoint: ", this%adjoint, "                             |"
            if (this%adjoint) then
                write (*,"(a40, l10, a30)")   "|                         save_forward: ", this%save_forward, "                             |"
                write (*,"(a40, l10, a30)")   "|                            inv_model: ", this%inv_model, "                             |"
                write (*,"(a40, l10, a30)")   "|                            adj_nstep: ", this%adj_nstep, "                             |"
            end if
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                Movie parameter                               |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                           Save movie: ", movie%movie, "                             |"
            if (movie%movie) then
                write (*,"(a40, i10, a30)")   "|   Number of timesteps for movieframe: ", movie%frame, "                             |"
                write (*,"(a40, l10, a30)")   "|           Create TriMesh-movie-files: ", movie%trimesh, "                             |"
                write (*,"(a40, l10, a30)")   "|            Create Points-movie-files: ", movie%points, "                             |"
                write (*,"(a40, l10, a30)")   "|                  Plot velocity field: ", movie%velocity, "                             |"
                write (*,"(a40, l10, a30)")   "|              Plot displacement field: ", movie%displacement, "                             |"
                write (*,"(a40, l10, a30)")   "|                    Plot stress field: ", movie%stress, "                             |"
                write (*,"(a40, l10, a30)")   "|    Plot pressure field for 1st fluid: ", movie%p1, "                             |"
                write (*,"(a40, l10, a30)")   "|    Plot velocity field for 1st fluid: ", movie%q1, "                             |"
                write (*,"(a40, l10, a30)")   "|    Plot pressure field for 2nd fluid: ", movie%p2, "                             |"
                write (*,"(a40, l10, a30)")   "|    Plot velocity field for 2nd fluid: ", movie%q2, "                             |"
            end if
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                         Parameter for timeintegration                        |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, i10, a30)")   "|                      Timeintegartion: ", this%timeint, "                             |"
            write (*,"(a40, i10, a30)")   "|                  Number of timesteps: ", this%nt, "                             |"
            write (*,"(a40, l10, a30)")   "|                               autodt: ", this%autodt, "                             |"
            write (*,"(a40, e10.4, a30)") "|                                   dt: ", this%dt, "                             |"
            write (*,"(a40, f10.1, a30)") "|                            cfl value: ", this%cfl, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                 PML parameter                                |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                          pml enabled: ", this%set_pml, "                             |"
            write (*,"(a40, f10.1, a30)") "|                                 xmin: ", this%xmin, "                             |"
            write (*,"(a40, f10.1, a30)") "|                                 xmax: ", this%xmax, "                             |"
            write (*,"(a40, f10.1, a30)") "|                                 zmin: ", this%zmin, "                             |"
            write (*,"(a40, f10.1, a30)") "|                                 zmax: ", this%zmax, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                        pml thickness: ", this%pml_delta, "                             |"
            write (*,"(a40, f10.3, a30)") "|           pml reclection coefficient: ", this%pml_rc, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                             pml kmax: ", this%pml_kmax, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                      factor for amax: ", this%pml_afac, "                             |"
            write (*,"(a40, i10, a30)")   "|                            lta window:", this%avg_window1, "                             |"
            write (*,"(a40, i10, a30)")   "|                            sta window:", this%avg_window2, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                    sta_lta threshold: ", this%sta_lta_trigger, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                          Parameter for attenuation                           |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   &
              "|                          Attenuation: ", this%attenuation, "                             |"
            if (this%attenuation) then
                write (*,"(a40, f10.1, a30)") "|                               f0_att: ", this%f0_att, "                             |"
                write (*,"(a40, f10.1, a30)") "|                               fr_att: ", this%fr_att, "                             |"
            end if
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                           Parameter for receiver                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, f10.1, a30)") "|                            rec_angle: ", this%rec_angle, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        end if
    end subroutine readParfile
end module parameterMod
