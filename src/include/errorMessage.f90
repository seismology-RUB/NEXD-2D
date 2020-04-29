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

!----------------------------------------------------------------------
!> \brief  Module to handle errors occuring in functions or subroutines
!> \par Description
!> Functions producing exceptions may return an error_message object to
!! specify the name of the function where the exception happened, to
!! provide a description of the error and to rate its severity
!! ranging from OK (no error), warning to error.
!<---------------------------------------------------------------------
 module errorMessage
      use iso_fortran_env
    use realloc
    use constantsMod
    implicit none
    interface new
        module procedure createErrorMessage
        module procedure createOKErrorMessage
    end interface
    interface add; module procedure updateErrorMessage; end interface
    interface dealloc; module procedure deallocErrorMessage; end interface
    interface print; module procedure printErrorMessage; end interface
    interface addTrace; module procedure addTraceErrorMessage; end interface
    interface operator (.level.); module procedure getLevelErrorMessage; end interface
    type error_message
        private
        integer :: level = 0                 !< error level: 0 = success, 1 = warning, 2 = error
        character (len=400), dimension(:), pointer :: messages => null()    !< error/warning messages
        character (len=132) :: fctname = ''    !< function name where error message was created
        character (len=132), dimension(:), pointer :: trace => null()  !< function names through which error was propagated
        character (len=132) :: filename = ''
    end type
!
 contains
!-----------------------------------------------------------------
!> \brief Create an error message
!
    subroutine createErrorMessage(this,level,message,fctname)
    type (error_message) :: this
    integer :: level
    character (len=*) :: message,fctname
    call deallocErrorMessage(this)
    this%level = level
    call addMessageErrorMessage(this,level,message,fctname)
    this%fctname = fctname
    end subroutine createErrorMessage
!------------------------------------------------------------------
!> \brief Create a default OK message
!
    subroutine createOKErrorMessage(this,fctname)
    type (error_message) :: this
    character (len=*) :: fctname
    call deallocErrorMessage(this)
    this%level = 0
    this%fctname = fctname
    end subroutine createOKErrorMessage
!-----------------------------------------------------------------
!> \brief Update an error message
!
    subroutine updateErrorMessage(this,level,message,fctname, filename)
    type (error_message) :: this
    integer :: level
    character (len=*) :: message,fctname
    character (len=*), optional :: filename
    ! error levels can only increase! (from success to warning or error, from warning to error)
    if(level>this%level) this%level = level
    if(this%fctname == '') then
        ! error message was actually not created yet
        this%fctname = fctname
    endif
    if (PRESENT(filename)) then
       this%filename = filename
    end if
    call addMessageErrorMessage(this,level,message,fctname)
    end subroutine updateErrorMessage
!-----------------------------------------------------------------
!> \brief Deallocate error message
!
    subroutine deallocErrorMessage(this)
    type (error_message) :: this
    if (associated(this%trace)) deallocate(this%trace)
    if (associated(this%messages)) deallocate(this%messages)
    this%level = 0
    this%fctname = ''
    end subroutine deallocErrorMessage
!------------------------------------------------------------------
!> \brief Print error message
!
    subroutine printErrorMessage(this, displayTrace)
    type (error_message), intent(in) :: this
    logical, optional :: displayTrace
    logical, parameter :: displayTrace_default = .false.
    logical :: displayTrace_local
    integer :: j
    write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    select case (this%level)
    case (0); write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>   SUCCESS   >>>>>>>>>>>>>>>'
    case (1); write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>   WARNING   >>>>>>>>>>>>>>>'
    case (2); write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>>   ERROR   >>>>>>>>>>>>>>>>'
    end select

    if (PRESENT(displayTrace)) then
       displayTrace_local = displayTrace
    else
       displayTrace_local = displayTrace_default
    end if

    if (associated(this%messages)) then
        write(OUTPUT_UNIT,*) '>>>> '
        do j=1,size(this%messages)
            write(OUTPUT_UNIT,*) '>>>> ',trim(this%messages(j))
        enddo
    endif

    if (this%filename == ' ') then
        write(OUTPUT_UNIT,*) '>>>> '
        write(OUTPUT_UNIT,*) '>>>> message created in --> ',trim(this%fctname)
    else
        write(OUTPUT_UNIT,*) '>>>> '
        write(OUTPUT_UNIT,*) '>>>> Affected file --> ',trim(this%filename)
    end if
    if (displayTrace_local) then
        if (associated(this%trace)) then
            do j=1,size(this%trace)
                write(OUTPUT_UNIT,*) '>>>>     passed through --> ',trim(this%trace(j))
            enddo
        endif
    end if
    write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    write(OUTPUT_UNIT,*) ''
    end subroutine printErrorMessage
!---------------------------------------------------------------
!> \brief Add a trace to error message
!
    subroutine addTraceErrorMessage(this,fctname)
    type (error_message) :: this
    character (len=*) :: fctname
    integer :: n
    if (associated(this%trace)) then
        n = size(this%trace)
        this%trace => reallocate(this%trace,n+1)
        this%trace(n+1) = trim(fctname)
    else
        allocate(this%trace(1))
        this%trace(1) = trim(fctname)
    endif
    end subroutine addTraceErrorMessage
!---------------------------------------------------------------
!> \brief Add a message to error message object
!
    subroutine addMessageErrorMessage(this,level,message,fctname)
    type (error_message) :: this
    integer :: level
    character (len=*) :: message,fctname
    integer :: n
    character (len=400) :: string
    select case (level)
    case(0); string = 'SUCCESS'
    case(1); string = 'WARNING'
    case(2); string = 'ERROR'
    end select

    string = trim(string)//' in '//trim(fctname)//' --> '//trim(message)
    if (associated(this%messages)) then
        n = size(this%messages)
        this%messages => reallocate(this%messages,n+1)
        this%messages(n+1) = trim(string)
    else
        allocate(this%messages(1))
        this%messages(1) = trim(string)
    endif
    end subroutine addMessageErrorMessage
!---------------------------------------------------------------
!> \brief Get level of error message
!
    integer function getLevelErrorMessage(this)
    type (error_message), intent(in) :: this
    getLevelErrorMessage = this%level
    end function getLevelErrorMessage

    subroutine pmlConflict(parameter, name, pml_minX, pml_minZ, pml_maxX, pml_maxZ, myname, filename, errmsg)
        !input
        !die subroutine muss noch erweitert werden. aktuell wird noch nicht ��berpr��pft ob die jeweilige seite wirklich eine pml ist.
        type(error_message) :: errmsg
        character(len=*) :: name, myname, filename
        real(kind = CUSTOM_REAL) :: parameter, pml_minX, pml_minZ, pml_maxX, pml_maxZ

        if (index(name, "x") > 0) then
            if (parameter < pml_minX) then
                call add(errmsg,2,"Parameter '"//trim(name)//"' is inside the left PML boundary!",myname, filename)
            end if
            if (parameter > pml_maxX) then
                call add(errmsg,2,"Parameter '"//trim(name)//"' is inside the right PML boundary!",myname, filename)
            end if
        else if (index(name, "z") > 0) then
            if (parameter < pml_minZ) then
                call add(errmsg,2,"Parameter '"//trim(name)//"' is inside the lower PML boundary!",myname, filename)
            end if
            if (parameter > pml_maxZ) then
                call add(errmsg,2,"Parameter '"//trim(name)//"' is inside the top PML boundary!",myname, filename)
            end if
        end if
    end subroutine

    subroutine modelConflict(parameter, name, minX, minZ, maxX, maxZ, myname, filename, errmsg)
        !input
        type(error_message) :: errmsg
        character(len=*) :: name, myname, filename
        real(kind = CUSTOM_REAL) :: parameter, minX, minZ, maxX, maxZ

        if (index(name, "x") > 0) then
            if (parameter < minX) then
                call add(errmsg,2,"Parameter '"//trim(name)//"' is outside the left boundary!",myname, filename)
            end if
            if (parameter > maxX) then
                call add(errmsg,2,"Parameter '"//trim(name)//"' is outside the right boundary!",myname, filename)
            end if
        else if (index(name, "z") > 0) then
            if (parameter < minZ) then
                call add(errmsg,2,"Parameter '"//trim(name)//"' is outside the lower boundary!",myname, filename)
            end if
            if (parameter > maxZ) then
                call add(errmsg,2,"Parameter '"//trim(name)//"' is outside the top boundary!",myname, filename)
            end if
        end if
    end subroutine
!
 end module errorMessage
