!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
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
module mpiMod
    !use mpi ! mpi mit echtem tab weil externe library
    use mpiincludeMod
    use constantsMod

    implicit none

    ! module to access the mpi routines
    contains

    ! call MPI_Init to start MPI
    subroutine init_mpi()
        implicit none
    
        integer :: ierr

        call MPI_Init(ierr)
        call check_ierr(ierr,"Error in Initializing MPI")
    end subroutine init_mpi

    ! call MPI_COMM_SIZE to get the world
    subroutine comm_size(size)
        implicit none
    
        integer :: size
        integer :: ierr
    
        call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
        call check_ierr(ierr,"Error in MPI_COMM_SIZE")
    end subroutine comm_size

    ! call MPI_COMM_RANK to geht the Rank of the node
    subroutine comm_rank(rank)
        implicit none
        integer :: ierr
        integer :: rank

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call check_ierr(ierr,"Error in MPI_COMM_RANK")
    end subroutine comm_rank
  
    ! call MPI_BARRIER to sync the processes
    subroutine sync_mpi()
        implicit none
        integer :: ierr
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_BARRIER")
    end subroutine sync_mpi

    ! call MPI_Isend to send the data (vektor)
    subroutine isendV_real(buf,count,dest,tag,req)
        implicit none

        integer :: count, dest,tag,req
        real(kind=CUSTOM_REAL), dimension(count) :: buf
        integer :: ierr
    
        call MPI_ISSEND(buf(1),count,MPI_REAL,dest,tag,MPI_COMM_WORLD,req,ierr)
        call check_ierr(ierr,"Error in MPI_ISSEND")
    end subroutine isendV_real

    ! call MPI_Irecv to receive the data (vektor)
    subroutine irecV_real(buf,count,dest,tag,req)
        implicit none

        integer :: count, dest,tag,req
        real(kind=CUSTOM_REAL), dimension(count) :: buf
        integer :: ierr

        call MPI_IRECV(buf(1),count,MPI_REAL,dest,tag,MPI_COMM_WORLD,req,ierr)
        call check_ierr(ierr,"Error in MPI_IREC")
    end subroutine irecV_real
    
    ! call mpi_wait to wait for the request
    subroutine wait_req(req)
        implicit none
        integer :: req
        integer, dimension(MPI_STATUS_SIZE) :: req_status
        integer :: ierr

        call mpi_wait(req,req_status,ierr)
        call check_ierr(ierr,"Error in MPI_WAIT")
    end subroutine wait_req

    ! call MPI_reduce to receive the maximum of a dataset
    subroutine maxval_real(send, rec)
        implicit none

        real :: send, rec
        integer :: ierr
    
        call MPI_REDUCE(send,rec,1,MPI_REAL, &
                      MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real")
    end subroutine maxval_real

    ! call MPI_broadcast to receive the maximum of a dataset
    subroutine maxval_real_all(send, rec)
        implicit none

        real :: send, rec
        integer :: ierr
    
        call MPI_ALLREDUCE(send,rec,1,MPI_REAL, &
                      MPI_MAX,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real")
    end subroutine maxval_real_all

    ! call MPI_reduce to receive the sum of a dataset
    subroutine sum_real(send, rec)
        implicit none

        real :: send, rec
        integer :: ierr

        call MPI_REDUCE(send,rec,1,MPI_REAL, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, sum_real")
    end subroutine sum_real

    ! call MPI_reduce to receive the sum of a dataset
    subroutine sum_all_real(send, rec)
        implicit none

        real :: send, rec
        integer :: ierr
    
        call MPI_ALLREDUCE(send,rec,1,MPI_REAL, &
                      MPI_SUM,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_ALLREDUCE, sum_all_real")
    end subroutine sum_all_real

    ! call MPI_reduce to receive the sum of a dataset
    subroutine sum_int(send, rec)
        implicit none

        integer :: send, rec
        integer :: ierr
    
        call MPI_REDUCE(send,rec,1,MPI_INTEGER, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, sum_int")
    end subroutine sum_int

    ! call MPI_finalize to stop MPI
    subroutine finalize_mpi()
        implicit none

        integer :: ierr
    
        call MPI_Finalize(ierr)
        if (ierr /= MPI_SUCCESS) then
           write(*,*) "Error in finalizing MPI"
        end if
    end subroutine finalize_mpi

    ! call MPI_about to stop all comunications
    subroutine stop_mpi()
        implicit none
        integer :: ierr
        integer :: errcode
    
        call MPI_Abort(MPI_COMM_WORLD, errcode,ierr)
        stop 'error, program stop in stop_mpi'
    end subroutine stop_mpi

    subroutine check_ierr(ierr,msg)
        implicit none
        integer :: ierr
        character(len=*) :: msg
    
        if (ierr /= MPI_SUCCESS) then
            write(*,*) ierr,trim(msg)
            call stop_mpi()
        end if
    end subroutine check_ierr
end module mpiMod
