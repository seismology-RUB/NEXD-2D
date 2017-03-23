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
module triangleNeighborMod
!
! code modified from John Burkards triangle neighbor routine to get the neighbors of a triangle
! modified by LL LL
! 
  implicit none
!
contains
!

subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end subroutine i4col_compare
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end subroutine i4col_sort_a
subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns 
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  N =    ', n
    stop

  end if

  if ( i == j ) then
    return
  end if

  col(1:m) = a(1:m,i)
  a(1:m,i) = a(1:m,j)
  a(1:m,j) = col(1:m)

  return
end subroutine i4col_swap

subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integer ( kind = 4 )s, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end subroutine sort_heap_external
subroutine triangulation_neighbor_triangles ( triangle_order, triangle_num, &
     triangle_node, triangle_neighbor )

  !*****************************************************************************80
  !
  !! TRIANGULATION_NEIGHBOR_TRIANGLES determines triangle neighbors.
  !
  !  Discussion:
  !
  !    A triangulation of a set of nodes can be completely described by
  !    the coordinates of the nodes, and the list of nodes that make up
  !    each triangle.  However, in some cases, it is necessary to know
  !    triangle adjacency information, that is, which triangle, if any,
  !    is adjacent to a given triangle on a particular side.
  !
  !    This routine creates a data structure recording this information.
  !
  !    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
  !    data items.
  !
  !    Note that ROW is a work array allocated dynamically inside this
  !    routine.  It is possible, for very large values of TRIANGLE_NUM,
  !    that the necessary amount of memory will not be accessible, and the
  !    routine will fail.  This is a limitation of the implementation of
  !    dynamic arrays in FORTRAN90.  One way to get around this would be
  !    to require the user to declare ROW in the calling routine
  !    as an allocatable array, get the necessary memory explicitly with
  !    an ALLOCATE statement, and then pass ROW into this routine.
  !
  !    Of course, the point of dynamic arrays was to make it easy to
  !    hide these sorts of temporary work arrays from the poor user!
  !
  !    This routine was revised to store the edge data in a column
  !    array rather than a row array.
  !
  !  Example:
  !
  !    The input information from TRIANGLE_NODE:
  !
  !    Triangle   Nodes
  !    --------   ---------------
  !     1         3      4      1
  !     2         3      1      2
  !     3         3      2      8
  !     4         2      1      5
  !     5         8      2     13
  !     6         8     13      9
  !     7         3      8      9
  !     8        13      2      5
  !     9         9     13      7
  !    10         7     13      5
  !    11         6      7      5
  !    12         9      7      6
  !    13        10      9      6
  !    14         6      5     12
  !    15        11      6     12
  !    16        10      6     11
  !
  !    The output information in TRIANGLE_NEIGHBOR:
  !
  !    Triangle  Neighboring Triangles
  !    --------  ---------------------
  !
  !     1        -1     -1      2
  !     2         1      4      3
  !     3         2      5      7
  !     4         2     -1      8
  !     5         3      8      6
  !     6         5      9      7
  !     7         3      6     -1
  !     8         5      4     10
  !     9         6     10     12
  !    10         9      8     11
  !    11        12     10     14
  !    12         9     11     13
  !    13        -1     12     16
  !    14        11     -1     15
  !    15        16     14     -1
  !    16        13     15     -1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 February 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) TRIANGLE_ORDER, the order of the triangles.
  !
  !    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, integer ( kind = 4 ) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
  !    the nodes that make up each triangle.
  !
  !    Output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the three
  !    triangles that are direct neighbors of a given triangle.  
  !    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1, 
  !    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative 
  !    if there is no neighbor on that side.  In this case, that side of the 
  !    triangle lies on the boundary of the triangulation.
  !
  implicit none

  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_order

  integer ( kind = 4 ) col(4,3*triangle_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) side1
  integer ( kind = 4 ) side2
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
  integer ( kind = 4 ) tri
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  integer ( kind = 4 ) tri1
  integer ( kind = 4 ) tri2
  !
  !  Step 1.
  !  From the list of nodes for triangle T, of the form: (I,J,K)
  !  construct the three neighbor relations:
  !
  !    (I,J,3,T) or (J,I,3,T),
  !    (J,K,1,T) or (K,J,1,T),
  !    (K,I,2,T) or (I,K,2,T)
  !
  !  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
  !
  do tri = 1, triangle_num

     i = triangle_node(1,tri)
     j = triangle_node(2,tri)
     k = triangle_node(3,tri)

     if ( i < j ) then
        col(1:4,3*(tri-1)+1) = (/ i, j, 3, tri /)
     else
        col(1:4,3*(tri-1)+1) = (/ j, i, 3, tri /)
     end if

     if ( j < k ) then
        col(1:4,3*(tri-1)+2) = (/ j, k, 1, tri /)
     else
        col(1:4,3*(tri-1)+2) = (/ k, j, 1, tri /)
     end if

     if ( k < i ) then
        col(1:4,3*(tri-1)+3) = (/ k, i, 2, tri /)
     else
        col(1:4,3*(tri-1)+3) = (/ i, k, 2, tri /)
     end if

  end do
  !
  !  Step 2. Perform an ascending dictionary sort on the neighbor relations.
  !  We only intend to sort on rows 1 and 2; the routine we call here
  !  sorts on rows 1 through 4 but that won't hurt us.
  !
  !  What we need is to find cases where two triangles share an edge.
  !  Say they share an edge defined by the nodes I and J.  Then there are
  !  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
  !  we make sure that these two columns occur consecutively.  That will
  !  make it easy to notice that the triangles are neighbors.
  !
  call i4col_sort_a ( 4, 3*triangle_num, col )
  !
  !  Step 3. Neighboring triangles show up as consecutive columns with
  !  identical first two entries.  Whenever you spot this happening,
  !  make the appropriate entries in TRIANGLE_NEIGHBOR.
  !
  triangle_neighbor(1:3,1:triangle_num) = -1

  icol = 1

  do

     if ( 3 * triangle_num <= icol ) then
        exit
     end if

     if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
        icol = icol + 1
        cycle
     end if

     side1 = col(3,icol)
     tri1 = col(4,icol)
     side2 = col(3,icol+1)
     tri2 = col(4,icol+1)

     triangle_neighbor(side1,tri1) = tri2
     triangle_neighbor(side2,tri2) = tri1

     icol = icol + 2

  end do

  return
end subroutine triangulation_neighbor_triangles
end module triangleNeighborMod
