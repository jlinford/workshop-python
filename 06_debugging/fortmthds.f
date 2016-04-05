      subroutine buildgrids(time)
      implicit none
      double precision, intent(in) :: time

      double precision :: dummy1,dummy2
      integer :: i,n

      dummy1 = 0.001
      dummy2 = 1.0

      n = 100000
      do i = 1,n
         dummy1 = dummy1 + dummy2
         dummy1 = dummy1 - dummy2
      end do

      return
      end
!
! ------------------------
!
      subroutine regrid(time)
      implicit none
      double precision, intent(in) :: time

      double precision :: dummy1,dummy2
      integer :: i,n

      dummy1 = 0.001
      dummy2 = 1.0

      n = 100000
      do i = 1,n
         dummy1 = dummy1 + dummy2
         dummy1 = dummy1 - dummy2
      end do

      return
      end

!
! ------------------------
!
      subroutine timestep(time,dt,x,len)
      implicit none
      double precision, intent(in) :: time,dt
      integer, intent(in) :: len
      double precision :: x(*)

      double precision :: dummy1,dummy2
      integer :: i,n

      dummy1 = 0.001
      dummy2 = 1.0

      n = 100000
      do i = 1,n
         dummy1 = dummy1 + dummy2
         dummy1 = dummy1 - dummy2
      end do

      !x(len+5) = 1

      return
      end


