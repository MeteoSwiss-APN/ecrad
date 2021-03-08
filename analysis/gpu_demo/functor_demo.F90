! Created on 03.02.21.
! Author: Mikhail Zhigun

module m_functor
 use :: iso_c_binding

 type, public :: cube_functor
    integer(kind=c_int32_t), dimension(:), pointer, private :: buf => null()
 contains
    procedure :: create => cube_functor_create
    procedure :: run => cube_functor_run
    final :: cube_functor_destroy
 end type cube_functor

contains
  subroutine cube_functor_create(this, vec_len)
    ! With OpenACC and mempool, this function will be called only once
    ! on host, to (stack-)allocate memory for the local variables on device
    class(cube_functor), intent(inout) :: this
    integer(kind=c_int32_t), intent(in) :: vec_len

    if (.not. associated(this%buf)) then
        allocate(this%buf(vec_len))
    end if
  end subroutine cube_functor_create

  subroutine cube_functor_destroy(this)
    type(cube_functor) :: this

    if (associated(this%buf)) then
      deallocate(this%buf)
    end if
  end subroutine cube_functor_destroy

  subroutine cube_functor_run(this, in, out)
    class(cube_functor), intent(inout) :: this
    integer(kind=c_int32_t), dimension(:), pointer, intent(in) :: in
    integer(kind=c_int32_t), dimension(:), pointer, intent(out) :: out
    integer(kind=c_int32_t) :: arr_len, i
    integer(kind=c_int32_t), dimension(:), pointer :: buf
    buf => this%buf
    arr_len = size(buf, 1)
    ! This is the downside of pointers in Fortran, they don't let you specify size
    ! But with automatic arrays there is no garantie that compiler will not silently accept
    ! array of lesser or bigger size
    if (size(in, 1) .ne. arr_len) then
        stop "array size mismatch"
    end if
    if (size(out, 1) .ne. arr_len) then
        stop "array size mismatch"
    end if

    do i = 1, arr_len
        buf(i) = in(i) * in(i)
        buf(i) = buf(i) * in(i)
        out(i) = buf(i)
    end do

  end subroutine cube_functor_run

  subroutine cube(arr_len, in, out)
    integer(kind=c_int32_t), intent(in)  :: arr_len
    !Auto-arrays will not play well with OpenACC, because of potential dynamic allocations
    !have to be replaced with pointers
    integer(kind=c_int32_t), intent(in) :: in(1:arr_len)
    integer(kind=c_int32_t), intent(out) :: out(1:arr_len)
    !Scalar variables are highly unlikely to be dynamically allocated, so they can remain
    integer(kind=c_int32_t) :: i
    !Local variables of array type have to be moved to the data section so that they can
    !be allocated before function is run
    integer(kind=c_int32_t) :: buf(1:arr_len)

    !Code can stay exactly as it was
    do i = 1, arr_len
        buf(i) = in(i) * in(i)
        buf(i) = buf(i) * in(i)
        out(i) = buf(i)
    end do

  end subroutine cube
end module m_functor

program functor_demo
    use m_functor
    integer(kind=c_int32_t), parameter :: size = 10, num_repeats = 10

    call test_cube(size, num_repeats)
    call test_cube_functor(size, num_repeats)
contains
    subroutine test_cube(size, num_repeats)
        integer(kind=c_int32_t) , intent(in):: size, num_repeats
        integer(kind=c_int32_t) :: x(1:size), r(1:size), i
        print *, "test_cube_start"

        x = 2
        r = 0
        do i = 1, num_repeats
            call cube(size, x, r)
        end do
        print *, "x:", x
        print *, "r:", r
        print *, "test_cube_end"
    end subroutine test_cube

    subroutine test_cube_functor(size, num_repeats)
        integer(kind=c_int32_t) , intent(in):: size, num_repeats
        integer(kind=c_int32_t), dimension(:), pointer :: x, r
        type(cube_functor) :: cube
        print *, "test_cube_functor_start"

        allocate(x(size))
        allocate(r(size))

        x = 2
        r = 0
        call cube%create(size)
        do i = 1, num_repeats
            call cube%run(x, r)
        end do
        print *, "x:", x
        print *, "r:", r
        print *, "test_cube_functor_end"
    end subroutine test_cube_functor



end program