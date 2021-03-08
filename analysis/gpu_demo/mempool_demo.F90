! Created on 02.02.21.
! Author: Mikhail Zhigun

module m_memory_pool
   use, intrinsic :: iso_c_binding

   implicit none

    private

    type, public  :: stack_memory_pool
        integer(kind=c_int8_t), dimension(:), pointer, public :: mem => null()
        integer(kind=c_int32_t), private        :: end_idx = 0
    contains
        procedure :: create => stack_memory_pool_create
        !procedure :: destroy => stack_memory_pool_destroy
        procedure :: allocate => stack_memory_pool_allocate
        procedure :: deallocate => stack_memory_pool_deallocate
        procedure :: initialized => stack_memory_pool_initialized
        final :: stack_memory_pool_destroy
    end type stack_memory_pool

  contains

    pure function stack_memory_pool_initialized(this) result(initialized)
        class(stack_memory_pool), intent(in) :: this
        logical :: initialized

        initialized = associated(this%mem)
    end function stack_memory_pool_initialized

    subroutine stack_memory_pool_create(this, byte_size)
        class(stack_memory_pool), intent(inout) :: this
        integer(kind=c_int32_t), intent(in) :: byte_size

        if (.not. this%initialized()) then
            allocate(this%mem(byte_size))
            this%end_idx = 0
        else
            stop "stack_memory_pool already initialized"
        end if
    end subroutine stack_memory_pool_create

    subroutine stack_memory_pool_destroy(this)
        type(stack_memory_pool) :: this

        if (this%initialized()) then
            deallocate(this%mem)
            this%end_idx = 0
        end if
    end subroutine stack_memory_pool_destroy

    function stack_memory_pool_allocate(this, byte_size, mem_idx) result(allocated_mem_ptr)
        class(stack_memory_pool), intent(inout) :: this
        integer(kind=c_int32_t), intent(in) :: byte_size
        integer(kind=c_int32_t), intent(out) :: mem_idx
        type(c_ptr) :: allocated_mem_ptr
        integer(kind=c_int32_t) :: used_mem_size
        integer(kind=c_int32_t) :: free_mem_size
        integer(kind=c_int32_t) :: total_mem_size

        if (this%initialized()) then
            total_mem_size = size(this%mem, 1, kind=c_int32_t)
            used_mem_size = this%end_idx
            free_mem_size = total_mem_size - used_mem_size
            if (free_mem_size >=  byte_size) then
                allocated_mem_ptr = c_loc(this%mem(1 + this%end_idx))
                mem_idx = this%end_idx
                this%end_idx = this%end_idx + byte_size
            else
                stop "not enough memory"
            end if
        else
            stop "stack_memory_pool not initialized"
        end if
    end function stack_memory_pool_allocate

    subroutine stack_memory_pool_deallocate(this, allocated_mem_idx)
        class(stack_memory_pool), intent(inout) :: this
        integer(kind=c_int32_t), intent(in) :: allocated_mem_idx

        if (this%initialized()) then
            if (allocated_mem_idx >= 0) then
                if (allocated_mem_idx < this%end_idx) then
                    this%end_idx = allocated_mem_idx
                else
                    !Unfortunately due to separated declaration and initialization
                    !stack deallocations can happen out of order in Fortran
                    !stop "address is within deallocated block"
                end if
            else
                stop "allocated_mem_ptr not in the memory pool"
            end if
        else
            stop "stack_memory_pool not initialized"
        end if
    end subroutine stack_memory_pool_deallocate

end module m_memory_pool

module m_array_int_1d
   use, intrinsic :: iso_c_binding
   use m_memory_pool

   implicit none

    interface sallocate_array
        module procedure sallocate_int_1d
    end interface sallocate_array

    interface sdeallocate_array
        module procedure sdeallocate_int_1d
    end interface sdeallocate_array

    public :: sallocate_array, sdeallocate_array

    type, public :: array_int_1d
         integer(kind=c_int32_t), dimension(:), pointer, private :: arr_ptr => null()
         integer(kind=c_int32_t) :: mempool_idx = 0
         type(stack_memory_pool), pointer, private :: stack_mem_pool => null()
    contains
        procedure :: create => array_int_1d_create
        procedure :: initialized => array_int_1d_initialized
        final :: array_int_1d_destroy
    end type array_int_1d

   contains

    function sallocate_int_1d(stack_mem_pool, size, mem_pool_idx) result(arr_ptr)
        class(stack_memory_pool), intent(inout) :: stack_mem_pool
        integer(kind=c_int32_t), dimension(:), pointer :: arr_ptr
        integer(kind=c_int32_t), intent(out) :: mem_pool_idx
        integer(kind=c_int32_t), intent(in) :: size
        type(c_ptr) :: ptr
        ptr = stack_mem_pool%allocate(size * c_int32_t, mem_pool_idx)
        call c_f_pointer(ptr, arr_ptr, shape=[size])
    end function sallocate_int_1d

    subroutine sdeallocate_int_1d(stack_mem_pool, mem_pool_idx)
        class(stack_memory_pool), intent(inout) :: stack_mem_pool
        integer(kind=c_int32_t), intent(in) :: mem_pool_idx
        integer(kind=c_intptr_t) :: ptr

        call stack_mem_pool%deallocate(mem_pool_idx)
    end subroutine sdeallocate_int_1d

    pure function array_int_1d_initialized(this) result(initialized)
        class(array_int_1d), intent(in) :: this
        logical :: initialized

        initialized = associated(this%arr_ptr)
    end function array_int_1d_initialized

    function array_int_1d_create(this, stack_mem_pool, size) result(arr_ptr)
        class(array_int_1d), intent(inout) :: this
        class(stack_memory_pool), target, intent(inout) :: stack_mem_pool
        integer(kind=c_int32_t), intent(in) :: size
        integer(kind=c_int32_t), dimension(:), pointer :: arr_ptr

        if (.not. this%initialized()) then
            this%arr_ptr => sallocate_array(stack_mem_pool, size, this%mempool_idx)
            arr_ptr => this%arr_ptr
            this%stack_mem_pool => stack_mem_pool
        else
            stop "Array already allocated"
        end if
    end function array_int_1d_create

    subroutine array_int_1d_destroy(this)
        type(array_int_1d) :: this

        if (this%initialized()) then
            call sdeallocate_array(this%stack_mem_pool, this%mempool_idx)
        else
            stop "Array not allocated"
        end if
    end subroutine array_int_1d_destroy

end module m_array_int_1d

module m_mempool_tests
use, intrinsic :: iso_c_binding

implicit none

integer(kind=c_int32_t), parameter :: KB = 1024, ARR_SIZE = 10, N_REPEATS = 1000

contains

subroutine test_auto_arrays()
    ! If compiler decides to allocate on the heap, there will be big performance penalty
    integer(kind=c_int32_t), dimension(ARR_SIZE):: a, b, c
    integer(kind=c_int32_t) :: i, j

    print *, "test_auto_arrays start"

    do i = 1, ARR_SIZE
        a(i) = i
        b(i) = ARR_SIZE - i
    end do

    c = a + b

    print *, a
    print *, b
    print *, c
    print *, "test_auto_arrays end"
end subroutine test_auto_arrays

subroutine test_stack_mempool_manual_dealloc(mem_pool)
    use m_memory_pool
    use m_array_int_1d

    type(stack_memory_pool), intent(inout) :: mem_pool
    integer(kind=c_int32_t), dimension(:), pointer :: a, b, c
    integer(kind=c_int32_t) :: a_mempool_idx, b_mempool_idx, c_mempool_idx
    integer(kind=c_int32_t) :: i, j

    print *, "test_stack_mempool_manual_dealloc start"

    a => sallocate_array(mem_pool, ARR_SIZE, a_mempool_idx)
    b => sallocate_array(mem_pool, ARR_SIZE, b_mempool_idx)
    c => sallocate_array(mem_pool, ARR_SIZE, c_mempool_idx)

    do i = 1, ARR_SIZE
        a(i) = i
        b(i) = ARR_SIZE - i
    end do

    c = a + b

    print *, a
    print *, b
    print *, c
    ! It is dangerous to leave deallocation up to the programmer,
    ! as it could be forgotten
    call sdeallocate_array(mem_pool, c_mempool_idx)
    call sdeallocate_array(mem_pool, b_mempool_idx)
    call sdeallocate_array(mem_pool, a_mempool_idx)

    print *, "test_stack_mempool_manual_dealloc end"
end subroutine test_stack_mempool_manual_dealloc

subroutine test_stack_mempool_auto_dealloc(mem_pool)
    use m_memory_pool
    use m_array_int_1d
    type(stack_memory_pool), intent(inout) :: mem_pool
    type(array_int_1d) :: a_, b_, c_
    integer(kind=c_int32_t), dimension(:), pointer :: a, b, c
    integer(kind=c_int32_t) :: i, j

    print *, "test_stack_mempool_auto_dealloc start"

    a => a_%create(mem_pool, ARR_SIZE)
    b => b_%create(mem_pool, ARR_SIZE)
    c => c_%create(mem_pool, ARR_SIZE)

    do i = 1, ARR_SIZE
        a(i) = i
        b(i) = ARR_SIZE - i
    end do

    c = a + b

    print *, a
    print *, b
    print *, c

    print *, "test_stack_mempool_auto_dealloc end"
end subroutine test_stack_mempool_auto_dealloc

subroutine test_stack_mempool_auto_dealloc_with_acc(mem_pool)
    use m_memory_pool
    use m_array_int_1d
    use openacc
    type(stack_memory_pool), intent(inout) :: mem_pool
    type(array_int_1d) :: a_, b_, c_
    integer(kind=c_int32_t), dimension(:), pointer :: a, b, c
    integer(kind=c_int32_t) :: i, j

    print *, "test_stack_mempool_auto_dealloc_with_acc start"
    ! copyin(mem_pool) should have worked, but it didn't
    ! TODO: find out why
    !$acc data present_or_copyin(mem_pool%mem)
    a => a_%create(mem_pool, ARR_SIZE)
    b => b_%create(mem_pool, ARR_SIZE)
    c => c_%create(mem_pool, ARR_SIZE)

    ! a, b, c point to locations inside memory pool, so there is no need to copy them
    !$acc data present(a, b, c)
    !print *, acc_deviceptr(mem_pool%mem)
    !print *, acc_deviceptr(b)
    !print *, acc_deviceptr(c)
    !$acc parallel default(none) num_gangs(1) num_workers(1) vector_length(128)
    !$acc loop vector
    do i = 1, ARR_SIZE
        a(i) = i
        b(i) = ARR_SIZE - i
    end do

    !$acc loop vector
    do i = 1, ARR_SIZE
        c(i) = a(i) + b(i)
    end do

    !$acc end parallel
    !$acc update host(a, b, c)
    !$acc end data
    !$acc end data

    print *, 'a:', a
    print *, 'b:',b
    print *, 'c:',c

    print *, "test_stack_mempool_auto_dealloc_with_acc end"
end subroutine test_stack_mempool_auto_dealloc_with_acc

end module m_mempool_tests


program mempool_test

use m_mempool_tests
use m_memory_pool

implicit none

type(stack_memory_pool) :: mem_pool
integer(kind=c_int32_t) :: i, j

!------------------------
do j = 1, N_REPEATS
    call test_auto_arrays()
end do
!------------------------
call mem_pool%create(1 * KB)
!------------------------
do j = 1, N_REPEATS
    call test_stack_mempool_manual_dealloc(mem_pool)
end do
!------------------------
do j = 1, N_REPEATS
    call test_stack_mempool_auto_dealloc(mem_pool)
end do
!------------------------
do j = 1, N_REPEATS
    call test_stack_mempool_auto_dealloc_with_acc(mem_pool)
end do

end program mempool_test