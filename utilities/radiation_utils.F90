module radiation_utils
    implicit none

    private :: transpose_real
    integer, parameter :: block_size = 64

    contains

    pure subroutine reorder_dims_1_3(input, output)
      use parkind1, only : jprb
      real(jprb), dimension(:,:,:), intent(in) :: input
      real(jprb), dimension(:,:,:), intent(out) :: output
      integer :: n(1:3), i

      n = shape(input)
      do i = 1, n(2)
        call transpose_real(input(:,i,:), output(:,i,:))
      end do
      !output = reshape(input, shape(output), order=(/3, 2, 1/))
    end subroutine reorder_dims_1_3

    pure subroutine reorder_dims_1_2(input, output)
      use parkind1, only : jprb
      real(jprb), dimension(:,:,:), intent(in) :: input
      real(jprb), dimension(:,:,:), intent(out) :: output
      integer :: n(1:3), i

      n = shape(input)
      do i = 1, n(3)
        call transpose_real(input(:,:,i), output(:,:,i))
      end do
      !output = reshape(input, shape(output), order=(/2, 1, 3/))
    end subroutine reorder_dims_1_2

    pure subroutine reverse_dim_2(arr)
      use parkind1, only : jprb
      real(jprb), dimension(:,:,:), intent(inout) :: arr
      integer, parameter :: buf_size = 8
      real(jprb), dimension(buf_size) :: buf
      real(jprb) :: val
      integer :: n_cols, n_rows, n_lvls, ci, ri, li, cbi, eri, eci

      n_cols = size(arr, 1)
      n_rows = size(arr, 2)
      n_lvls = size(arr, 3)
      do li = 1, n_lvls
          ri = 1
          eri = n_rows
          do while (eri .gt. ri)
            ci = 1
            eci = ci + buf_size - 1
            do while(eci .le. n_cols)
              call swap(arr(ci:eci, ri, li),  arr(ci:eci, eri, li), buf)
              ci = ci + buf_size
              eci = eci + buf_size
            end do
            do ci = ci, n_cols
              val = arr(ci, ri, li)
              arr(ci, ri, li) = arr(ci, eri, li)
              arr(ci, eri, li) = val
            end do
            ri = ri + 1
            eri = eri - 1
          end do
      end do
    contains
      pure subroutine swap(arr1, arr2, buf)
        real(jprb), dimension(buf_size), intent(inout) :: arr1, arr2, buf
        buf = arr1
        arr1 = arr2
        arr2 = buf
      end subroutine swap
    end subroutine reverse_dim_2

    !This subroutine performs cache-friendly transposition by processing the matrix in "small" blocks of 64x64
    pure subroutine transpose_real(input, output)
      use parkind1, only : jprb
      real(jprb), dimension(:,:), intent(in) :: input
      real(jprb), dimension(:,:), intent(out) :: output
      real(jprb), dimension(block_size,block_size) :: from, to
      integer :: ri ,ci, rbi, cbi, n_r_blocks, n_c_blocks, n_rows, n_cols, i, rem_rows, rem_cols

      n_rows = size(input, 2)
      n_cols = size(input, 1)
      n_r_blocks = n_rows / block_size
      rem_rows = mod(n_rows, block_size)
      n_c_blocks = n_cols / block_size
      rem_cols = mod(n_cols, block_size)
      ri = 1
      do rbi = 1, n_r_blocks
        ci = 1
        do cbi = 1, n_c_blocks
          from = input(ci:ci+block_size-1, ri:ri+block_size-1)
          to = transpose(from)
          output(ri:ri+block_size-1, ci:ci+block_size-1) = to
          ci = ci + block_size
        end do
        if (rem_cols > 0) then
            from(1:rem_cols, :) = input(ci:n_cols, ri:ri+block_size-1)
            to(:, 1:rem_cols) = transpose(from(1:rem_cols, :))
            output(ri:ri+block_size-1, ci:n_cols) = to(:, 1:rem_cols)
        end if
        ri = ri + block_size
      end do
      if (rem_rows > 0) then
        ci = 1
        do cbi = 1, n_c_blocks
          from(:, 1:rem_rows) = input(ci:ci+block_size-1, ri:n_rows)
          to(1:rem_rows, :) = transpose(from(:, 1:rem_rows))
          output(ri:n_rows, ci:ci+block_size-1) = to(1:rem_rows, :)
          ci = ci + block_size
        end do
        if (rem_cols > 0) then
            from(1:rem_cols, 1:rem_rows) = input(ci:n_cols, ri:n_rows)
            to(1:rem_rows, 1:rem_cols) = transpose(from(1:rem_cols, :))
            output(ri:n_rows, ci:n_cols) = to(1:rem_rows, 1:rem_cols)
        end if
      end if

!      do ri = 1, n_rows
!        do ci = 1, n_cols
!          output(ri, ci) = input(ci, ri)
!        end do
!      end do
      !output = transpose(input)
    end subroutine transpose_real

end module radiation_utils
