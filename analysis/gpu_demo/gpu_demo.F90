! Created by  on 28.01.21.
! Author: Mikhail Zhigun

!#define DEBUG_ON_CPU

module m_gpu_demo
  use openacc
  use cudafor
    integer, parameter :: KB = 1024
    integer, parameter :: MB = 1024 * KB
    integer, parameter :: GB = 1024 * MB
    integer, parameter :: NUM_COLS = 20000
	integer, parameter :: NUM_LEVELS = 137
	integer, parameter :: NUM_G_SW = 112
	integer, parameter :: NUM_G_LW = 140
	integer, parameter :: INPUT_BYTE_SIZE = 1 * GB
	integer, parameter :: WORKING_SET_BYTE_SIZE_SW = 2.5 * MB
	integer, parameter :: WORKING_SET_BYTE_SIZE_LW = 2.5 * MB
    integer, parameter :: SIZEOF_DOUBLE = 8
    integer, parameter :: INPUT_DATA_SIZE_PER_COL = INPUT_BYTE_SIZE / (NUM_COLS * NUM_LEVELS * SIZEOF_DOUBLE)
    integer, parameter :: WORKING_SET_SW_DATA_SIZE_PER_COL = WORKING_SET_BYTE_SIZE_SW / (NUM_G_SW * NUM_LEVELS * SIZEOF_DOUBLE)
    integer, parameter :: WORKING_SET_LW_DATA_SIZE_PER_COL = WORKING_SET_BYTE_SIZE_LW / (NUM_G_LW * NUM_LEVELS * SIZEOF_DOUBLE)
contains
  subroutine run()
    !-----------------------------------------------------------------------------------------
    integer :: threads_per_block_sw, threads_per_block_lw
    integer :: num_thread_blocks_sw, num_thread_blocks_lw
    integer :: num_warps_per_block_sw, num_warps_per_block_lw
    integer :: num_thread_blocks_per_mp_sw, num_thread_blocks_per_mp_lw
    integer :: num_cols_per_sw_thread_block, num_cols_per_lw_thread_block
    real, dimension(:, :, :), allocatable :: input ! INPUT_DATA_SIZE_PER_COL x levels x cols
    real, dimension(:, :, :, :), allocatable :: working_set_sw ! NUM_G_SW x WORKING_SET_SW_DATA_SIZE_PER_COL x NUM_LEVELS x num_thread_blocks_sw
    real, dimension(:, :, :, :), allocatable :: working_set_lw ! NUM_G_LW x WORKING_SET_LW_DATA_SIZE_PER_COL x NUM_LEVELS x num_thread_blocks_lw
#ifndef DEBUG_ON_CPU
    !$acc declare device_resident(input, working_set_sw, working_set_lw)
#endif
    integer, dimension(:), allocatable :: start_col_sw ! num_thread_blocks_sw
    integer, dimension(:), allocatable :: end_col_sw ! num_thread_blocks_sw
    integer, dimension(:), allocatable :: start_col_lw ! num_thread_blocks_lw
    integer, dimension(:), allocatable :: end_col_lw ! num_thread_blocks_lw
    integer, dimension(:, :), allocatable :: thread_block_idx_by_mp_sw ! num_thread_blocks_per_mp_sw x num_mp
    integer, dimension(:, :), allocatable :: thread_block_idx_by_mp_lw ! num_thread_blocks_per_mp_lw x num_mp
    integer :: i, mp_i, blk_i, t_blk_i, start_col, end_col, col_i, lvl_i, data_i, g_i, r_i, col_wset_i
    !-----------------------------------------------------------------------------------------
    integer, parameter :: DEV_NUM = 0 ! Always pick up first found device
    integer :: num_dev, istat
    character*100 :: dev_name
    integer(8) :: gpu_mem_size, gpu_mem_size_mb, num_mp, max_threads_per_mp, max_threads_per_block, warp_size
    integer(kind=cuda_count_kind) :: used_gpu_mem, total_gpu_mem, free_gpu_mem
    integer :: cuda_status 
    type (cudaDeviceProp) :: prop

    num_dev = acc_get_num_devices(ACC_DEVICE_NVIDIA)
    if (num_dev .lt. 1) then
      print *, 'Error: There are no devices available on this host'
      stop
    endif
    write(*,"(' Number of Nvidia GPU devices: ',i0)") num_dev
    call acc_get_property_string(DEV_NUM, ACC_DEVICE_NVIDIA, ACC_PROPERTY_NAME, dev_name)
    write(*,"(' Selecting device:             ',i0,' ""',a,'""')") DEV_NUM, dev_name
    call acc_set_device_num(DEV_NUM, acc_device_nvidia)
    gpu_mem_size = acc_get_property(DEV_NUM, ACC_DEVICE_NVIDIA, ACC_PROPERTY_MEMORY)
    gpu_mem_size_mb = gpu_mem_size / (1024*1024)
    write(*,"(' Device memory size:           ',i0,'mb')") gpu_mem_size_mb
    istat = cudaGetDeviceProperties(prop, DEV_NUM)
    num_mp = prop%multiProcessorCount
    max_threads_per_block = prop%maxThreadsPerBlock
    warp_size = prop%warpSize
    max_threads_per_mp = prop%maxThreadsPerMultiProcessor
    write(*,"(' Number of multiprocessors:    ',i0)") num_mp
    write(*,"(' Max threads per MP:           ',i0)") max_threads_per_mp
    write(*,"(' Max threads per block:        ',i0)") max_threads_per_block
    write(*,"(' Warp (OpenACC vector/num_workers) size:   ',i0)") warp_size
    write(*,"('---------------------------------------')")
    threads_per_block_sw = round_up_to_mod(NUM_G_SW,  int(warp_size))
    threads_per_block_lw = round_up_to_mod(NUM_G_LW, int(warp_size))
    num_warps_per_block_sw = threads_per_block_sw / warp_size
    num_warps_per_block_lw = threads_per_block_lw / warp_size
    num_thread_blocks_per_mp_sw = (max_threads_per_mp / threads_per_block_sw)
    num_thread_blocks_per_mp_lw = (max_threads_per_mp / threads_per_block_lw)
    num_thread_blocks_sw = num_mp * num_thread_blocks_per_mp_sw
    num_thread_blocks_lw = num_mp * num_thread_blocks_per_mp_lw
    num_cols_per_sw_thread_block = NUM_COLS / num_thread_blocks_sw ! Imprecise, some cols in the end maybe left out
    num_cols_per_lw_thread_block = NUM_COLS / num_thread_blocks_lw ! last thread block should get it share of cols + rounded of calls in the end
    write(*,"(' threads_per_block_sw:    ',i0)") threads_per_block_sw
    write(*,"(' num_warps_per_block_sw:    ',i0)") num_warps_per_block_sw
    write(*,"(' num_thread_blocks_sw:    ',i0)") num_thread_blocks_sw
    write(*,"(' num_thread_blocks_per_mp_sw:    ',i0)") num_thread_blocks_per_mp_sw
    write(*,"(' num_cols_per_sw_thread_block:    ',i0)") num_cols_per_sw_thread_block
    write(*,"(' threads_per_block_lw:    ',i0)") threads_per_block_lw
    write(*,"(' num_warps_per_block_lw:    ',i0)") num_warps_per_block_lw
    write(*,"(' num_thread_blocks_lw:    ',i0)") num_thread_blocks_lw
    write(*,"(' num_thread_blocks_per_mp_lw:    ',i0)") num_thread_blocks_per_mp_lw
    write(*,"(' num_cols_per_lw_thread_block:    ',i0)") num_cols_per_lw_thread_block
    write(*,"(' INPUT_DATA_SIZE_PER_COL:    ',i0)") INPUT_DATA_SIZE_PER_COL
    write(*,"(' WORKING_SET_SW_DATA_SIZE_PER_COL:    ',i0)") WORKING_SET_SW_DATA_SIZE_PER_COL
    write(*,"(' WORKING_SET_LW_DATA_SIZE_PER_COL:    ',i0)") WORKING_SET_LW_DATA_SIZE_PER_COL
    !-----------------------------
    cuda_status = cudaMemGetInfo(free_gpu_mem, total_gpu_mem)
    used_gpu_mem = total_gpu_mem - free_gpu_mem
    write(*,"(' Used GPU memory before allocation:    ',i0,'mb')") used_gpu_mem / MB
    !-----------------------------
    !These allocations should happen directly on GPU due to "device resident" directive
    allocate(input(INPUT_DATA_SIZE_PER_COL, NUM_LEVELS, NUM_COLS))
    allocate(working_set_sw(NUM_G_SW, WORKING_SET_SW_DATA_SIZE_PER_COL, NUM_LEVELS, num_thread_blocks_sw))
    allocate(working_set_lw(NUM_G_LW, WORKING_SET_LW_DATA_SIZE_PER_COL, NUM_LEVELS, num_thread_blocks_lw))
    !-----------------------------
    cuda_status = cudaMemGetInfo(free_gpu_mem, total_gpu_mem)
    used_gpu_mem = total_gpu_mem - free_gpu_mem
    write(*,"(' Used GPU memory after allocation:    ',i0,'mb')") used_gpu_mem / MB
    !-----------------------------
    !These allocations happen on CPU, so data should be copied after
    allocate(start_col_sw(num_thread_blocks_sw))
    allocate(end_col_sw(num_thread_blocks_sw))
    allocate(start_col_lw(num_thread_blocks_lw))
    allocate(end_col_lw(num_thread_blocks_lw))
    allocate(thread_block_idx_by_mp_sw(num_thread_blocks_per_mp_sw, num_mp))
    allocate(thread_block_idx_by_mp_lw(num_thread_blocks_per_mp_lw, num_mp))
    !-------------------------------------------------------
    do i = 1, num_thread_blocks_sw
      start_col_sw(i) = num_cols_per_sw_thread_block  * (i - 1) + 1
      end_col_sw(i) = start_col_sw(i) + num_cols_per_sw_thread_block - 1
    end do
    end_col_sw(num_thread_blocks_sw) = num_cols
    do i = 1, num_thread_blocks_lw
      start_col_lw(i) = num_cols_per_lw_thread_block  * (i - 1) + 1
      end_col_lw(i) = start_col_lw(i) + num_cols_per_lw_thread_block - 1
    end do
    end_col_lw(num_thread_blocks_lw) = num_cols
    i = 1
    DO mp_i = 1, num_mp
      DO blk_i = 1, num_thread_blocks_per_mp_sw
        thread_block_idx_by_mp_sw(blk_i, mp_i) = i
        i = i + 1
	  END DO
	END DO
    i = 1
    DO mp_i = 1, num_mp
      DO blk_i = 1, num_thread_blocks_per_mp_lw
        thread_block_idx_by_mp_lw(blk_i, mp_i) = i
        i = i + 1
	  END DO
	END DO
    call solver_sw(blk_i, col_i, col_wset_i, end_col, end_col_sw, g_i, input, lvl_i, mp_i, &
                   num_mp, threads_per_block_sw, num_warps_per_block_sw, num_thread_blocks_sw, &
                   num_thread_blocks_per_mp_sw, r_i, start_col, start_col_sw, t_blk_i, thread_block_idx_by_mp_sw, &
                   working_set_sw)
    call solver_lw(blk_i, col_i, col_wset_i, end_col, end_col_lw, g_i, input, lvl_i, mp_i, &
                   num_mp, threads_per_block_lw, num_warps_per_block_lw, num_thread_blocks_lw, &
                   num_thread_blocks_per_mp_lw, r_i, start_col, start_col_lw, t_blk_i, thread_block_idx_by_mp_lw, &
                   working_set_lw)
  end subroutine run

  subroutine solver_sw(blk_i, col_i, col_wset_i, end_col, end_col_sw, g_i, input, lvl_i, mp_i, &
                       num_mp, threads_per_block_sw, num_warps_per_block_sw, num_thread_blocks_sw, &
                       num_thread_blocks_per_mp_sw, r_i, start_col, start_col_sw, t_blk_i, thread_block_idx_by_mp_sw, &
                       working_set_sw)
      implicit none
      integer :: blk_i
      integer :: col_i
      integer :: col_wset_i
      integer :: end_col
      integer, allocatable :: end_col_sw(:)
      integer :: g_i
      real, allocatable :: input(:, :, :)
      integer :: lvl_i
      integer :: mp_i
      integer(8) :: num_mp
      integer :: threads_per_block_sw
      integer :: num_warps_per_block_sw
      integer :: num_thread_blocks_sw
      integer :: num_thread_blocks_per_mp_sw
      integer :: r_i
      integer :: start_col
      integer, allocatable :: start_col_sw(:)
      integer :: t_blk_i, t_blk_i_counter
      integer, allocatable :: thread_block_idx_by_mp_sw(:, :)
      real, allocatable :: working_set_sw(:, :, : , :)
      integer :: gang_i

      t_blk_i_counter = 1
#ifndef DEBUG_ON_CPU
      !$acc data present(input, working_set_sw) copyin(thread_block_idx_by_mp_sw, start_col_sw, end_col_sw, t_blk_i)
      !$acc parallel default(none) num_gangs(num_thread_blocks_sw) num_workers(num_warps_per_block_sw) copyin(t_blk_i_counter) &
      !$acc& vector_length(threads_per_block_sw)
      !$acc loop independent gang
#endif
      DO gang_i = 1, num_thread_blocks_sw
          ! For unknown reason simple "t_blk_i = gang_i" does not work, because gang_i is always equal to 1
          t_blk_i = allocate_gang_unique_idx(t_blk_i_counter) ! this can be replaced by non-portable __pgi_gangidx() + 1
          start_col = start_col_sw(t_blk_i)
          end_col = end_col_sw(t_blk_i)
          !visited(t_blk_i) = t_blk_i
          !print *, t_blk_i, start_col, end_col
          !write(*,"(' t_blk_i:   ',i0)") t_blk_i
          !write(*,"('   start_col: ',i0)") start_col
          !write(*,"('   end_col:   ',i0)") end_col
          !Input IO
          !$acc loop independent seq
          DO col_i = start_col, end_col
            !$acc loop seq
            DO lvl_i = 1, NUM_LEVELS
              ! This should be in parallel
              input(:, lvl_i, col_i) = 1.
              input(:, lvl_i, col_i) = input(:, lvl_i, col_i) + 2.
            END DO
          END DO
          !Computation
          !$acc loop independent seq
          DO col_i = start_col, end_col
            col_wset_i = col_i - start_col + 1
            !print *, shape(working_set_sw)
            !print *, ' ', col_wset_i
            !$acc loop seq
            DO lvl_i = 1, NUM_LEVELS
              !$acc loop seq
              DO r_i = 1, WORKING_SET_SW_DATA_SIZE_PER_COL
                ! This should be in parallel
                working_set_sw(:, r_i, lvl_i, t_blk_i) = 1.
                !$acc loop vector
                do g_i = 1, NUM_G_SW
                  working_set_sw(g_i, r_i, lvl_i, t_blk_i) = exp(working_set_sw(g_i, r_i, lvl_i, t_blk_i))
                end do
              END DO
            END DO
          END DO
          !Input IO
          !$acc loop independent seq
          DO col_i = start_col, end_col
            !$acc loop seq
            DO lvl_i = 1, NUM_LEVELS
              ! This should be in parallel
              input(:, lvl_i, col_i) = 2.
            END DO
          END DO
      END DO
#ifndef DEBUG_ON_CPU
      !$acc end parallel
      !$acc end data
 #endif
      !print *, visited
  end subroutine

  subroutine solver_lw(blk_i, col_i, col_wset_i, end_col, end_col_lw, g_i, input, lvl_i, mp_i, &
                       num_mp, threads_per_block_lw, num_warps_per_block_lw, num_thread_blocks_lw, &
                       num_thread_blocks_per_mp_lw, r_i, start_col, start_col_lw, t_blk_i, thread_block_idx_by_mp_lw, &
                       working_set_lw)
      implicit none
      integer :: blk_i
      integer :: col_i
      integer :: col_wset_i
      integer :: end_col
      integer, allocatable :: end_col_lw(:)
      integer :: g_i
      real, allocatable :: input(:, :, :)
      integer :: lvl_i
      integer :: mp_i
      integer(8) :: num_mp
      integer :: threads_per_block_lw
      integer :: num_warps_per_block_lw
      integer :: num_thread_blocks_lw
      integer :: num_thread_blocks_per_mp_lw
      integer :: r_i
      integer :: start_col
      integer, allocatable :: start_col_lw(:)
      integer :: t_blk_i, t_blk_i_counter
      integer, allocatable :: thread_block_idx_by_mp_lw(:, :)
      real, allocatable :: working_set_lw(:, :, : , :)
      integer :: gang_i

      t_blk_i_counter = 1
#ifndef DEBUG_ON_CPU
      !$acc data present(input, working_set_lw) copyin(thread_block_idx_by_mp_lw, start_col_lw, end_col_lw, t_blk_i)
      !$acc parallel num_gangs(num_thread_blocks_lw) num_workers(num_warps_per_block_lw) copyin(t_blk_i_counter) &
      !$acc& vector_length(threads_per_block_lw)
      !$acc loop independent gang
#endif
      DO gang_i = 1, num_thread_blocks_lw
          ! For unknown reason simple "t_blk_i = gang_i" does not work, because gang_i is always equal to 1
          t_blk_i = allocate_gang_unique_idx(t_blk_i_counter) ! this can be replaced by non-portable __pgi_gangidx() + 1
          start_col = start_col_lw(t_blk_i)
          end_col = end_col_lw(t_blk_i)
          !visited(t_blk_i) = t_blk_i
          !print *, t_blk_i, start_col, end_col
          !write(*,"(' t_blk_i:   ',i0)") t_blk_i
          !write(*,"('   start_col: ',i0)") start_col
          !write(*,"('   end_col:   ',i0)") end_col
          !Input IO
          !$acc loop independent seq
          DO col_i = start_col, end_col
            !$acc loop seq
            DO lvl_i = 1, NUM_LEVELS
              ! This should be in parallel
              input(:, lvl_i, col_i) = 1.
              input(:, lvl_i, col_i) = input(:, lvl_i, col_i) + 2.
            END DO
          END DO
          !Computation
          !$acc loop independent seq
          DO col_i = start_col, end_col
            col_wset_i = col_i - start_col + 1
            !print *, shape(working_set_lw)
            !print *, ' ', col_wset_i
            !$acc loop seq
            DO lvl_i = 1, NUM_LEVELS
              !$acc loop seq
              DO r_i = 1, WORKING_SET_lw_DATA_SIZE_PER_COL
                ! This should be in parallel
                working_set_lw(:, r_i, lvl_i, t_blk_i) = 1.
                !$acc loop vector
                do g_i = 1, NUM_G_lw
                  working_set_lw(g_i, r_i, lvl_i, t_blk_i) = exp(working_set_lw(g_i, r_i, lvl_i, t_blk_i))
                end do
              END DO
            END DO
          END DO
          !Input IO
          !$acc loop independent seq
          DO col_i = start_col, end_col
            !$acc loop seq
            DO lvl_i = 1, NUM_LEVELS
              ! This should be in parallel
              input(:, lvl_i, col_i) = 2.
            END DO
          END DO
      END DO
#ifndef DEBUG_ON_CPU
      !$acc end parallel
      !$acc end data
 #endif
      !print *, visited
  end subroutine

  pure function round_up_to_mod(val, mod) result(res)
    integer(4), intent(in) :: val, mod
    integer(4) :: res
    res = ((val + mod - 1) / mod) * mod
  end function round_up_to_mod

  integer function allocate_gang_unique_idx(counter_var) result(idx)
      !$acc routine vector
      integer, intent(inout) :: counter_var
      integer :: idx, v_i
      !$acc loop vector private(v_i)
      do v_i = 1, 1
          !$acc atomic capture
          idx = counter_var
          counter_var = counter_var + 1
          !$acc end atomic
      end do
  end function allocate_gang_unique_idx

end module m_gpu_demo

program gpu_demo
    use m_gpu_demo
#ifdef _OPENACC
  use m_gpu_demo
implicit none
  call run()
#else
  print *, "Error: compiled without OpenACC"
  stop
#endif
end program gpu_demo
