!Version 1.0.1 Oct 05, 2010


!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data Structure
!!!!!!!!!!!!!!!!!!!!!!!!!
#include "PRESET.INI"

module data_structure


    implicit none
    include 'mpif.h'
    !!!!!!!!!!!!!!!!
    ! Atoms and Chains Data Structure
    !!!!!!!!!!!!!!!
    integer :: DIM    !Dimension
    integer :: NON    !Number of Nodes (CPU)
    integer :: NOS    !Max Number of atoms in each node
    integer :: NOSEA  !Max Number of atoms to be sent
    integer :: NOREA  !Max Number of atoms to be received
    integer :: NOA    !Actual number of atoms in each node
    integer :: NOAA   !Actual number of active atoms in each node
    integer :: NOPA   !Actual number of pinned atoms  
    integer :: NOTA   !Total Number of atoms in the system
    integer :: NOTAA  !Total Number of active atoms in the system 
    integer, parameter :: Nlist_RAM = 2147483640 !2147483640
    double precision :: PI = 3.1415926535897932384626433832795d0
    
    !Dynamical Variables
    double precision :: dt, dt2, dt3, dt4, dt5        !time increment and 0.5dt^2
    double precision :: dt2d2, dt3d6, dt4d24, dt5d120
    double precision :: f_cut            !forcdt4d4,e_cutoff distance
    double precision :: f_cut_2
    double precision :: f_cut_sq         !force_cutoff ^ 2
    double precision :: f_cut_vect(3)
    double precision :: f_cut_vect_2(3)
    double precision :: f_cut_sq_max(0:4)
    
    !Quick Square Root Table
    double precision :: qsqrt_max, qsqrt_min, qsqrt_inc, qsqrt_n
    double precision :: qtbl1(500), qtbl2(500)
    
    !Run Time Variables
    integer :: time_index, end_time, movie_time, sync_time   !number of cycles to be run, rate at which movies are recorded
    integer, allocatable, dimension(:), save :: status
    logical :: free_pin
    
    !Run Time Batch Mode Variable
    logical :: batch_mode
    logical :: cmd_new_job, cmd_quit
    integer :: my_job_id, my_job_index, my_job_version, my_job_force_mode
    integer :: myfile_scr01
    integer :: root_len, frame_counter
    integer :: n_para_int, n_para_double
    logical :: write_stress
    character*100 :: root_path
    double precision :: para_vect(3)
    logical, allocatable, dimension(:), save :: job_id
    integer, allocatable, dimension(:), save :: para_int
    double precision, allocatable, dimension(:), save :: para_double
    
        
    !parameters associate with an atom
    double precision, allocatable,  dimension (:, :), save :: r0, r1, r2, r3, r4, r5   !Position and its derivitive
    double precision, allocatable,  dimension (:), save :: r1_ave
    double precision :: a0, a1, a3, a4, a5
    integer, allocatable, dimension (:), save  :: t, s, id     !t : type of atom, s : status of atom
    integer, allocatable, dimension (:), save  :: s1_list, s1_list_inv
    
    !Dynamic XP EAM Energy parameters
    integer, allocatable, dimension (:), save  :: EAM_time_stamp
    double precision, allocatable,  dimension (:), save :: rho_sum_xp, Vi_sum_xp
    double precision, allocatable,  dimension (:, :), save :: dFi_sum_xp, dVi_sum_xp
    
    !Dynamic parameters
    double precision, allocatable,  dimension (:, :), save :: dr2          !difference in force
    double precision, allocatable,  dimension (:, :), save :: dr2_ex       !applied stress 
    double precision, allocatable,  dimension (:, :), save :: dr0_ex       !applied strain
    double precision, allocatable,  dimension (:, :), save :: Force        !Force Vector
    double precision, allocatable,  dimension (:), save :: Energy           !Energy for each atom
    double precision, allocatable,  dimension (:), save :: Iv1, Iv2, Iv3, Jv2, Jv3
    double precision, allocatable,  dimension (:, :, :), save :: Stress, Stress_ave     !Stress
    double precision :: f_cen(3)
    double precision :: m_d3k(2)                                           !m/3k from 1/2 m v^2 = 3/2 k T (in Angstrom & femtosecond)
    double precision :: m_atom(2)
    double precision :: m_dq(2), q_dm(2)
    double precision :: T_current, T_wanted, T_ratio_sq
    double precision :: time_exp
    
    !Force & Energy Spline parameters
    integer :: NOP_Spline  !Number of data point of Spline interpolation
    double precision, allocatable, dimension (:, :, :), save :: pot_V, for_V
    double precision, allocatable, dimension (:, :, :), save :: pot_F, for_F
    double precision, allocatable, dimension (:, :, :), save :: rho, drho
    double precision, allocatable, dimension (:), save :: v_min, v_max, f_max, rho_min, rho_max  !the domain range of the functions
    double precision, allocatable, dimension (:), save :: v_inc, f_inc, rho_inc !increment sizes of the domains
    double precision, allocatable, dimension (:), save :: v_inc_inv, f_inc_inv, rho_inc_inv  !reciprocal of increments
    
    !chain list parameters
    integer, allocatable, dimension (:), save :: next, prev  !record the index of the next & previous element in a chain, 0=start, -1=end
    integer, allocatable, dimension (:), save :: ch_empty    !record which element in the array can be replaced
    integer, allocatable, dimension (:), save :: ch_list    !record which element in the array can be replaced
    integer :: s1_count, sn0_count
    integer :: ch_start, ch_end, ch_len     !start, end and length of the chain
    integer :: ch_empty_count               !number of elements in the array which are empty
    
    !Neighbor List Parameters
    integer, allocatable, dimension (:,:,:), save :: nbox_start !keep record of the head of the chain list
    integer, allocatable, dimension (:,:,:), save :: nbox_count !keep a record of how many elements in that cell
    integer, allocatable, dimension (:,:,:), save :: nbox_time  !keep track of the latest update time of that cell
    integer, allocatable, dimension (:), save :: nbox_next      !next element of the neighbor list chain
    double precision nbox_width(3)
    
    !MPI Parameters
    integer myid, mygpid, MPI_MD_GROUP, MPI_MD_NEW_GROUP, MPI_MD_WORLD         ! The id assigned to each node
    integer, allocatable, dimension (:), save :: block_size3, block_size1, block_index
    integer, allocatable, dimension (:), save :: gpid, wdid
    double precision vect_one(3)
    
    !receive paremters
    double precision , allocatable, dimension (:, :), save ::rshadow_r0, ratom_r0, ratom_r1, ratom_r2,ratom_r3,ratom_r4,ratom_r5
    integer, allocatable, dimension (:), save  :: rshadow_t, rshadow_id, ratom_t, ratom_s, ratom_id !list of received atoms
    integer rshadow_count,  ratom_count
    logical should_run_forward, should_run_reverse
    
    !send parameters
    integer, allocatable, dimension (:), save :: sshadow_list, satom_list, image_list !list of atoms pending to be sent
    integer sshadow_count, satom_count, image_count
    
    !!!!!!!!!!!!!!!!!
    ! Box and Node Data Structure
    !!!!!!!!!!!!!!!!!
    double precision, allocatable, dimension (:, :, :), save :: br !Each Box's Coordinates
    double precision, allocatable, dimension (:, :), save :: sr, sr2 !Space's Coordinates, sr2 = 2*sr
    double precision, allocatable, dimension (:, :), save :: srw     !Space Width
    integer NOD(3)      !Number of division is each subbox
    integer b_type(3,2) !Boundary type
    double precision mybox_width(3) !My Box dimension
        
    !!!!!!!!!!!!!!!!
    !Timing Data
    !!!!!!!!!!!!!!!!
    double precision init_time, final_time
    double precision operation_time, communication_time, time_runtime
    double precision time_remove_shadow, time_redistribute_atoms
    double precision time_nlist_renew
    double precision time_predict, time_EAM, time_correct
    double precision time_boundary_folding
    double precision time_temperature_control
    double precision time_collect_atoms
    double precision time_remove_atom, time_shadow_unfold
    real time_sync_start, time_sync_end
  
    !!!!!!!!!!!!!!!!
    !Size of memories
    !!!!!!!!!!!!!!!!
    integer size_of_real, size_of_int
    integer size_of_shadow, size_of_atom
  
    !!!!!!!!!!!!!!!!
    !Qsort Parameter
    !!!!!!!!!!!!!!!!
    integer Short_Array_Limit
  
    !!!!!!!!!!!!!!!!
    !Debug Mode Parameter
    !!!!!!!!!!!!!!!
    integer stage_index
  
    !!!!!!!!!!!!!!!!!
    !Hard Particle Mode
    !!!!!!!!!!!!!!!!!
#ifdef HARD_PARTICLE_MODE
    integer HPM_NOPPA(HPM_NOP) !number of pinned particle atom, in particle 1 and 2
    double precision HPM_force(3, HPM_NOP), HPM_torque(HPM_NOP)
    double precision HPM_cm(3, HPM_NOP), HPM_mass(HPM_NOP), HPM_inertia(HPM_NOP)
    double precision HPM_dr0(3,HPM_NOP), HPM_dtheta(HPM_NOP)
#endif
  
end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Quick Sort Atom Function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module qsort_library
    use data_structure
    contains
    !Quick Sort Atoms
    recursive subroutine Qsort_Atoms(array_start, array_end, level)
      integer, intent(IN) :: array_start, array_end, level
      integer :: split, my_level
      my_level = level + 1
     
      if(array_end-array_start > Short_Array_Limit) then
         call Partition(array_start, array_end, split)       
         
         call Qsort_Atoms(array_start,split-1, my_level)
         call Qsort_Atoms(split, array_end, my_level)

      else if (array_end - array_start > 0) then

        call Isort_Atoms(array_start, array_end)
        
      end if
     
    end subroutine
    
    !Insert Sorting for Short Array
    subroutine Isort_Atoms(array_start, array_end)
        integer, intent(IN) :: array_start, array_end
        integer swap_counter
        integer i1, i2
        
        swap_counter = 10
        
        do while (swap_counter .ne. 0)
            
            swap_counter = 0
            
            do i1 = array_start, array_end-1
            
                i2 = i1 + 1
            
                if (atom_sort_compare(i1, i2) > 0) then
                    swap_counter = swap_counter + 1
                    call swap_atom(i1, i2)
                end if
            
            end do
        end do 
        
    end subroutine
    
    !Qsort Partition Function
    subroutine Partition(array_start, array_end, marker)
      integer, intent(IN) :: array_start, array_end
      integer, intent(OUT) :: marker
      integer :: left, right, pivot
      

      pivot = Find_pivot(array_start, array_end)
      left = array_start - 1 
      right = array_end + 1
     
      do while (left < right)
         right = right - 1
         do while (atom_sort_compare(right,pivot)> 0)
            right = right-1
         end do
         left = left + 1
         do while (atom_sort_compare(left,pivot) < 0)
            left = left + 1
         end do
         if (left < right) then 
         
            ! Swap Pivot Point
            if (left == pivot) then
                pivot = right
            else if (right == pivot) then
                pivot = left
            end if
            
            ! Swap Atom
            call swap_atom(left, right)
            
         end if
      end do

      if (left == right) then
         marker = left + 1
      else
         marker = left
      end if
     
    end subroutine
    
    !Find Best Pivot
    function Find_Pivot(array_start, array_end)
        integer, intent (IN) :: array_start, array_end
        integer pivot_temp(3), swap_temp
        integer Find_Pivot

        pivot_temp(1) = int((array_end - array_start)/4) + array_start
        pivot_temp(2) = int((array_end - array_start)/2) + array_start
        pivot_temp(3) = int((array_end - array_start)/4d0 *3d0) + array_start
          
        if (atom_sort_compare(pivot_temp(1),pivot_temp(2))> 0) then
            swap_temp = pivot_temp(1)
            pivot_temp(1) = pivot_temp(2)
            pivot_temp(2) = swap_temp
        end if
        
        if (atom_sort_compare(pivot_temp(2),pivot_temp(3))> 0) then
            swap_temp = pivot_temp(2)
            pivot_temp(2) = pivot_temp(3)
            pivot_temp(3) = swap_temp
        
        end if
        
        if (atom_sort_compare(pivot_temp(1),pivot_temp(2))> 0) then
            swap_temp = pivot_temp(1)
            pivot_temp(1) = pivot_temp(2)
            pivot_temp(2) = swap_temp
        end if
    
        Find_Pivot = pivot_temp(2)
        
        return
    
    end function 
           
    !Obselete Partitition
    subroutine Partition_with_limit(array_start, array_end, marker)
      integer, intent(IN) :: array_start, array_end
      integer, intent(OUT) :: marker
      integer :: left, right, pivot

      
     
      !pivot_r0(:) =  r0_chunk(:,int((1+size(r0_chunk,2))/2)) ! Average of first and last elements to prevent quadratic 
      pivot = int((array_end - array_start)/2) + array_start
      left = array_start - 1                         ! behavior with sorted or reverse sorted data
      right = array_end + 1
     
      do while ((left < right) .and. ((right-array_start) > short_array_limit) .and. ((array_end-left) > short_array_limit))
         right = right - 1
         do while ((atom_sort_compare(right,pivot)> 0) .and. ((right - array_start) > short_array_limit))
            right = right-1
         end do
         left = left + 1
         do while ((atom_sort_compare(left,pivot) < 0) .and. ((array_end - left) > short_array_limit))
            left = left + 1
         end do
         if ((left < right) .and. ((right-array_start) > short_array_limit) .and. ((array_end-left) > short_array_limit)) then
         
            ! Swap Pivot Point
            if (left == pivot) then
                pivot = right
            else if (right == pivot) then
                pivot = left
            end if
            
            ! Swap Atom
            call swap_atom(left, right)
            
         end if
      end do
      if ((right-array_start) <= short_array_limit) then
         marker = array_start + short_array_limit
      else if ((array_end-left) <= short_array_limit) then
         marker = array_end - short_array_limit  
      else if (left == right) then
         marker = left + 1
      else
         marker = left
      end if
     
    end subroutine
    
    !Rebuild Chain
    subroutine Qsort_Rebuild_Chain()
        integer i1
    
        do i1 = 1, ch_len
            next(i1) = i1+1
            prev(i1) = i1-1
        end do
        
        prev(1) = 0
        next(ch_len) = -1
        ch_start = 1
        ch_end = ch_len
    
        do i1 = 1, ch_empty_count
            ch_empty(i1) = NOS - i1 + 1
        end do
    end subroutine            
    
    !Swap two atoms
    subroutine swap_atom(left, right)
        integer, intent(IN) :: left, right
        double precision rt0(3), rt1(3), rt2(3), rt3(3), rt4(3), rt5(3), drt2_ex(3)
        integer id1, t1, s1, nt1, pt1
    
        rt0(:) = r0(:, right)
        rt1(:) = r1(:, right)
        rt2(:) = r2(:, right)
        rt3(:) = r3(:, right)
        rt4(:) = r4(:, right)
        rt5(:) = r5(:, right) 
        drt2_ex(:) = dr2_ex(:, right)
        id1 = id(right)
        t1  =  t(right)
        s1  =  s(right)
        pt1 = prev(right)
        nt1 = next(right)
    
        r0(:, right) = r0(:, left)
        r1(:, right) = r1(:, left)
        r2(:, right) = r2(:, left)
        r3(:, right) = r3(:, left)
        r4(:, right) = r4(:, left)
        r5(:, right) = r5(:, left)
        dr2_ex(:, right) = dr2_ex(:, right)
        id(right) = id(left)
        t(right) = t(left)
        s(right) = s(left)
        prev(right) = prev(left)
        next(right) = next(left)        
        
        
        r0(:, left) = rt0(:)
        r1(:, left) = rt1(:)
        r2(:, left) = rt2(:)
        r3(:, left) = rt3(:)
        r4(:, left) = rt4(:)
        r5(:, left) = rt5(:)
        dr2_ex(:, left) = drt2_ex(:)
        id(left) = id1
        t(left) = t1
        s(left) = s1
        prev(left) = pt1
        next(left) = nt1        
        
    end subroutine
    
    !Compare Atoms
    function atom_sort_compare(cur, pivot)
        integer, intent (IN) :: cur, pivot
        integer atom_sort_compare
        integer ix(3), ip(3), ix_ip(3), ix_ip_3(3), ix_ip_dir(3)
        integer i_diff(3)
        logical not_chain_cur, not_chain_pivot
        
        not_chain_cur = (s(cur) == -999)
        not_chain_pivot = (s(pivot) == -999)
        
        !if (((next(cur) == -1) .and. (prev(cur) == -1)) .and. ((next(pivot) == -1) .and. (prev(pivot) == -1))) then
        if (cur .eq. pivot) then
        
            atom_sort_compare = 0
        
        else if (not_chain_cur .and. not_chain_pivot) then
        
            atom_sort_compare = 0

        else if (not_chain_cur) then
        
            atom_sort_compare = 1
        
        else if (not_chain_pivot) then
        
            atom_sort_compare = -1
            
        else
       
            ix(:) = qsort_to_int3( (r0(:, cur) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))
            
            ip(:) = qsort_to_int3( (r0(:, pivot) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))      
        
            ix_ip(:) = ix(:) - ip(:)
            
            ix_ip_3(:) = (ix(:)/3) - (ip(:)/3)
            
            !Direction Snake
            ix_ip_dir(1) = (mod(ix(1)/3,2)*2 - 1)*(-1)
            ix_ip_dir(3) = (mod(ix(3)/3,2)*2 - 1)*(-1)
            
            if (ix_ip_3(3) .ne. 0) then
            
                atom_sort_compare = ix_ip_3(3)
            
            else if (ix_ip_3(1).ne. 0) then
            
                atom_sort_compare = ix_ip_3(1) * ix_ip_dir(3)
            
            else if (ix_ip(2) .ne. 0) then
            
                atom_sort_compare = ix_ip(2) * ix_ip_dir(3) * ix_ip_dir(1)
        
            else if (ix_ip(3) .ne. 0) then
                
                atom_sort_compare = ix_ip(3)
                
            else if (ix_ip(1) .ne. 0) then
             
                atom_sort_compare = ix_ip(1)
            
            else 
                i_diff(:) = qsort_to_int3((r0(:, cur) - r0(:, pivot)) * 1000d0)
                
                if (i_diff(3) .ne. 0) then
                
                    atom_sort_compare = i_diff(3)
                
                else if (i_diff(1) .ne. 0) then
                
                    atom_sort_compare = i_diff(1)
                
                else
                
                    atom_sort_compare = i_diff(2)
                    
                end if
            end if
        
        end if
        
        
        return
    end function
    
    !Qsort To Int (3)   
    function qsort_to_int3(rt1)
        integer qsort_to_int3(3)
        double precision rt1(3)
        
        qsort_to_int3(1) = int(rt1(1))
        qsort_to_int3(2) = int(rt1(2))
        qsort_to_int3(3) = int(rt1(3))
        
        return 
    
    end function

    !Dump Sorted Atoms
    subroutine MPI_Dump_Atoms(array_start, array_end, split, tag)
        integer i1, cur
        integer myfile, ix(3)
        integer array_start, array_end, split, tag
        character*100 str_node
        character*200 mystring
        double precision r_p(3)
        double precision mybox_width_ex(3)
        
        myfile = 1500 + myid + array_start + array_end + level
       
        write (str_node, "(I2.2,4I8.8)") myid, array_start, array_end, split, tag
        str_node = "datoms"//str_node(1:34)//".out"
        open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/"//str_node(1:44), STATUS="replace")
        write (mystring, "(9ES14.5)") br(myid, 1, :), nbox_width(:), f_cut_vect_2(:)
        do i1 = array_start, array_end
            ix(:) = qsort_to_int3( (r0(:, i1) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))

            write (mystring, "(3ES14.5,4I8, 3I6)") r0(:, i1), id(i1), prev(i1), next(i1), i1, ix(:)
            write (myfile, '(A)') TRIM(ADJUSTL(mystring))    
            
            
        end do
        close (myfile)
    end subroutine
    
    !Debug Check Sorted Array
    subroutine MPI_Debug_Check_Sort(array_start, array_end, split, level, target_id, target_time)
        integer, intent(IN) :: array_start, array_end, split, level, target_id, target_time
        integer i1, ix1(3), ix2(3)
        
        !if ((myid == target_id))then
        !    print *, "Checking", array_start, array_end, split, level, target_id, target_time
            do i1 = array_start, array_end - 1
                if (atom_sort_compare(i1, i1+1) > 0) then
                    ix1(:) = qsort_to_int3( (r0(:, i1) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))
                    ix2(:) = qsort_to_int3( (r0(:, i1+1) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))
                    print *, "Sort Error (Compare)", myid, i1, ix1(:), ix2(:), atom_sort_compare(i1, i1+1)
                end if
            end do
        !    print *, "Check Done"
            
            
            if ((array_start == 1) .and. (array_end == NOS)) then
                do i1 = array_start, array_end
                
                    if (((i1 <= ch_len) .and. (id(i1) == -999)) .or. ((i1 > ch_len) .and. (id(i1) .ne. -999)))then
                        print *, "Sort Error (Chain)", i1, id(i1), ch_len
                    end if
                end do
            end if
        !end if
        
        
    end subroutine
    
end module

! Neighbor_Maximizing_Sort_Algorithms
module NM_Sort_library
    use data_structure
    contains
    integer NM_nlist_max
    integer NM_nlist[40][200]
    integer NM_nlist_cur[40]
    integer NM_pick_list
    integer NM_nlist_count
    integer NM_corner_id
    integer, allocatable, dimension (:) :: NM_order
    
    subroutine NM_Data_Init()
    
    end subroutine
    
    subroutine Find_Corner()
    
    end subroutine
    
    subroutine NM_assign_order()
    
    end subroutine
    
    subroutine NM_sort_from_order()
    
    subroutine NM_max_nlist()
    
    
end module

!!!!!!!!!!!!!!!
! Module Containing all functions
!!!!!!!!!!!!!!!

module function_library
    !adapt data_structure
    use data_structure
    implicit none
    integer size, ierr, rowtype
    
    contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Data Structure Functions
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
       
    ! Remove atoms in a chain
    subroutine remove_atom(n)
        integer n
        if (n .eq. ch_end) then
            ch_end = prev(n)
            next(prev(n)) = -1
        else if (n .eq. ch_start) then
            ch_start = next(n)
            prev(next(n)) = 0        
        else 
            next(prev(n))=next(n)
            prev(next(n))=prev(n)
        end if
        if (s(n) == 1) NOAA = NOAA - 1
        if (s(n) >= 2) NOPA = NOPA - 1
        next(n) = -1
        prev(n) = -1
        s(n) =  -999
        id(n) = -999
        t(n) =  -999
        ch_empty_count = ch_empty_count + 1
        ch_empty(ch_empty_count) = n
        ch_len = ch_len - 1
        NOA = NOA - 1
        
        
    end subroutine
    
    !Add an atom to the chain
    subroutine add_atom(rt0, rt1, rt2, rt3, rt4, rt5, t1, s1, id1)
        double precision rt0(3), rt1(3), rt2(3), rt3(3), rt4(3), rt5(3)
        integer t1, id1, s1 
        integer n
        
        n = ch_empty(ch_empty_count)
        id(n) = id1  !Copy information of atom to the chain
        r0(:, n) = rt0(:)
        r1(:, n) = rt1(:)
        r2(:, n) = rt2(:)
        r3(:, n) = rt3(:)
        r4(:, n) = rt4(:)
        r5(:, n) = rt5(:)
        t(n) = t1
        s(n) = s1
        id(n) = id1
        Energy(n) = 0
        
        NOA = NOA + 1
        next(n) = -1   !Link the chain
        prev(n) = ch_end
        next(ch_end) = n
        
        if (s1 == 1) NOAA = NOAA + 1
        if (s1 >= 2) NOPA = NOPA + 1
        
        ch_end = n      !Update Chain info
        ch_len = ch_len + 1
        ch_empty_count = ch_empty_count - 1
        
        
    
    end subroutine
     
    !Add an atom to the nbox
    subroutine add_atom_nbox(cur)
        integer cur, ix(3)
        
        ix(:) = to_int3( (r0(:, cur) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))
        
        if (ix(1) < 0) then 
            ix(1) = 0
        end if
        if (ix(1) > NOD(1)) then 
            ix(1) = NOD(1)
        end if
        if (ix(2) < 0) then 
            ix(2) = 0
        end if
        if (ix(2) > NOD(2)) then 
            ix(2) = NOD(2)
        end if
        if (ix(3) < 0) then 
            ix(3) = 0
        end if
        if (ix(3) > NOD(3)) then 
            ix(3) = NOD(3)
        end if


        if ((nbox_time(ix(1), ix(2), ix(3)) == time_index)) then
            nbox_count(ix(1),ix(2),ix(3)) = nbox_count(ix(1),ix(2),ix(3)) + 1
            nbox_next (cur) = nbox_start(ix(1),ix(2),ix(3))
            nbox_start(ix(1),ix(2),ix(3)) = cur
        else 
            nbox_count(ix(1),ix(2),ix(3)) = 1
            nbox_start(ix(1),ix(2),ix(3)) = cur
            nbox_next (cur) = 0
            nbox_time (ix(1), ix(2), ix(3)) = time_index
        end if
    
    end subroutine
       
    !Print atoms to the screen, for Debug puposes only
    subroutine print_list()
        integer i
        
        do i = 1, NOS
            print *, i, id(i), prev(i), next(i), ch_empty(i)

        end do
        print *, ch_empty_count, ch_start, ch_end, ch_len
    end subroutine
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Border Implementation functions
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    function is_my_atom(rt1, myid1)
        logical is_my_atom(2)
        logical in(6)
        double precision rt1(3)
        integer myid1
        
        ! Check is the atom inside my extended area (extended by the force cutoff distance)
        in(1) = (rt1(1) .ge. br(myid1, 1, 1) - f_cut) .and. (rt1(1) .lt. br(myid1, 2, 1) + f_cut)
        in(2) = (rt1(2) .ge. br(myid1, 1, 2) - f_cut) .and. (rt1(2) .lt. br(myid1, 2, 2) + f_cut)
        in(3) = (rt1(3) .ge. br(myid1, 1, 3) - f_cut) .and. (rt1(3) .lt. br(myid1, 2, 3) + f_cut)              
        
        is_my_atom(1) = in(1).and.in(2).and.in(3)
        
        ! Check is the atom inside my area
        in(1) = (rt1(1) .ge. br(myid1, 1, 1)) .and. (rt1(1) .lt. br(myid1, 2, 1))
        in(2) = (rt1(2) .ge. br(myid1, 1, 2)) .and. (rt1(2) .lt. br(myid1, 2, 2))
        in(3) = (rt1(3) .ge. br(myid1, 1, 3)) .and. (rt1(3) .lt. br(myid1, 2, 3))

        is_my_atom(2) = in(1).and.in(2).and.in(3)
        
        return
    
    end function
    
    function in_private_area(rt1)
        implicit none
        logical in_private_area, in(3)
        double precision rt1(3)
        
        in(1) = (rt1(1) .ge. br(myid, 1, 1) + f_cut) .and. (rt1(1) .lt. br(myid, 2, 1) - f_cut)
        in(2) = (rt1(2) .ge. br(myid, 1, 2) + f_cut) .and. (rt1(2) .lt. br(myid, 2, 2) - f_cut)
        in(3) = (rt1(3) .ge. br(myid, 1, 3) + f_cut) .and. (rt1(3) .lt. br(myid, 2, 3) - f_cut)              
       
        in_private_area = in(1).and.in(2).and.in(3)
        return
    end function
    
    function in_domain(rt1)
        implicit none
        logical in_domain, in(3)
        double precision rt1(3)
        
        in(1) = (rt1(1) .ge. br(myid, 1, 1)) .and. (rt1(1) .lt. br(myid, 2, 1))
        in(2) = (rt1(2) .ge. br(myid, 1, 2)) .and. (rt1(2) .lt. br(myid, 2, 2))
        in(3) = (rt1(3) .ge. br(myid, 1, 3)) .and. (rt1(3) .lt. br(myid, 2, 3))              
       
        in_domain = in(1).and.in(2).and.in(3)
        return    
    end function
    
    function in_public_area(rt1)
        implicit none
        logical in_public_area
        double precision rt1(3)
    
        in_public_area = in_domain(rt1) .and. (.not. in_private_area(rt1))
        return    
    end function
    
    function to_int3(rt1)
        integer to_int3(3)
        double precision rt1(3)
        
        to_int3(1) = int(rt1(1))
        to_int3(2) = int(rt1(2))
        to_int3(3) = int(rt1(3))
        
        return 
    
    end function
    
    function fold_boundary(rt0)
        double precision rt0(3)
        double precision fold_boundary(3)
    
        if (sr(1,1) - rt0(1) < 1e-6) then
            rt0(1) = sr(1,1) + abs(sr(1,1) - rt0(1))
        end if
        
        if (sr(2,1) - rt0(1) < 1e-6) then
            rt0(1) = sr(2,1) - abs(sr(2,1) - rt0(1))
        end if
        
        if (sr(1,2) - rt0(2) < 1e-6) then
            rt0(2) = sr(1,2) + abs(sr(1,2) - rt0(2))
        end if
        
        if (sr(2,2) - rt0(2) < 1e-6) then
            rt0(2) = sr(2,2) - abs(sr(2,2) - rt0(2))
        end if
        
        if (sr(1,3) - rt0(3) < 1e-6) then
            rt0(3) = sr(1,3) + abs(sr(1,3) - rt0(3))
        end if
        
        if (sr(2,3) - rt0(3) < 1e-6) then
            rt0(3) = sr(2,3) - abs(sr(2,3) - rt0(3))
        end if
        
        fold_boundary(:) = rt0(:)


        return
    
       
    end function
    
    subroutine update_space_width()
        integer i
        double precision temp(3) /0.0, 0.0, 0.0/
        
        do i = 1, 3
            srw(i, :) = temp(:)
            srw(i, i) = sr(2,i) - sr(1,i)
        end do
    
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Spline Extraction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    function get_rho(t1, rt1)
        integer t1, in1
        double precision rt1, dr1, get_rho
        
        if (rt1 < rho_min(t1)) then
            get_rho = rho(t1, 0, 1)
            return
        else if (rt1 > rho_max(t1)) then
            get_rho = 0d0
            return
        else 
            in1 = int((rt1-rho_min(t1))*rho_inc_inv(t1))
            dr1 = rt1 - (in1 * rho_inc(t1) + rho_min(t1))
            get_rho = rho(t1, in1, 1) + dr1*rho(t1, in1, 2)
            return        
        end if
    end function
    
    function get_drho(t1, rt1)
        integer t1, in1
        double precision rt1, dr1, get_drho
        
        if (rt1 < rho_min(t1)) then
            get_drho = drho(t1, 0, 1)
            return
        else if (rt1 > rho_max(t1)) then
            get_drho = 0d0 !drho(t1, in1, 1) + dr1*drho(t1, in1, 2)
            return
        else 
            in1 = int((rt1-rho_min(t1))*rho_inc_inv(t1))
            dr1 = rt1 - (in1 * rho_inc(t1) + rho_min(t1))
            get_drho = drho(t1, in1, 1) + dr1*drho(t1, in1, 2)
            return        
        end if
    end function
    
    function get_v(ta, tb, rt1)
        integer t1, ta, tb, in1
        double precision rt1, drt1, get_V
        
        if (ta == tb) then
            t1 = ta
        else 
            t1 = 3
        end if
        
        if (rt1 < v_min(t1)) then
            get_v = pot_V(t1, 0, 1)*2
            return
        else if (rt1 > v_max(t1)) then
            get_v = 0d0
            return
        else 
            in1 = int((rt1-v_min(t1))*v_inc_inv(t1))
            drt1 = rt1 - (in1 * v_inc(t1) + v_min(t1))
            get_v = pot_V(t1, in1, 1) + drt1*pot_V(t1, in1, 2)
            return        
        end if
    
    end function
    
    function get_dv(ta, tb, rt1)
        integer t1, ta, tb,  in1
        double precision rt1, drt1, get_dv
        
        if (ta == tb) then
            t1 = ta
        else 
            t1 = 3
        end if
        
        if (rt1 < v_min(t1)) then
            get_dv = for_V(t1, 0, 1)
            !get_dv = 50
            return
        else if (rt1 > v_max(t1)) then
            get_dv = 0d0
            return
        else 
            in1 = int((rt1-v_min(t1))*v_inc_inv(t1))
            drt1 = rt1 - (in1 * v_inc(t1) + v_min(t1))
            get_dv = for_V(t1, in1, 1) + drt1*for_V(t1, in1, 2)
            return        
        end if
    
    end function
    
    function get_f(t1, rho1)
        integer t1, in1
        double precision rho1, drho1, get_f
        
        if (rho1 < 0d0) then
            print *, "error, rho negative", rho1
            get_f = pot_F(t1, 0, 1)
            return
        else if (rho1 > f_max(t1)) then
            get_f = 0d0
            return
        else 
            in1 = int(rho1*f_inc_inv(t1))
            drho1 = rho1 - (in1 * f_inc(t1))
            get_f = pot_F(t1, in1, 1) + drho1*pot_F(t1, in1, 2)
            return        
        end if
    
    end function
    
    function get_df(t1, rho1)
        integer t1, in1
        double precision rho1, drho1, get_df
        
        if (rho1 < 0d0) then
            print *, "error, rho negative", rho1
            get_df = for_F(t1, 0, 1)
            !get_df = 50
            return
        else if (rho1 > f_max(t1)) then
            get_df = 0d0
            return
        else 
            in1 = int((rho1)*f_inc_inv(t1))
            drho1 = rho1 - (in1 * f_inc(t1))
            get_df = for_F(t1, in1, 1) + drho1*for_F(t1, in1, 2)
            return        
        end if
    
    end function
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FILE IO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Read Border File
    subroutine read_border()
        integer i, bindex, auto_divide
        double precision x_width
        double precision srt(2,3)
        double precision sr_offset(3)
        
        open (200,FILE=root_path(1:root_len)//"in/input_border.txt", STATUS="OLD")
        read (200, *) 
        read (200, *) bindex
        read (200, *)
        read (200, *)
        read (200, *) srt(1,1), srt(1,2), srt(1,3), srt(2,1), srt(2,2), srt(2,3)
        read (200, *) 
        read (200, *) 
        read (200, *) sr_offset(1), sr_offset(2), sr_offset(3)
        read (200, *)
        read (200, *)
        read (200, *) b_type(1, 1), b_type(1,2)
        read (200, *) b_type(2, 1), b_type(2,2)
        read (200, *) b_type(3, 1), b_type(3,2)
        read (200, *) 
        read (200, *) 
        read (200, *) auto_divide         
        read (200, *) 
        read (200, *)

        sr(1,:) = srt(1,:) - sr_offset(:)
        sr(2,:) = srt(2,:) + sr_offset(:)        
        
        if (auto_divide == 1) then
            x_width = srt(2,1) - srt(1,1)
            do i = 1, NON-1
                br(i, 1, 1) = dble(i-1)*(x_width/dble(NON-1)) + srt(1,1)
                br(i, 1, 2) = sr(1,2)
                br(i, 1, 3) = sr(1,3)
                br(i, 2, 1) = dble(i)*(x_width/dble(NON-1)) + srt(1,1)
                br(i, 2, 2) = sr(2,2)
                br(i, 2, 3) = sr(2,3)
            end do
            br(1,1,1) = sr(1,1)
            br(NON-1,2,1) = sr(2,1)
        else 
            do i = 1, NON-1
                read(200,*) bindex, br(i,1,1), br(i,1,2), br(i,1,3), br(i,2,1), br(i,2,2), br(i,2,3)
            end do  
        end if

        
        close (200)
    end subroutine
    
    !Read Input Parameters     
    subroutine read_parameters()
        integer NON1
        open (UNIT=300,FILE=root_path(1:root_len)//"in/input_values.txt", STATUS="old")
        read (300, *)
        read (300, *) DIM
        read (300, *) 
        read (300, *) 
        read (300, *) NON1
        read (300, *) 
        read (300, *) 
        read (300, *) NOS
        read (300, *) 
        read (300, *) 
        read (300, *) NOSEA
        read (300, *) 
        read (300, *) 
        read (300, *) NOREA
        read (300, *)
        read (300, *) 
        read (300, *) f_cut
        read (300, *) 
        read (300, *)
        read (300, *) T_wanted
        read (300, *) 
        read (300, *)
        read (300, *) dt
        read (300, *)
        read (300, *)
        read (300, *) end_time
        read (300, *)
        read (300, *)
        read (300, *) movie_time
        read (300, *)
        read (300, *)
        read (300, *) sync_time
        close (300)
        

    end subroutine

    !Read Input Atoms
    subroutine read_atoms()
        integer i, id1, t1, s1
        logical should_add(2)
        double precision rt0(DIM), rt1(DIM), rt2(DIM), rt3(DIM), rt4(DIM),rt5(DIM)
        open (UNIT=100,FILE=root_path(1:root_len)//"in/input_atom.txt", STATUS="OLD")        
        read(100, *) NOTA
       
        do i = 1, NOTA
            read(100, *) id1, t1, s1, rt0(1), rt0(2), rt0(3), rt1(1), rt1(2), rt1(3), rt2(1), rt2(2), rt2(3), &
                        &rt3(1), rt3(2), rt3(3), rt4(1), rt4(2), rt4(3), rt5(1), rt5(2), rt5(3)
            should_add = is_my_atom(rt0, myid)
            if (should_add(1)) then
                if (should_add(2)) then 
                    s1 = 1
                else 
                    s1 = 0
                end if
                call add_atom(rt0, rt1, rt2, rt3, rt4, rt5, t1, s1, id1)
            end if
        end do
        close (100)
    end subroutine
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! MPI Implemented Subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    !Initialize MPI and Data Structure
    subroutine MPI_MD_init()
        integer i, member(0:100)
        integer mygpid1, myid1
        
        
        
        call MPI_INIT( ierr ) 
        call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) 
        call MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr ) 
        call MPI_COMM_GROUP( MPI_COMM_WORLD, MPI_MD_GROUP, ierr ) 
        
        
                
        allocate (status(MPI_STATUS_SIZE))

        allocate (gpid(size), wdid(0:size-2))
        
        !Create New World with Root Excluded        
        do i=1, size
            member(i-1)=i
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call MPI_GROUP_INCL(MPI_MD_GROUP, size - 1, member, MPI_MD_NEW_GROUP, ierr)
        call MPI_COMM_CREATE(MPI_COMM_WORLD, MPI_MD_NEW_GROUP, MPI_MD_WORLD, ierr)
        if (myid .ne. 0) then
            call MPI_COMM_RANK( MPI_MD_WORLD, mygpid, ierr )
            do i =0, size-2
                if (mygpid == i)then
                    call MPI_Bcast(myid, 1, MPI_INTEGER, i, MPI_MD_WORLD, ierr)
                    call MPI_Bcast(mygpid, 1, MPI_INTEGER, i, MPI_MD_WORLD, ierr)
                    gpid(myid) = mygpid
                    wdid(mygpid) = myid
                else
                    call MPI_Bcast(myid1, 1, MPI_INTEGER, i, MPI_MD_WORLD, ierr)
                    call MPI_Bcast(mygpid1, 1, MPI_INTEGER, i, MPI_MD_WORLD, ierr)
                    gpid(myid1) = mygpid1
                    wdid(mygpid1) = myid1
                end if

            end do
            print *, mygpid, " in ", myid, " has started"
        else 
        print *, myid, " the root has started"
        end if
        
        if (myid .eq. 0) then
            allocate (job_id(10000))
            do i = 1, 1000
                job_id(i) = .false.
            end do
        end if
        cmd_quit = .false.
        my_job_force_mode = 1
        batch_mode = .true.
        NON = size
        
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

    end subroutine
    
    ! Initiate Starting Instruction
    function MPI_MD_designate_instruction()
        logical MPI_MD_designate_instruction
        integer i, j, tmp_job_id, tmp_job_index, tmp_job_version, tmp_job_force_mode, tmp_n_int, tmp_n_double
        integer tmp_para_int(10)
        double precision tmp_para_double(10)
        integer n_counter
        real time_start_read, time_now
        double precision temp
        integer myfile1
        character*10 str_index, str_version 

        myfile1 = 1201
        cmd_new_job = .false.
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        n_counter = 0
        
        if (batch_mode) then
            if (myid .eq. 0) then
         
            
                do while ((.not. cmd_new_job) .and. (.not. cmd_quit))
                    print *, "Opening file"
                    open (UNIT=myfile1,FILE="batch/instruction.txt", STATUS="OLD")
                    print *, "Reading file"
                    do while (.true.)
                        read(myfile1, *, end=1000) tmp_job_id, tmp_job_index, tmp_job_version, tmp_job_force_mode, &
                                                   tmp_n_int, tmp_n_double, &
                                                  (tmp_para_int(j),j=1,tmp_n_int),(tmp_para_double(j),j=1,tmp_n_double)
                        
                        if (.not. job_id(tmp_job_id)) then
                            if (.not. cmd_new_job) then 
                                my_job_id = tmp_job_id
                                my_job_index = tmp_job_index
                                my_job_version = tmp_job_version
                                my_job_force_mode = tmp_job_force_mode
                                n_para_int = tmp_n_int
                                n_para_double = tmp_n_double
                                
                                if (n_para_int .ne. 0) then
                                    allocate (para_int(n_para_int))                                    
                                    do j = 1, n_para_int
                                        para_int(j) = tmp_para_int(j)
                                    end do
                                end if
                                
                                if (n_para_double .ne. 0) then
                                    allocate (para_double(n_para_double))
                                    do j = 1, n_para_double
                                        para_double(j) = tmp_para_double(j)
                                    end do
                                end if
                                
                                
                                job_id(my_job_id) = .true.
                                
                                cmd_new_job = .true.
                                
                                call cpu_time(time_sync_start)
                                
                                if (my_job_id == 9999) then
                                    cmd_quit = .true.
                                end if
                            end if
                        end if
                        
                    
                    end do
                    1000 continue
                    close (myfile1)
                    
                    ! Wait 5 minutes
                    if ((.not. cmd_new_job) .and. ( mod (n_counter, 60) .eq. 0))then
                        print *, "No New Instruction Found"
                        Print *, "Wait 5 Minutes"
                    end if
                    n_counter = n_counter + 1

                    time_start_read = secnds(0.0)
                    do while ((secnds(time_start_read) < 300) .and. (.not. cmd_new_job) .and. (.not. cmd_quit))

                        do i = 1, 32760
                            temp = sin (2.0 * PI * cos(2.0 * PI * temp))
                        end do
                    end do        
                end do
            end if
            
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_Bcast(my_job_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(my_job_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(my_job_version, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(my_job_force_mode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(n_para_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(n_para_double, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            if (myid .ne. 0) then
                if (n_para_int .ne. 0) then
                    allocate (para_int(n_para_int))
                    call MPI_Bcast(para_int(1), n_para_int, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                    
                end if
            
                if (n_para_double .ne. 0) then
                    allocate (para_double(n_para_double))
                    call MPI_Bcast(para_double(1), n_para_double, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                end if
            else 
                if (n_para_int .ne. 0) then
                    call MPI_Bcast(para_int(1), n_para_int, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                     
                end if
            
                if (n_para_double .ne. 0) then
                    call MPI_Bcast(para_double(1), n_para_double, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                end if            
            
            end if                
            
            
            
            
            if (my_job_id == 9999) then
                cmd_quit = .true.
            else
                write (str_index, "(I4.4)") my_job_index
                write (str_version, "(I1.1)") my_job_version
#ifdef SharcNet
                root_path = "../MD/BATCH.MD.MPI."//str_index(1:4)//"v" // str_version(1:1)//"/"
                root_len = 26
#else
                root_path = "../BATCH.MD.MPI."//str_index(1:4)//"v" // str_version(1:1)//"/"
                root_len = 23
#endif                
            end if
            MPI_MD_designate_instruction = (.not. cmd_quit)
            
        else
            root_path = ""
            root_len = 0
            MPI_MD_designate_instruction = (.not. cmd_quit)
            cmd_quit = .true.
            
        end if
        
        if (myid .eq. 0) then
            print *, "My Job", my_job_id, "has started"
            if (my_job_id .ne. 1) then
                close (myfile_scr01)
            end if
            
            if (.not. cmd_quit) then
                myfile_scr01 = 35000 + my_job_id
                print *, "My Root Path = " // root_path
                print *, "My Parameter (INT)   : ", n_para_int, (para_int(j), j=1,n_para_int)
                print *, "My Parameter (Double): ", n_para_double, (para_double(j), j=1,n_para_double)
                open (UNIT=myfile_scr01,FILE=root_path(1:root_len)//"status_log.txt", STATUS="REPLACE")
                write (myfile_scr01, *) "My Job ID", my_job_id
                write (myfile_scr01, *) "My Job Index", my_job_index
                write (myfile_scr01, *) "My Job Version", my_job_version
                write (myfile_scr01, *) "My Force Mode",  my_job_force_mode
                write (myfile_scr01, *) "My Parameter (INT)   : ", n_para_int, (para_int(j), j=1,n_para_int)
                write (myfile_scr01, *) "My Parameter (Double): ", n_para_double, (para_double(j), j=1,n_para_double)
            end if
            
        end if

        return
        
    end function
    
    ! Distribute Parameters
    subroutine MPI_MD_distribute_parameters()      
                
        !Read Parameters
        if (myid .eq. 0) then 
            call read_parameters()
        end if
        
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call MPI_Bcast(DIM, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(NON, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(NOS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(NOSEA, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(NOREA, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        
        call MPI_Bcast(f_cut, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)  !! kind change caution
        call MPI_Bcast(T_wanted, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)  !! kind change caution
        call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)   !! kind change caution
        call MPI_Bcast(end_time, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(movie_time, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(sync_time, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    end subroutine
    
    ! Initiate MPI routines and Data structures
    subroutine MPI_MD_data_init()
        integer i, ierr_list(17)
        integer array_int(2)
        double precision array_real(2)
        
        ! Determind Platform determinded variable
        size_of_int = LOC(array_int(2)) - LOC(array_int(1))
        size_of_real = LOC(array_real(2)) - LOC(array_real(1))
        size_of_shadow = 3*size_of_real + size_of_int
        size_of_atom = 6*3*size_of_real + 3*size_of_int
        
        ! Allocate memories        
                                            !Obseleted Code: allocate (r(DIM, NOS), v(DIM, NOS), a(DIM, NOS))
        allocate (r0(DIM, NOS), r1(DIM, NOS), r2(DIM, NOS), r3(DIM, NOS), r4(DIM, NOS), r5(DIM, NOS), r1_ave(DIM),STAT=ierr_list(1))
        allocate (Force(DIM, NOS), STAT=ierr_list(2))
        allocate (Energy(NOS), Stress(DIM, DIM, NOS),Stress_ave(DIM, DIM, NOS), STAT=ierr_list(3))
        allocate (dr2(DIM, NOS), STAT=ierr_list(4))
        allocate (dr2_ex(DIM, NOS), STAT=ierr_list(4))
        allocate (dr0_ex(DIM, NOS), STAT=ierr_list(4))
        allocate (t(NOS), s(NOS), s1_list(NOS), s1_list_inv(NOS), EAM_time_stamp(NOS), id(NOS), STAT=ierr_list(5))
        allocate (next(NOS), prev(NOS), ch_list(NOS), STAT=ierr_list(6))
        allocate (ch_empty(NOS), STAT=ierr_list(7))
        allocate (block_size3(NOSEA), block_size1(NOSEA), block_index(NOSEA), STAT=ierr_list(8))
        allocate (rshadow_r0(DIM, NOREA), STAT=ierr_list(9))
        allocate (ratom_r0(DIM, NOREA), ratom_r1(DIM, NOREA), ratom_r2(DIM, NOREA), STAT=ierr_list(10))
        allocate (ratom_r3(DIM, NOREA), ratom_r4(DIM, NOREA), ratom_r5(DIM, NOREA), STAT=ierr_list(11))
        allocate (rshadow_t(NOREA),rshadow_id(NOREA),ratom_t(NOREA),ratom_s(NOREA),ratom_id(NOREA),STAT=ierr_list(12))
        allocate (sshadow_list(NOSEA), satom_list(NOSEA), image_list(NOSEA*6), STAT=ierr_list(13))
        allocate (br(NON, 2, DIM), STAT=ierr_list(14))
        allocate (sr(2, DIM), sr2(2,DIM), STAT=ierr_list(15))
        allocate (srw(3, DIM), STAT=ierr_list(16))
        allocate (dFi_sum_xp(3, NOS), dVi_sum_xp(3, NOS), rho_sum_xp(NOS), Vi_sum_xp(NOS), STAT=ierr_list(17))
        allocate (Iv1(NOS), Iv2(NOS), Iv3(NOS), Jv2(NOS), Jv3(NOS), STAT=ierr_list(18))
        do i = 1, 18
            !print *, ierr_list(i)
            if (ierr_list(i) .ne. 0) then
                print *, "Allocation Error", myid, i
            end if
        end do
        
        !Initiate Parameters
        do i = 1, NOS
            id(i) = -999
            t(i) = -999
            s(i) = -999
        end do
        
        
        !Initiate Redistribute Module
        if (NON == 2) then
            should_run_forward = .true.
            should_run_reverse = .true.
        else if (NON == 3) then
            should_run_forward = .true.
            should_run_reverse = .true.
        else
            should_run_forward = .true.
            should_run_reverse = .true.
        end if

        ! Initialize Chain        
        do i = 1, NOS
            prev(i) = -1
            next(i) = -1
            ch_empty(i) = NOS - i + 1
            EAM_time_stamp(i)=-1
        end do
        
        do i = 1, NOSEA
            block_size3(i) = 3
            block_size1(i) = 1
        end do
        ch_empty_count = NOS
        ch_start = 1
        ch_end = 0
        ch_len = 0
        next(1) = 0
        NOA = 0
        NOAA = 0
        NOPA = 0
        
        ! Null Time Variable
        time_index = 0
        
        operation_time = 0d0
        communication_time = 0d0
        time_runtime = 0d0
        
        time_remove_shadow = 0d0
        time_redistribute_atoms = 0d0
        time_nlist_renew = 0d0
        time_predict = 0d0
        time_EAM = 0d0
        time_correct = 0d0
        time_boundary_folding = 0d0
        time_temperature_control = 0d0
        time_collect_atoms = 0d0
        time_remove_atom = 0d0
        time_shadow_unfold = 0d0
        
        
        !F_cut Initialization
        f_cut_2 = f_cut*2
        f_cut_sq = f_cut**2
        f_cut_vect(1) = f_cut
        f_cut_vect(2) = f_cut
        f_cut_vect(3) = f_cut
        f_cut_vect_2(1) = f_cut_2
        f_cut_vect_2(2) = f_cut_2
        f_cut_vect_2(3) = f_cut_2
        
        !Vector Class Initialization
        vect_one(1) = 1d0
        vect_one(2) = 1d0
        vect_one(3) = 1d0
        
        !Time Variable
        dt2 = dt**2
        dt3 = dt**3
        dt4 = dt**4
        dt5 = dt**5
               
        !Constants for Gear Predictor
        dt2d2   = dt2/2d0
        dt3d6   = dt3/6d0
        dt4d24  = dt4/24d0
        dt5d120 = dt5/120d0
        
        !Temperature Variables
        T_ratio_sq = 1d0
        
        !Gear Correction Constants
        a0 = 3d0/20d0     !3/16 if force is dependent on velocity  [Com. Sim. of Liquid, APPENDIX E]
        a1 = 251d0/360d0
        a3 = 11d0/18d0
        a4 = 1d0/6d0
        a5 = 1d0/60d0
        
        !Free the pinned atoms for relaxation
        free_pin = .false.
        
        !Set Frame_counter 
        frame_counter = 0
        
        !Set Write_Stress
        write_stress = .true.
        
        !Set Stage_index
        stage_index = 0
        
        !Set Short_Array_Limit
        Short_Array_Limit = 20
        
    end subroutine
    
    ! Neighbor List Data Initialization
    subroutine MPI_MD_nlist_init()
        integer i1, i
        integer ierr_list(4)
        time_index = 0
        if (myid .ne. 0 ) then
        
            i1 = 1

            mybox_width(:) = br(myid,2,:) - br(myid,1,:)
            
            NOD(1) = int((mybox_width(1) + 1.001*f_cut)/f_cut)+3
            NOD(2) = int((mybox_width(2) + 1.001*f_cut)/f_cut)+3
            
            nbox_width(1) = f_cut
            nbox_width(2) = f_cut
            
            do
                nbox_width(3) = f_cut*i1
                NOD(3) = int((mybox_width(3) + f_cut*2d0 + 1.001*f_cut)/(f_cut*i1)) + 1
                if (NOD(1)*NOD(2)*NOD(3)*4*size_of_int < Nlist_RAM) exit
                i1 = i1 + 1
            end do
            print *, "allocate nlist", NOD(1), NOD(2), NOD(3)
            
            allocate (nbox_next(NOS),STAT=ierr_list(1))
            allocate (nbox_start(0:NOD(1), 0:NOD(2), 0:NOD(3)),STAT=ierr_list(2))
            allocate (nbox_count(0:NOD(1), 0:NOD(2), 0:NOD(3)),STAT=ierr_list(3))
            allocate (nbox_time( 0:NOD(1), 0:NOD(2), 0:NOD(3)),STAT=ierr_list(4))
            
            do i = 1, 4
                if (ierr_list(i) .ne. 0) then
                    print *, "NBOX Allocation Error", myid, i
                end if
            end do
            
            call MPI_MD_nlist_renew()
        
        end if
        
    end subroutine
    
    subroutine MPI_MD_nlist_init_XP()
        integer i1, i, j, k
        integer ierr_list(4)
        time_index = 0
        if (myid .ne. 0 ) then
        
            i1 = 1

            mybox_width(:) = br(myid,2,:) - br(myid,1,:)
            
            NOD(1) = CEILING (mybox_width(1)/f_cut)+4
            NOD(2) = CEILING (mybox_width(2)/f_cut)+4
            NOD(3) = CEILING (mybox_width(3)/f_cut)+4

            nbox_width(1) = f_cut
            nbox_width(2) = f_cut
            nbox_width(3) = f_cut

            print *, "allocate nlist", NOD(1), NOD(2), NOD(3)
            
            allocate (nbox_next(NOS),STAT=ierr_list(1))
            allocate (nbox_start(0:NOD(1), 0:NOD(2), 0:NOD(3)),STAT=ierr_list(2))
            allocate (nbox_count(0:NOD(1), 0:NOD(2), 0:NOD(3)),STAT=ierr_list(3))
            allocate (nbox_time( 0:NOD(1), 0:NOD(2), 0:NOD(3)),STAT=ierr_list(4))
            
            do i = 1, 4
                if (ierr_list(i) .ne. 0) then
                    print *, "NBOX Allocation Error", myid, i
                end if
            end do
            
            do k = 0, NOD(3)
            do j = 0, NOD(2)
            do i = 0, NOD(1)
                nbox_time(i,j,k) = -1
                nbox_count(i,j,k) = 0
                nbox_start(i,j,k) = 0
            end do
            end do
            end do
            
            !call MPI_MD_nlist_renew()
        
        end if
        
    end subroutine
    
    ! Distribute Borders Information
    subroutine MPI_MD_distribute_border()
        if (myid .eq. 0) then
            call read_border()
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call MPI_TYPE_CONTIGUOUS(3, MPI_DOUBLE_PRECISION, rowtype, ierr)  !! kind change caution
        call MPI_TYPE_COMMIT(rowtype, ierr)
        call MPI_Bcast(sr(1,1), 6, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(br(1,1,1), 6*NON, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(b_type(1,1), 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call update_space_width()
        sr2 = 2.0*sr
        
   
    end subroutine
    
    !Read and Distribute EAM
    subroutine MPI_MD_distribute_EAM()
        integer read1, read2, i1, i2, in1
        integer i, ierr_list(6)
        double precision read3, read4
        
        if (myid .eq. 0)then
            open (UNIT=400,FILE=root_path(1:root_len)//"in/input_force.txt", STATUS="OLD")
            read(400, *)
            read(400, *) read1
            read(400, *)
            read(400, *)
            read(400, *) read3, read4
            read(400, *)
            read(400, *)
            read(400, *) NOP_Spline
            
            call MPI_Bcast(NOP_Spline, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            m_d3k(1) = read3 * 400907.375046506076143587961766d0
            m_d3k(2) = read4 * 400907.375046506076143587961766d0
            m_atom(1) = read3 * 1.6605387831627258978902255305417d-24   !atomic mass / mole
            m_atom(2) = read4 * 1.6605387831627258978902255305417d-24
            m_dq(1) =  read3 * 103.64269009187202125602479083619d0
            m_dq(2) =  read4 * 103.64269009187202125602479083619d0
            q_dm(1) = 1d0/m_dq(1)
            q_dm(2) = 1d0/m_dq(2)
            
            allocate(pot_V(3, 0:NOP_Spline-1, 2), for_V(3, 0:NOP_Spline-1, 2),STAT=ierr_list(1))
            allocate(pot_F(2, 0:NOP_Spline-1, 2), for_F(2, 0:NOP_Spline-1, 2),STAT=ierr_list(2))
            allocate(rho(2, 0:NOP_Spline-1, 2), drho(2, 0:NOP_Spline-1, 2),STAT=ierr_list(3))
            allocate(v_min(3), v_max(3), f_max(2), rho_max(2), rho_min(2),STAT=ierr_list(4))
            allocate(v_inc(3), f_inc(2), rho_inc(2),STAT=ierr_list(5))
            allocate(v_inc_inv(3), f_inc_inv(2), rho_inc_inv(2),STAT=ierr_list(6))
            
            do i = 1, 6
                if (ierr_list(i) .ne. 0) then
                    print *, "Spline Allocation Error", myid, i
                end if
            end do
            
            
            do i1 = 1, 3
                read(400, *)
                read(400, *)
                read(400, *)
                
                read(400, *) read1, read2
                if (read1 == read2) then
                    in1 = read1
                else 
                    in1 = 3
                end if
                
                read(400, *)
                read(400, *)
                
                read(400, *) v_max(in1), v_min(in1)
                v_inc(in1) = ( v_max(in1) - v_min(in1))/dble(NOP_Spline)
                v_inc_inv(in1) = 1d0/(v_inc(in1))
                                
                read(400, *)
                read(400, *)
                
                do i2 = 0, NOP_Spline-1
                    read(400, *) read3, pot_V(in1, i2, 1), pot_V(in1, i2, 2)
                end do
                
                read(400, *)
                read(400, *)
                
                do i2 = 0, NOP_Spline-1
                    read(400, *) read3, for_V(in1, i2, 1), for_V(in1, i2, 2)
                end do
            end do

            do i1 = 1, 2
                read(400, *)
                read(400, *)
                read(400, *)
                
                read(400, *) read1
                in1 = read1
                
                read(400, *)
                read(400, *)
                
                read(400, *) f_max(in1)
                f_inc(in1) = f_max(in1)/dble(NOP_Spline)
                f_inc_inv(in1) = 1d0/(f_inc(in1))
                
                read(400, *)
                read(400, *)
                
                do i2 = 0, NOP_Spline-1
                    read(400, *) read3, pot_F(in1, i2, 1), pot_F(in1, i2, 2)
                end do
                
                read(400, *)
                read(400, *)
                
                do i2 = 0, NOP_Spline-1
                    read(400, *) read3, for_F(in1, i2, 1), for_F(in1, i2, 2)
                end do
            
            end do
            
            do i1 = 1, 2
                read(400, *)
                read(400, *)
                read(400, *)
                
                read(400, *) read1
                in1 = read1
                
                read(400, *)
                read(400, *)
                
                read(400, *) rho_max(in1), rho_min(in1)
                rho_inc(in1) =  (rho_max(in1)-rho_min(in1))/dble(NOP_Spline)
                rho_inc_inv(in1) = 1d0/(rho_inc(in1))
                
                read(400, *)
                read(400, *)
                
                do i2 = 0, NOP_Spline-1
                    read(400, *) read3, rho(in1, i2, 1), rho(in1, i2, 2)
                end do
                                
                read(400, *)
                read(400, *)
                
                do i2 = 0, NOP_Spline-1
                    read(400, *) read3, drho(in1, i2, 1), drho(in1, i2, 2)
                end do
           
            end do
            call MPI_Bcast(m_d3k, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(m_atom, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)            
            call MPI_Bcast(m_dq, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(q_dm, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(v_max, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(v_min, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(f_max, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(rho_max, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(rho_min, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(pot_V, int(3*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(for_V, int(3*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(pot_F, int(2*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(for_F, int(2*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

            call MPI_Bcast(rho, int(2*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(drho, int(2*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(v_inc, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(f_inc, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(rho_inc, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(v_inc_inv, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(f_inc_inv, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(rho_inc_inv, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            
            deallocate(pot_V, for_V,STAT=ierr_list(1))
            deallocate(pot_F, for_F,STAT=ierr_list(2))
            deallocate(rho, drho,STAT=ierr_list(3))
            deallocate(v_min, v_max, f_max, rho_max, rho_min,STAT=ierr_list(4))
            deallocate(v_inc, f_inc, rho_inc,STAT=ierr_list(5))
            deallocate(v_inc_inv, f_inc_inv, rho_inc_inv,STAT=ierr_list(6))
            
            do i = 1, 6
                if (ierr_list(i) .ne. 0) then
                    print *, "My Root Spline Deallocation Error", i
                end if
            end do
            close (400)

        else 
            
            call MPI_Bcast(NOP_Spline, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            allocate(pot_V(3, 0:NOP_Spline-1, 2), for_V(3, 0:NOP_Spline-1, 2), STAT=ierr_list(1))
            allocate(pot_F(2, 0:NOP_Spline-1, 2), for_F(2, 0:NOP_Spline-1, 2), STAT=ierr_list(2))
            allocate(rho(2, 0:NOP_Spline-1, 2), drho(2, 0:NOP_Spline-1, 2), STAT=ierr_list(3))
            allocate(v_min(3), v_max(3), f_max(2), rho_max(2), rho_min(2), STAT=ierr_list(4))
            allocate(v_inc(3), f_inc(2), rho_inc(2), STAT=ierr_list(5))
            allocate(v_inc_inv(3), f_inc_inv(2), rho_inc_inv(2), STAT=ierr_list(6))

            do i = 1, 6
                if (ierr_list(i) .ne. 0) then
                    print *, "Spline Allocation Error", myid, i
                end if
            end do
            
            
            
            call MPI_Bcast(m_d3k, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(m_atom, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(m_dq, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(q_dm, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(v_max, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(v_min, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(f_max, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(rho_max, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(rho_min, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
                        
            call MPI_Bcast(pot_V, int(3*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(for_V, int(3*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(pot_F, int(2*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(for_F, int(2*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

            call MPI_Bcast(rho, int(2*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(drho, int(2*NOP_Spline*2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            
            call MPI_Bcast(v_inc, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(f_inc, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(rho_inc, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_Bcast(v_inc_inv, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(f_inc_inv, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(rho_inc_inv, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
            f_cut_sq_max(0)= f_cut_sq
            f_cut_sq_max(1)= f_cut_sq
            f_cut_sq_max(2)= dmax1(v_max(1),rho_max(1))**2
            f_cut_sq_max(3)= dmax1(v_max(1),v_max(2), v_max(3),rho_max(1),rho_max(2), rho_max(3))**2
            f_cut_sq_max(4)= dmax1(v_max(2),rho_max(2))**2
    
            !do i1 = 0, 4
            !    print *, f_cut_sq_max(i1)
            !end do
    
        end if
    end subroutine   
       
    ! Distribute Atoms
    subroutine MPI_MD_distribute_atoms()
        integer i, id1, t1, s1, s2, edge_count
        logical should_add(2)
        double precision rt0(DIM), rt1(DIM), rt2(DIM), rt3(DIM), rt4(DIM), rt5(DIM)
        
        NOTAA = 0
        edge_count = 0
        
        if (myid .eq. 0)then
            open (UNIT=110,FILE=root_path(1:root_len)//"in/input_atom.txt", STATUS="OLD")
            read(110, *) NOTA
        end if
        call MPI_Bcast(NOTA, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        !print *, NOTA
       
        do i = 1, NOTA

            if (myid .eq. 0)then 
                read(110, *) id1, t1, s1, rt0(1), rt0(2), rt0(3), rt1(1), rt1(2), rt1(3), rt2(1), rt2(2), rt2(3), &
                            & rt3(1), rt3(2), rt3(3), rt4(1), rt4(2), rt4(3), rt5(1), rt5(2), rt5(3)

                rt0(:) = fold_boundary(rt0(:))
                if (s1 == 1) NOTAA = NOTAA + 1
 
                call MPI_Bcast(rt0(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt1(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt2(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt3(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt4(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt5(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(id1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call MPI_Bcast(t1,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call MPI_Bcast(s1,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
 
            else 
                call MPI_Bcast(rt0(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt1(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt2(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt3(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt4(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(rt5(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)!! kind change caution
                call MPI_Bcast(id1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call MPI_Bcast(t1,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                call MPI_Bcast(s1,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                !print *,"Myid", myid, id1, t1, s1, rt0(1), rt0(2), rt0(3), rt1(1), rt1(2), rt1(3), rt2(1), rt2(2), rt2(3), &
                !                & rt3(1), rt3(2), rt3(3), rt4(1), rt4(2), rt4(3), rt5(1), rt5(2), rt5(3)
                            
                should_add = is_my_atom(rt0, myid)
                !if ((id1 == 485) .and. (myid .ne. 0)) print *, "Receiving 485", myid

                if (should_add(2)) then

                    call add_atom(rt0, rt1, rt2, rt3, rt4, rt5, t1, s1, id1)

                end if
            
            end if
            
        end do
        
        call MPI_Bcast(NOTAA,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        
        if (myid .eq. 0) then
            close (110)
        else
            print *, "Edge_atoms ", edge_count, "Myid ", myid
        end if
        
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine
 

#ifdef HARD_PARTICLE_MODE
    subroutine MPI_MD_hard_particle_init()
        integer i1, i2, cur
        double precision sum_rm(3, HPM_NOP), sum_m(HPM_NOP), sum_mr2(HPM_NOP)
        double precision sum_rm_tmp(3, HPM_NOP), sum_m_tmp(HPM_NOP), sum_mr2_tmp(HPM_NOP)

        cur = ch_start
        
        do i2 = 1, HPM_NOP
            sum_rm_tmp(:, i2) = 0d0*vect_one(:)
            sum_mr2_tmp(i2) = 0d0
            sum_m_tmp(i2) = 0d0
        end do
        
        
        do i1 = 1, ch_len
            do i2 = 1, HPM_NOP
                if (s(cur) == 200 + i2) then
                    sum_rm_tmp(:, i2) = sum_rm_tmp(:, i2)  + r0(:, cur)*m_atom(t(cur))
                    sum_mr2_tmp(i2) = sum_mr2_tmp(i2) + m_atom(t(cur)) *(r0(1, cur)**2 + r0(2, cur)**2) !XY plane
                    sum_m_tmp(i2) = sum_m_tmp(i2) + m_atom(t(cur))                
                end if
            end do
        end do
           
        call MPI_Allreduce(sum_rm(1), sum_rm_tmp(1), 3*HPM_NOP, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MD_WORLD, ierr)
        call MPI_Allreduce(sum_mr2(1), sum_mr2_tmp(1), HPM_NOP, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MD_WORLD, ierr)        
        call MPI_Allreduce(sum_m(1), sum_m_tmp(1),     HPM_NOP, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MD_WORLD, ierr)
        
        
        do i2 = 1, HPM_NOP
            HPM_cm(:, i2) = sum_rm(:, i2)/sum_m(i2)
        end do
        HPM_mass(:) = sum_m(:)
        HPM_inertia(:) = sum_mr2(:) - sum_m(:)*(HPM_cm(1,:)**2 + HPM_cm(2,:)**2)
        
    end subroutine
#endif
 
    ! Sort Atoms into nlist
    subroutine MPI_MD_nlist_renew()
        integer i1, cur
        double precision start_time1, stop_time1, elasped_time1
        start_time1 = MPI_Wtime()
        cur = ch_start
        
                
        s1_count = 0
        
        do i1 = 1, ch_len
            !Force(:, cur) = 0d0*vect_one(:)
            
            call add_atom_nbox(cur)
!            if (s(cur) == 1) then
!                s1_count = s1_count + 1
!                s1_list(s1_count) = cur
!                s1_list_inv(cur) = s1_count
!            else if ((s(cur) == 2) .and. free_pin) then
!                s1_count = s1_count + 1
!                s1_list(s1_count) = cur
!                s1_list_inv(cur) = s1_count
!            end if

            cur = next(cur)
        end do
        
           
        
        stop_time1 = MPI_Wtime()
        elasped_time1 = stop_time1 - start_time1
        time_nlist_renew = time_nlist_renew + elasped_time1
        operation_time = operation_time + elasped_time1
    end subroutine
    
    !Redistribute Atoms for Linear Mode [new]
    subroutine MPI_MD_redist_atoms_Linear_XP()
        integer j, k, si, send_recv1, send_recv2, NON_MD
        integer indexed_real_type, indexed_int_type, sign1
        integer buffer_index
        logical should_add(2)
        double precision cache_zero(3) /0.0, 0.0, 0.0/
        integer rshadow_image_count
        double precision rshadow_image_r(3,8)
        integer send_count(2), recv_count(2), send_size, recv_size
        integer sshadow_size, rshadow_size, satom_size, ratom_size
        integer stotal_count, rtotal_count
        double precision start_time1, stop_time1, elasped_time1
        double precision start_time2, stop_time2, elasped_time2
        character, allocatable, dimension(:), save :: sbuffer, rbuffer
        
        !Debug
        !logical id_cache(11)
        !do i = 1, 11
        !id_cache(i) = .false.
        !end do
        
        start_time1 = MPI_wtime()
        image_count = 0
        NON_MD = NON - 1
        
        !!!!!!!!!!!!!!!!!!!
        ! Send Target
        !!!!!!!!!!!!!!!!!!!
        if (mod(myid,2) == 0) then
            sign1 = 1
            send_recv1 = mod(myid , NON_MD) + 1
            send_recv2 = mod((myid + NON_MD - 2), NON_MD) + 1
        else
            sign1 = -1    
            send_recv1 = mod((myid + NON_MD - 2), NON_MD) + 1
            send_recv2 = mod(myid , NON_MD) + 1
        end if
        

            
        if (should_run_forward) then
        
        call MD_generate_send_list_Linear_XP(sign1)
        
        !!!!!!!!!!!!!!!!!!!!!!
        !   Debug
        !!!!!!!!!!!!!!!!!!!!!
                
        !if (time_index == end_time) then
        !    call MPI_MD_debug_send_data(myid)
        !end if
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Pack and Send Shadow atoms
        !!!!!!!!!!!!!!!!!!!!!!!!!! 
        
        send_count(1) = sshadow_count
        send_count(2) = satom_count
        
        if (b_type(1,1) .ne. 1) then
            if ((send_recv1 == NON_MD) .and. (myid == 1) .and. (NON_MD .ne. 2)) then 
                send_count(1) = 0
                send_count(2) = 0
            end if
        end if
               
        call MPI_sendrecv(send_count,2,MPI_INTEGER,send_recv1,1,recv_count,2,MPI_INTEGER,send_recv1,1,MPI_COMM_WORLD,status,ierr)
        
        rshadow_count = recv_count(1)
        ratom_count = recv_count(2)
        stotal_count = sshadow_count + satom_count
        rtotal_count = rshadow_count + ratom_count


        !!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate Buffer
        !!!!!!!!!!!!!!!!!!!!!! 
        sshadow_size = sshadow_count * size_of_shadow + 1
        satom_size   = satom_count   * size_of_atom   + 1
        rshadow_size = rshadow_count * size_of_shadow + 1
        ratom_size   = ratom_count   * size_of_atom   + 1
        send_size = sshadow_size + satom_size
        recv_size = rshadow_size + ratom_size
            
        allocate (sbuffer(send_size))
        allocate (rbuffer(recv_size))

        if (stotal_count > 0) then

            buffer_index = 1            
        
            
            !!!!!!!!!!!!!!!!
            ! Pack the Shadow
            !!!!!!!!!!!!!!!!

            if (sshadow_count > 0) then
                call MPI_TYPE_indexed (sshadow_count, block_size3, 3*sshadow_list, MPI_DOUBLE_PRECISION, indexed_real_type, ierr)
                call MPI_TYPE_COMMIT(indexed_real_type, ierr)
                            
                call MPI_TYPE_indexed (sshadow_count, block_size1, sshadow_list, MPI_INTEGER, indexed_int_type, ierr)
                call MPI_TYPE_COMMIT(indexed_int_type, ierr)
                
                call MPI_Pack(r0(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(t(1),   1, indexed_int_type,  sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                
                call MPI_TYPE_FREE(indexed_real_type, ierr)
                call MPI_TYPE_FREE(indexed_int_type, ierr)
            end if
            
            !!!!!!!!!!!!!!!!
            ! Pack the Atoms
            !!!!!!!!!!!!!!!
            
            if (satom_count > 0) then
                call MPI_TYPE_indexed (satom_count, block_size3, 3*satom_list, MPI_DOUBLE_PRECISION, indexed_real_type, ierr)
                call MPI_TYPE_COMMIT(indexed_real_type, ierr)
                    
                call MPI_TYPE_indexed (satom_count, block_size1, satom_list, MPI_INTEGER, indexed_int_type, ierr)
                call MPI_TYPE_COMMIT(indexed_int_type, ierr)
                    
                call MPI_Pack(r0(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r1(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r2(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r3(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r4(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r5(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(t(1),   1, indexed_int_type,  sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(s(1),   1, indexed_int_type,  sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(id(1),  1, indexed_int_type,  sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                    
                call MPI_TYPE_FREE(indexed_real_type, ierr)
                call MPI_TYPE_FREE(indexed_int_type, ierr)  
            end if
                     
        end if
          
        !!!!!!!!!!!!!!
        ! Send and Receive
        !!!!!!!!!!!!!

        if ((stotal_count > 0) .and. (rtotal_count > 0)) then
          call MPI_SENDRECV(sbuffer,send_size,MPI_PACKED,send_recv1,2, &
                           &rbuffer,recv_size,MPI_PACKED,send_recv1,2,MPI_COMM_WORLD,status,ierr)
        else if (stotal_count > 0) then
          call MPI_SSEND(sbuffer, send_size, MPI_PACKED, send_recv1, 2, MPI_COMM_WORLD, ierr)  
        else if (rtotal_count > 0) then
          call MPI_RECV(rbuffer(1), recv_size, MPI_PACKED, send_recv1, 2, MPI_COMM_WORLD, status, ierr)  
        end if
        

        !!!!!!!!!!!!!!!!!!!
        ! Unpack Data
        !!!!!!!!!!!!!!!!!!!
        if (rtotal_count > 0) then
            buffer_index = 1
            if (rshadow_count > 0) then
call MPI_Unpack(rbuffer, recv_size, buffer_index, rshadow_r0, rshadow_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, rshadow_t, rshadow_count,   MPI_INTEGER, MPI_COMM_WORLD, ierr)
            end if
            if (ratom_count > 0) then
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r0, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r1, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r2, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r3, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r4, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r5, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)           
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_t, ratom_count,   MPI_INTEGER, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_s, ratom_count,   MPI_INTEGER, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_id, ratom_count,  MPI_INTEGER, MPI_COMM_WORLD, ierr)
            end if
        end if

        deallocate (sbuffer) 
        deallocate (rbuffer)

                 
        start_time2 = MPI_Wtime()
         
        
        ! Atom Adding
        if ((sign1 == 1) .and. (myid == NON_MD)) then
            do j = 1, rshadow_count
                rshadow_r0(1, j) = rshadow_r0(1, j) + srw(1, 1)
                if (rshadow_r0(1,j) > (br(myid, 2, 1) + f_cut)) then
                    call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                    &cache_zero, cache_zero, cache_zero, rshadow_t(j), -1, 0)
                else
                    call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                    &cache_zero, cache_zero, cache_zero, rshadow_t(j), 0, 0)
                end if
            
            end do
            
            do j = 1, ratom_count
                ratom_r0(1, j) = ratom_r0(1, j) + srw(1, 1)
                call add_atom(ratom_r0(:,j), ratom_r1(:,j), ratom_r2(:,j),&
                            &ratom_r3(:,j), ratom_r4(:,j), ratom_r5(:,j), ratom_t(j), ratom_s(j), ratom_id(j))                
            end do 
            
            
        else if ((sign1 == -1) .and. (myid == 1)) then
            do j = 1, rshadow_count
                rshadow_r0(1, j) = rshadow_r0(1, j) - srw(1, 1)
                if (rshadow_r0(1,j) < (br(myid, 1, 1) - f_cut)) then
                    call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                    &cache_zero, cache_zero, cache_zero, rshadow_t(j), -1, 0)
                else
                    call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                    &cache_zero, cache_zero, cache_zero, rshadow_t(j), 0, 0)                
                end if
            end do
            
            do j = 1, ratom_count
                ratom_r0(1, j) = ratom_r0(1, j) - srw(1, 1)
                call add_atom(ratom_r0(:,j), ratom_r1(:,j), ratom_r2(:,j),&
                              &ratom_r3(:,j), ratom_r4(:,j), ratom_r5(:,j), ratom_t(j), ratom_s(j), ratom_id(j))
            end do
            
        else            
            do j = 1, rshadow_count
                if (sign1 == 1) then
                    if (rshadow_r0(1,j) < (br(myid, 2, 1) + f_cut)) then
                        call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, rshadow_t(j), 0, 0)
                    else
                        call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, rshadow_t(j), -1, 0)
                    end if
                else
                    if (rshadow_r0(1,j) > (br(myid, 1, 1) - f_cut)) then
                        call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, rshadow_t(j), 0, 0)
                    else
                        call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, rshadow_t(j), -1, 0)
                    end if
                end if                 
            end do
            
            do j = 1, ratom_count
                call add_atom(ratom_r0(:,j), ratom_r1(:,j), ratom_r2(:,j),&
                             &ratom_r3(:,j), ratom_r4(:,j), ratom_r5(:,j), ratom_t(j), ratom_s(j), ratom_id(j))
            end do  
            
            
                         
        end if
        
        
        
        end if
        
        
        ! Calculate Folding Time
        stop_time2 = MPI_Wtime()
        elasped_time2 = stop_time2 - start_time2
        time_shadow_unfold = time_shadow_unfold + elasped_time2
        operation_time = operation_time + elasped_time2
        communication_time = communication_time - elasped_time2
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Communicate the reverse side
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (should_run_reverse) then
        
        sign1 = -1*sign1
        call MD_generate_send_list_Linear_XP(sign1)
                
!        if (time_index == end_time) then
            !call MPI_MD_debug_send_data(myid+10)
!        end if

        !if ((time_index == 56) .and. (satom_count > 0)) then
        !        print *, "s", j, r(1,satom_list(1)+1), r(2,satom_list(1)+1),r(3,satom_list(1)+1)
        !end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Pack and Send Shadow atoms
        !!!!!!!!!!!!!!!!!!!!!!!!!! 
        send_count(1) = sshadow_count
        send_count(2) = satom_count   
        
        if (b_type(1,2) .ne. 1) then
            if ((send_recv2 == 1) .and. (myid == NON_MD) .and. (NON_MD .ne. 2)) then 
                send_count(1) = 0
                send_count(2) = 0
            end if
        end if
                
        call MPI_sendrecv(send_count,2,MPI_INTEGER,send_recv2,1,recv_count,2,MPI_INTEGER,send_recv2,1,MPI_COMM_WORLD,status,ierr)
        rshadow_count = recv_count(1)
        ratom_count = recv_count(2)
        
        stotal_count = sshadow_count + satom_count
        rtotal_count = rshadow_count + ratom_count

        !!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate Buffer
        !!!!!!!!!!!!!!!!!!!!!!

        sshadow_size = sshadow_count * size_of_shadow + 1
        satom_size   = satom_count   * size_of_atom   + 1
        rshadow_size = rshadow_count * size_of_shadow + 1
        ratom_size   = ratom_count   * size_of_atom   + 1
        send_size = sshadow_size + satom_size
        recv_size = rshadow_size + ratom_size
            
        allocate (sbuffer(send_size))
        allocate (rbuffer(recv_size))

        if (stotal_count > 0) then

            buffer_index = 1            
        
            !!!!!!!!!!!!!!!!
            ! Pack the Shadow
            !!!!!!!!!!!!!!!!

            if (sshadow_count > 0) then
                call MPI_TYPE_indexed (sshadow_count, block_size3, 3*sshadow_list, MPI_DOUBLE_PRECISION, indexed_real_type, ierr)
                call MPI_TYPE_COMMIT(indexed_real_type, ierr)
                            
                call MPI_TYPE_indexed (sshadow_count, block_size1, sshadow_list, MPI_INTEGER, indexed_int_type, ierr)
                call MPI_TYPE_COMMIT(indexed_int_type, ierr)
                
                call MPI_Pack(r0(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(t(1),   1, indexed_int_type,  sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                
                call MPI_TYPE_FREE(indexed_real_type, ierr)
                call MPI_TYPE_FREE(indexed_int_type, ierr)
            end if
            
            !!!!!!!!!!!!!!!!
            ! Pack the Atoms
            !!!!!!!!!!!!!!!
            
            if (satom_count > 0) then
                call MPI_TYPE_indexed (satom_count, block_size3, 3*satom_list, MPI_DOUBLE_PRECISION, indexed_real_type, ierr)
                call MPI_TYPE_COMMIT(indexed_real_type, ierr)
                    
                call MPI_TYPE_indexed (satom_count, block_size1, satom_list, MPI_INTEGER, indexed_int_type, ierr)
                call MPI_TYPE_COMMIT(indexed_int_type, ierr)
                    
                call MPI_Pack(r0(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r1(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r2(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r3(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r4(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(r5(1,1), 1, indexed_real_type, sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(t(1),   1, indexed_int_type,  sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(s(1),   1, indexed_int_type,  sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                call MPI_Pack(id(1),  1, indexed_int_type,  sbuffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
                    
                call MPI_TYPE_FREE(indexed_real_type, ierr)
                call MPI_TYPE_FREE(indexed_int_type, ierr)  
            end if
            
        end if
  
        !!!!!!!!!!!!!!
        ! Send and Receive
        !!!!!!!!!!!!!
        if ((stotal_count > 0) .and. (rtotal_count > 0)) then
          call MPI_SENDRECV(sbuffer,send_size,MPI_PACKED,send_recv2,2,&
                           &rbuffer,recv_size,MPI_PACKED,send_recv2,2,MPI_COMM_WORLD,status,ierr)  
        else if (stotal_count > 0) then     
          call MPI_SSEND(sbuffer, send_size, MPI_PACKED, send_recv2, 2, MPI_COMM_WORLD, ierr)  
        else if (rtotal_count > 0) then
          call MPI_RECV(rbuffer, recv_size, MPI_PACKED, send_recv2, 2, MPI_COMM_WORLD, status, ierr)  
        end if
               
        if (rtotal_count > 0) then
            buffer_index = 1
            if (rshadow_count > 0) then
call MPI_Unpack(rbuffer, recv_size, buffer_index, rshadow_r0, rshadow_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, rshadow_t, rshadow_count,   MPI_INTEGER, MPI_COMM_WORLD, ierr)
            end if
            if (ratom_count > 0) then
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r0, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r1, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r2, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r3, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r4, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_r5, ratom_count*3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_t, ratom_count,   MPI_INTEGER, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_s, ratom_count,   MPI_INTEGER, MPI_COMM_WORLD, ierr)
call MPI_Unpack(rbuffer, recv_size, buffer_index, ratom_id, ratom_count,   MPI_INTEGER, MPI_COMM_WORLD, ierr)
            end if
        end if

        deallocate (sbuffer) 
        deallocate (rbuffer)
        
  !      if (time_index == 56)then
  !          print *, "ab", myid, send_recv2, satom_count, ratom_count
  !      end if
            
                 
        start_time2 = MPI_Wtime()
                                       
        ! Self Shadow Adding
        do j = 1, image_count
            si = image_list(j)
            call add_atom(r0(:,si), cache_zero,cache_zero,&
                            &cache_zero, cache_zero, cache_zero, t(si), 0, 0)
            call remove_atom(si)
            
        end do
        
                                              
        ! Shadow Adding                                      
        do j = 1, rshadow_count
            if (sign1 == 1) then
                if (rshadow_r0(1,j) < (br(myid, 2, 1) + f_cut)) then
                    call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                             &cache_zero, cache_zero, cache_zero, rshadow_t(j), 0, 0)
                else
                    call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                             &cache_zero, cache_zero, cache_zero, rshadow_t(j), -1, 0)
                end if
            else
                if (rshadow_r0(1,j) > (br(myid, 1, 1) - f_cut)) then
                    call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                             &cache_zero, cache_zero, cache_zero, rshadow_t(j), 0, 0)
                else
                    call add_atom(rshadow_r0(:,j),cache_zero,cache_zero,&
                             &cache_zero, cache_zero, cache_zero, rshadow_t(j), -1, 0)
                end if
            end if                 
        end do
        
        do j = 1, ratom_count
            call add_atom(ratom_r0(:,j), ratom_r1(:,j), ratom_r2(:,j),&
                         &ratom_r3(:,j), ratom_r4(:,j), ratom_r5(:,j), ratom_t(j), ratom_s(j), ratom_id(j))
        end do  
        
        end if
        
        !i = ch_start
        !do j = 1, ch_len
        !    if (.not. id_cache(id(i))) then
        
        !    id_cache(id(i)) = .true.
        !    else 
        !        if (id(i) .ne. 0) then
        !        print *, "B", time_index, myid, id(i)
        !        end if
        !    end if
            !print *, time_index, myid, id(i), s(i)
        !    i = next(i)
        !end do
        
        ! Calculate Folding Time
        stop_time2 = MPI_Wtime()
        elasped_time2 = stop_time2 - start_time2
        time_shadow_unfold = time_shadow_unfold + elasped_time2
        operation_time = operation_time + elasped_time2
        communication_time = communication_time - elasped_time2

        !Calculate Communication Time
        stop_time1 = MPI_wtime()
        elasped_time1 = stop_time1 - start_time1
        time_redistribute_atoms = time_redistribute_atoms + elasped_time1
        communication_time = communication_time + elasped_time1        
    end subroutine
    
    !Collect atoms in Binar Mode
    subroutine MPI_MD_collect_atoms_binary()
        integer i, cur
        integer Stripe(0:NON-1)
        integer TNOAA(0:NON-1)
        integer buffer_index
        integer buffer_size, buffer_size_sum
        integer rbuffer_size(0:NON-1)
        integer myfile
        double precision T_ave, r1_ave_tmp(3)
        double precision start_time1, stop_time1, elasped_time1
        character, allocatable, dimension (:) :: buffer
        
        start_time1 = MPI_wtime()
        
        
        if (myid == 0) then
            
            if (time_index == 0)then
                write (myfile_scr01, *) "Output File Initiated", end_time, sync_time, dt
                open (UNIT=50,FILE=root_path(1:root_len)//"out/atoms.out", FORM="unformatted", STATUS="replace")
                write (50) end_time
                write (50) sync_time
                write (50) movie_time
                write (50) dt
                write (50) NOTA
            end if
           
            call MPI_Gather(NOAA + NOPA, 1, MPI_INTEGER, TNOAA(0), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            buffer_size_sum = 0
            Stripe(0) = 0
            T_ave = 0
            
            if (mod(time_index, sync_time*movie_time) == 0) write (50) TNOAA
         
            
            do i = 1, NON-1
                rbuffer_size(i) = buffer_size_sum
                Stripe(i) = TNOAA(i)*size_of_atom + (TNOAA(i)+1)*size_of_real + 2 
                buffer_size_sum = buffer_size_sum + Stripe(i)
                
                !print *, Stripe(i), rbuffer_size(i)
            end do
            call cpu_time(time_sync_end)
            write(myfile_scr01,*)"Synchronized:",time_index," Buffer Size",buffer_size_sum,"Time ",(time_sync_end-time_sync_start)
            write(*,*)"Synchronized:",time_index," Buffer Size",buffer_size_sum,"Time ",(time_sync_end-time_sync_start)
            time_sync_start = time_sync_end
            allocate (buffer(buffer_size_sum))
            
            call MPI_GatherV(i, 0, MPI_INTEGER, buffer, Stripe,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

            
            if (mod(time_index, sync_time*movie_time) == 0) then
                write (myfile_scr01, *) "write ", time_index
#ifdef Write_Atom_Binary
                write (50) Stripe, buffer_size_sum
                write (50) time_index, time_index*dt
                write (50) buffer
#endif                
            end if

            deallocate (buffer)
            
        else

            call MPI_Gather(NOAA + NOPA, 1, MPI_INTEGER, TNOAA, NON, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            r1_ave(:) = 0d0*vect_one(:)           
            buffer_size = (NOAA + NOPA) * size_of_atom + (NOAA + NOPA + 1) * size_of_real + 2
            allocate (buffer(buffer_size))
            buffer_index = 1
            cur = ch_start
            call MPI_Pack(T_current,   1, MPI_DOUBLE_PRECISION,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
            
            do i = 1, ch_len
                if (s(cur) > 0) then

!print *, myid, cur, id(cur), Energy(cur), r0(1, cur), r0(2, cur), r0(3, cur)                    
                    
                    
                    call MPI_Pack(id(cur),   1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Energy(cur),1,MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r0(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r1(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r2(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r3(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r4(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r5(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(t(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(s(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    
                    if (s(cur) == 1) then
                        r1_ave(:) = r1_ave(:) + r1(:,cur)
                    end if
                    
                end if
                cur = next(cur)
            end do
            
            call MPI_GatherV(buffer, buffer_size, MPI_PACKED, buffer, TNOAA,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
            call MPI_Allreduce(r1_ave(1), r1_ave_tmp(1), 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MD_WORLD, ierr)
            r1_ave(:) = r1_ave_tmp(:)/dble(NOTAA)
            
            deallocate (buffer)
        
        end if
        stop_time1 = MPI_wtime()
        elasped_time1 = stop_time1 - start_time1
        time_collect_atoms = time_collect_atoms + elasped_time1
        communication_time = communication_time + elasped_time1  
        
    end subroutine
    
    subroutine MPI_MD_write_frame()
        integer i, j, cur
        integer Stripe(0:NON-1)
        integer TNOAA(0:NON-1)
        integer buffer_index
        integer buffer_size, buffer_size_sum
        integer rbuffer_size(0:NON-1)
        integer myfile
        integer wid, wt, ws        
        double precision wr0(3), wr1(3), wr2(3), wr3(3), wr4(3), wr5(3)
        double precision T_ave, wPE
        character*500 mystring
        character, allocatable, dimension (:) :: buffer
        
        
        if (myid == 0) then
            myfile = 10001

            open (UNIT=myfile,FILE=root_path(1:root_len)//"out/intput_atom.txt", STATUS="replace")

            write (mystring, *) NOTA
            write (myfile, '(A)') TRIM(ADJUSTL(mystring))
            !print *, "b2"

            
            !print *, "b3"
            call MPI_Gather(NOAA+NOPA, 1, MPI_INTEGER, TNOAA(0), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            !print *, "b4"
            buffer_size_sum = 0
            Stripe(0) = 0
            T_ave = 0
            !print *, "b5"
            do i = 1, NON-1
                !print *, "b6", i
                rbuffer_size(i) = buffer_size_sum
                Stripe(i) = TNOAA(i)*size_of_atom + (TNOAA(i)+1)*size_of_real + 2 
                buffer_size_sum = buffer_size_sum + Stripe(i)
            end do
            
            allocate (buffer(buffer_size_sum))
            
            call MPI_GatherV(i, 0, MPI_INTEGER, buffer, Stripe,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
            
            buffer_index = 1

            do i = 1, NON-1
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, T_current, 1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD, ierr)
                T_ave = T_ave + T_current * TNOAA(i)/ dble(NOTA)
                do j = 1, TNOAA(i)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wid, 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wPE, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr0, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr1, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr2, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr3, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr4, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr5, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wt , 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, ws , 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                
                write (mystring, *) wid, wt, ws, wr0(1), wr0(2), wr0(3), wr1(1), wr1(2), wr1(3), &
                                             & wr2(1), wr2(2), wr2(3), wr3(1), wr3(2), wr3(3), &
                                             & wr4(1), wr4(2), wr4(3), wr5(1), wr5(2), wr5(3)
                write (myfile, '(A)') TRIM(ADJUSTL(mystring))
            
                end do
                buffer_index = buffer_index + 2
            end do
            T_current = T_ave
            close (myfile)
            deallocate (buffer)
            
        else
            
            call MPI_Gather(NOAA + NOPA, 1, MPI_INTEGER, TNOAA, NON, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            buffer_size = (NOAA + NOPA) * size_of_atom + (NOAA + NOPA + 1) * size_of_real + 2
            allocate (buffer(buffer_size))
            buffer_index = 1
            cur = ch_start
            call MPI_Pack(T_current,   1, MPI_DOUBLE_PRECISION,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
            do i = 1, ch_len
                if (s(cur) > 0) then
                    
                    call MPI_Pack(id(cur),   1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Energy(cur),1,MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r0(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r1(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r2(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r3(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r4(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r5(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(t(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(s(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    
                end if
                cur = next(cur)
            end do
            
            call MPI_GatherV(buffer, buffer_size, MPI_PACKED, buffer, TNOAA,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
            
            deallocate (buffer)
        
        end if
        
    end subroutine
    
    subroutine MPI_MD_write_cfg()
        integer i, j, cur
        integer Stripe(0:NON-1)
        integer TNOAA(0:NON-1)
        integer buffer_index
        integer buffer_size, buffer_size_sum
        integer rbuffer_size(0:NON-1)
        integer myfile1
        integer wid, wt, ws        
        double precision wr0(3), wr1(3), wr2(3), wr3(3), wr4(3), wr5(3)
        double precision T_ave, wPE
        double precision r_p(3)
        double precision mybox_width_ex(3)
        character*30 str_node
        character*500 mystring
        character, allocatable, dimension (:) :: buffer
        
                
        if (myid == 0) then
            myfile1 = 15000 + frame_counter
            
            write (str_node, "(I4.4)") frame_counter
            str_node = "RunTime_"//str_node(1:4)//".cfg"
            open (UNIT=myfile1,FILE=root_path(1:root_len)//"out/"//str_node(1:16), STATUS="replace")
                      
            
            mybox_width_ex(:) = sr(2,:) - sr(1,:)

            write (mystring, *) "Number of particles = ", NOTA
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "A = 1 Angstorm"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,1) = ", mybox_width_ex(1), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,2) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,3) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,1) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,2) = ", mybox_width_ex(2), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,3) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,1) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,2) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,3) = ", mybox_width_ex(3), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) ".NO_VELOCITY."
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "entry_count = 7"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[0] = PE [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[1] = ID [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[2] = Type [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[3] = Status [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))    
            write (mystring, *) "27.0"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "Al"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))            
            !print *, "b3"
            call MPI_Gather(NOAA+NOPA, 1, MPI_INTEGER, TNOAA(0), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            !print *, "b4"
            buffer_size_sum = 0
            Stripe(0) = 0
            T_ave = 0
            !print *, "b5"
            do i = 1, NON-1
                !print *, "b6", i
                rbuffer_size(i) = buffer_size_sum
                Stripe(i) = TNOAA(i)*(3*size_of_int) + (TNOAA(i)*4 + 1)*size_of_real + 2 
                buffer_size_sum = buffer_size_sum + Stripe(i)
            end do
            
            allocate (buffer(buffer_size_sum))
            
            call MPI_GatherV(i, 0, MPI_INTEGER, buffer, Stripe,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
            
            buffer_index = 1
            do i = 1, NON-1
                write (myfile_scr01, *) "ID :", i, "TNOAA :", TNOAA(i)
            end do

            do i = 1, NON-1
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, T_current, 1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD, ierr)
                T_ave = T_ave + T_current * TNOAA(i)/ dble(NOTA)
                do j = 1, TNOAA(i)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wid, 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wPE, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr0, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                !call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr1, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                !call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr2, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                !call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr3, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                !call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr4, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                !call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr5, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wt , 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, ws , 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                
                r_p(:) = (wr0(:)-sr(1,:))/mybox_width_ex(:)
                write (mystring, "(4ES14.5,I8,2I4)") r_p(:), wPE, wid, wt, ws
                write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            
                end do
                buffer_index = buffer_index + 2
            end do
            T_current = T_ave
            
            deallocate (buffer)
            close (myfile1)
        else
            
            call MPI_Gather(NOAA + NOPA, 1, MPI_INTEGER, TNOAA, NON, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            buffer_size = (NOAA + NOPA) * (3*size_of_int) + (4*(NOAA + NOPA) + 1) * size_of_real + 2
            allocate (buffer(buffer_size))
            buffer_index = 1
            cur = ch_start
            call MPI_Pack(T_current,   1, MPI_DOUBLE_PRECISION,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
            do i = 1, ch_len
                if (s(cur) > 0) then
                    
                    call MPI_Pack(id(cur),   1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Energy(cur),1,MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r0(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    !call MPI_Pack(r1(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    !call MPI_Pack(r2(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    !call MPI_Pack(r3(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    !call MPI_Pack(r4(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    !call MPI_Pack(r5(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(t(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(s(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    
                end if
                cur = next(cur)
            end do
            
            call MPI_GatherV(buffer, buffer_size, MPI_PACKED, buffer, TNOAA,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
            
            deallocate (buffer)
        
        end if
        
        frame_counter = frame_counter + 1
        
    end subroutine

    ! Write CFG with Stress
    subroutine MPI_MD_write_stress_cfg()
        integer i, j, cur
        integer Stripe(0:NON-1)
        integer TNOAA(0:NON-1)
        integer buffer_index
        integer buffer_size, buffer_size_sum
        integer rbuffer_size(0:NON-1)
        integer myfile1  !, myfile2
        integer wid, wt, ws        
        double precision wr0(3), wr1(3), wr2(3), wr3(3), wr4(3), wr5(3)
        double precision T_ave, wPE, wStress(9), wStress_ave(9), wIv1, wIv2, wIv3, wJv2, wJv3
        double precision r_p(3)
        double precision mybox_width_ex(3)
        character*30 str_node
        character*500 mystring
        character, allocatable, dimension (:) :: buffer
        
                
        if (myid == 0) then
            myfile1 = 15000 + frame_counter
            !myfile2 = 25000 + frame_counter
            
            write (str_node, "(I4.4)") frame_counter
            str_node = "Runtime_"//str_node(1:4)//".cfg"
            open (UNIT=myfile1,FILE=root_path(1:root_len)//"out/"//str_node(1:16), STATUS="replace")
            !open (UNIT=myfile2,FILE=root_path(1:root_len)//"out/"//str_node(1:11)//".aux", STATUS="replace")
                      
            
            mybox_width_ex(:) = sr(2,:) - sr(1,:)

            write (mystring, *) "Number of particles = ", NOTA
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "A = 1 Angstorm"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,1) = ", mybox_width_ex(1), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,2) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,3) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,1) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,2) = ", mybox_width_ex(2), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,3) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,1) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,2) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,3) = ", mybox_width_ex(3), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) ".NO_VELOCITY."
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "entry_count = 18"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[0] = PE [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[1] = ID [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[2] = Type [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[3] = Status [reduced unit]"
            
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[4] = I1 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[5] = I2 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[6] = I3 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[7] = J2 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[8] = J3 [reduced unit]"            
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))  
              
            write (mystring, *) "auxiliary[9] = S11 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[10] = S12 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[11] = S13 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))  
            write (mystring, *) "auxiliary[12] = S22 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[13] = S23 [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "auxiliary[14] = S33 [reduced unit]"

            write (myfile1, '(A)') TRIM(ADJUSTL(mystring)) 
            write (mystring, *) "27.0"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "Al"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))            
            !print *, "b3"
            call MPI_Gather(NOAA+NOPA, 1, MPI_INTEGER, TNOAA(0), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            !print *, "b4"
            buffer_size_sum = 0
            Stripe(0) = 0
            T_ave = 0
            !print *, "b5"
            do i = 1, NON-1
                !print *, "b6", i
                rbuffer_size(i) = buffer_size_sum
                Stripe(i) = TNOAA(i)*3*size_of_int + (18*TNOAA(i)+1)*size_of_real + 2 
                buffer_size_sum = buffer_size_sum + Stripe(i)
            end do
            
            allocate (buffer(buffer_size_sum))
            
            call MPI_GatherV(i, 0, MPI_INTEGER, buffer, Stripe,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
            
            buffer_index = 1
            do i = 1, NON-1
                write (myfile_scr01, *) "ID :", i, "TNOAA :", TNOAA(i)
            end do

            do i = 1, NON-1
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, T_current, 1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD, ierr)
                T_ave = T_ave + T_current * TNOAA(i)/ dble(NOTA)
                do j = 1, TNOAA(i)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wid, 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wPE, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr0, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
!                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr1, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
!                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr2, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
!                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr3, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
!                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr4, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
!                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr5, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wt , 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, ws , 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wStress, 9, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wIv1, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wIv2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wIv3, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wJv2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wJv3, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                
                r_p(:) = (wr0(:)-sr(1,:))/mybox_width_ex(:)
                write (mystring, "(4ES14.5,I8,2I4,11ES14.5)") r_p(:), wPE, wid, wt, ws, wIv1, wIv2, wIv3, wJv2, wJv3, &
                & wStress(1:3),wStress(5:6),wStress(9)
                write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
                !write (mystring, "(14ES14.5)") wStress(:), wIv1, wIv2, wIv3, wJv2, wJv3
                !write (myfile2, '(A)') TRIM(ADJUSTL(mystring))
                end do
                buffer_index = buffer_index + 2
            end do
            T_current = T_ave
            
            deallocate (buffer)
            close (myfile1)
            !close (myfile2)
        else
            
            call MPI_Gather(NOAA + NOPA, 1, MPI_INTEGER, TNOAA, NON, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            buffer_size = (NOAA + NOPA) * 3 * size_of_int + ((NOAA + NOPA)* 18 + 1) * size_of_real + 2
            allocate (buffer(buffer_size))
            buffer_index = 1
            cur = ch_start
            call MPI_Pack(T_current,   1, MPI_DOUBLE_PRECISION,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
            do i = 1, ch_len
                if (s(cur) > 0) then
                    
                    call MPI_Pack(id(cur),    1,MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Energy(cur),1,MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r0(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
!                    call MPI_Pack(r1(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
!                    call MPI_Pack(r2(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
!                    call MPI_Pack(r3(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
!                    call MPI_Pack(r4(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
!                    call MPI_Pack(r5(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(t(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(s(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
         call MPI_Pack(Stress_ave(1,1,cur), 9, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Iv1(cur), 1, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Iv2(cur), 1, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Iv3(cur), 1, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Jv2(cur), 1, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Jv3(cur), 1, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    
                end if
                cur = next(cur)
            end do
            
            call MPI_GatherV(buffer, buffer_size, MPI_PACKED, buffer, TNOAA,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
            
            deallocate (buffer)
        
        end if
        
        frame_counter = frame_counter + 1
        
    end subroutine    
        
    subroutine MPI_MD_write_input()
        integer i, j, cur
        integer Stripe(0:NON-1)
        integer TNOAA(0:NON-1)
        integer buffer_index
        integer buffer_size, buffer_size_sum
        integer rbuffer_size(0:NON-1)
        integer myfile
        integer wid, wt, ws        
        double precision wr0(3), wr1(3), wr2(3), wr3(3), wr4(3), wr5(3)
        double precision T_ave, wPE
        double precision r_p(3)
        double precision mybox_width_ex(3)
        character*30 str_node
        character*500 mystring
        character, allocatable, dimension (:) :: buffer
        
        call debug_print_stage(0)
        if (myid == 0) then
            myfile = 20001
            !write (str_node, "(I3.3)") nint(time_index/2500d0) 
            str_node = "input_atom_final.txt"
            open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/"//str_node(1:20), STATUS="replace")

            write (mystring, '(I8)') NOTA
            write (myfile, '(A)') TRIM(ADJUSTL(mystring))   
            
call debug_print_stage(stage_index)!01
            
            call MPI_Gather(NOAA+NOPA, 1, MPI_INTEGER, TNOAA(0), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            buffer_size_sum = 0
            Stripe(0) = 0
            T_ave = 0

            do i = 1, NON-1
                !print *, "b6", i
                rbuffer_size(i) = buffer_size_sum
                Stripe(i) = TNOAA(i)*size_of_atom + (TNOAA(i)+1)*size_of_real + 2 
                buffer_size_sum = buffer_size_sum + Stripe(i)
            end do
            
call debug_print_stage(stage_index)!02
            
            allocate (buffer(buffer_size_sum))
            
call debug_print_stage(stage_index)!03

            call MPI_GatherV(i, 0, MPI_INTEGER, buffer, Stripe,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

            buffer_index = 1
call debug_print_stage(stage_index)!04

            write (myfile_scr01, *) time_index, TNOAA(:)
call debug_print_stage(stage_index)!05
            do i = 1, NON-1
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, T_current, 1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD, ierr)
                T_ave = T_ave + T_current * TNOAA(i)/ dble(NOTA)
                do j = 1, TNOAA(i)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wid, 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wPE, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr0, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr1, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr2, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr3, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr4, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wr5, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, wt , 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, ws , 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                
                write (mystring, "(I8,2I3,18ES14.5)") wid, wt, ws, wr0(:), wr1(:), wr2(:), wr3(:), wr4(:), wr5(:) 
                write (myfile, '(A)') TRIM(ADJUSTL(mystring))
            
                end do
                buffer_index = buffer_index + 2
            end do
            T_current = T_ave
call debug_print_stage(stage_index)!06            
            deallocate (buffer)
            close (myfile)
call debug_print_stage(stage_index)!07         
        else
        
call debug_print_stage(stage_index)!*1
            call MPI_Gather(NOAA + NOPA, 1, MPI_INTEGER, TNOAA, NON, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            buffer_size = (NOAA + NOPA) * size_of_atom + (NOAA + NOPA + 1) * size_of_real + 2
            allocate (buffer(buffer_size))
call debug_print_stage(stage_index)!*2            
            buffer_index = 1
            cur = ch_start
            call MPI_Pack(T_current,   1, MPI_DOUBLE_PRECISION,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
call debug_print_stage(stage_index)!*3
            do i = 1, ch_len
                if (s(cur) > 0) then
                    
                    call MPI_Pack(id(cur),   1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Energy(cur),1,MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r0(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r1(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r2(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r3(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r4(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(r5(1,cur), 3, MPI_DOUBLE_PRECISION, buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(t(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(s(cur),    1, MPI_INTEGER,  buffer, buffer_size, buffer_index, MPI_COMM_WORLD, ierr)
                    
                end if
                cur = next(cur)
            end do
call debug_print_stage(stage_index)!*4
            call MPI_GatherV(buffer, buffer_size, MPI_PACKED, buffer, TNOAA,rbuffer_size, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
            deallocate (buffer)
call debug_print_stage(stage_index)!*5            
        
        end if
        
    end subroutine    
    
    !Record Time
    subroutine MPI_MD_record_time()
        integer i
        CHARACTER(LEN=480) :: FMT1 = "(1A30, 1E20.14)"

        final_time = MPI_Wtime()
        time_runtime = final_time - init_time
        if (myid == 0) then
            open (UNIT=800,FILE=root_path(1:root_len)//"out/time.txt", STATUS="replace")
            write (800, *) "Root Run Time:", time_runtime
            write (800, *) "Root Collect Atom Time:", time_collect_atoms
            
        end if
       
        
        if (myid == 0)then
            do i = 1, NON-1
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_RECV(operation_time, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(communication_time, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(time_runtime, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)

            call MPI_RECV(time_remove_shadow, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)          
            call MPI_RECV(time_redistribute_atoms, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(time_nlist_renew, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)

            call MPI_RECV(time_predict, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(time_EAM, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(time_correct, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)

            call MPI_RECV(time_boundary_folding, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(time_temperature_control, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(time_collect_atoms, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)

            call MPI_RECV(time_remove_atom, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(time_shadow_unfold, 1, MPI_DOUBLE_PRECISION , i, 0, MPI_COMM_WORLD, status, ierr)

            write (800, *) "Myid", i
            write (800, fmt1) "Operation Time:", operation_time
            write (800, fmt1) "Communication Time:", communication_time
            write (800, fmt1) "Total Run Time:", time_runtime

            write (800, fmt1) "Remove Shadow Time:", time_remove_shadow
            write (800, fmt1) "Redistribute Atoms Time:", time_redistribute_atoms - time_shadow_unfold
            write (800, fmt1) "Nlist Renew time:", time_nlist_renew

            write (800, fmt1) "Gear Predict Time:", time_predict
            write (800, fmt1) "EAM Time:", time_EAM
            write (800, fmt1) "Gear Correct Time:", time_correct

            write (800, fmt1) "Boundary Folding Time:", time_boundary_folding
            write (800, fmt1) "Temperature Control Time:", time_temperature_control
            write (800, fmt1) "Collect Atoms Time:", time_collect_atoms

            write (800, fmt1) "Remove Atom Time:", time_remove_atom
            write (800, fmt1) "Shadow Unfold Time:", time_shadow_unfold
            write (800, *)
            end do
            close (800)
        else 
            do i = 1, NON-1
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            if(myid == i) then
            call MPI_SEND(operation_time, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(communication_time, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(time_runtime, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_SEND(time_remove_shadow, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(time_redistribute_atoms, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(time_nlist_renew, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_SEND(time_predict, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(time_EAM, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(time_correct, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_SEND(time_boundary_folding, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(time_temperature_control, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(time_collect_atoms, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            
            call MPI_SEND(time_remove_atom, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(time_shadow_unfold, 1, MPI_DOUBLE_PRECISION , 0, 0, MPI_COMM_WORLD, ierr)
            
            end if
            end do
        end if
        
    end subroutine
   
    !Finalize
    subroutine MPI_MD_finalize()
        integer i, ierr_list(50), i1

        deallocate (r0, r1, r2, r3, r4, r5, r1_ave, STAT=ierr_list(1))       
        deallocate (Force, STAT=ierr_list(2))               
        deallocate (Energy, STAT=ierr_list(3))
        deallocate (dr2, STAT=ierr_list(4))                
        deallocate (dr2_ex, STAT=ierr_list(5))             
        deallocate (dr0_ex, STAT=ierr_list(6))                  
        deallocate (t, s, id, STAT=ierr_list(7))        
        deallocate (s1_list, s1_list_inv, STAT=ierr_list(8))        
        deallocate (EAM_time_stamp, STAT=ierr_list(9))        
        deallocate (next, prev, STAT=ierr_list(10))        
        deallocate (ch_list, STAT=ierr_list(11))                
        deallocate (ch_empty, STAT=ierr_list(12))        
        deallocate (block_size3, block_size1, block_index, STAT=ierr_list(13))        
        deallocate (rshadow_r0, ratom_r0, ratom_r1, ratom_r2, ratom_r3, ratom_r4, ratom_r5, STAT=ierr_list(14))        
        deallocate (rshadow_t, rshadow_id, ratom_t, ratom_s, ratom_id, STAT=ierr_list(15))        
        deallocate (sshadow_list, satom_list, image_list, STAT=ierr_list(16))        
        deallocate (br, STAT=ierr_list(17))        
        deallocate (sr, sr2, STAT=ierr_list(18))        
        deallocate (srw, STAT=ierr_list(19))         
        deallocate (dFi_sum_xp, dVi_sum_xp, STAT=ierr_list(20))                         
        deallocate (rho_sum_xp, Vi_sum_xp, STAT=ierr_list(21))
        if (n_para_int .ne. 0) then
        deallocate (para_int, STAT=ierr_list(22))
        end if
        if (n_para_double .ne. 0) then          
        deallocate (para_double, STAT=ierr_list(23))
        end if
        
        
!        deallocate (act_list) 
        
        
        if (myid == 0) then
            print *, "Closing File"
            
            close (UNIT=50, IOSTAT=ierr_list(23))
            
            do i = 1, 23
                if (ierr_list(i) .ne. 0) then
                    print *, "Deallocation Error", myid, i
                end if
            end do

        else
            deallocate (pot_V, for_V, STAT=ierr_list(22))    
            deallocate (pot_F, for_F, STAT=ierr_list(23))
            deallocate (rho, drho, STAT=ierr_list(24))
            deallocate (v_min, v_max, f_max, rho_min, rho_max, STAT=ierr_list(25))              
            deallocate (v_inc, f_inc, rho_inc, STAT=ierr_list(26))            
            deallocate (v_inc_inv, f_inc_inv, rho_inc_inv, STAT=ierr_list(27))            
            deallocate (nbox_next, STAT=ierr_list(28))            
            deallocate (nbox_start, STAT=ierr_list(29))            
            deallocate (nbox_count, STAT=ierr_list(30))            
            deallocate (nbox_time, STAT=ierr_list(31))

            
            do i = 1, 31
                if (ierr_list(i) .ne. 0) then
                    print *, "Deallocation Error", myid, i
                end if
            end do
            
            
        end if            
        print *, "My ID", myid, "Finalized"
            
    end subroutine
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! MD implementation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Apply external strain and stress                                                       wwjd
    subroutine apply_strain_stress(cur)
            integer cur
            double precision dr_cen(3), dr_norm_sq, dr_norm_inv, sign_x, mag

            
            !Strain field
            !dr0_ex(:, cur) = 0d0*vect_one(:)
            
            if (my_job_force_mode == 1) then            
                if (time_index > 1000) then
                    r0(1, cur) = r0(2,cur)*para_double(1) + r0(1, cur)
                end if
            else if (my_job_force_mode == 2) then
                if (time_index > 1) then
                    r0(2, cur) = (r0(1,cur)-(sr(1,1)+sr(2,1))/2d0)*para_double(1) + r0(2, cur)
                end if
            else if (my_job_force_mode == 3) then
                if (time_index > 1000) then
                    r0(1, cur) = r0(2,cur)*para_double(1) + r0(1, cur)
                    r0(2, cur) = (r0(1,cur)-(sr(1,1)+sr(2,1))/2d0)*para_double(2) + r0(2, cur)
                end if
            else if (my_job_force_mode == 4) then
                if (time_index > 1000) then
                    r0(1, cur) = r0(2,cur)*2d-6 + r0(1, cur)
                    r0(2, cur) = (r0(1,cur)-(sr(1,1)+sr(2,1))/2d0)*(-2d-6) + r0(2, cur)
                end if
            else if (my_job_force_mode == 5) then
                if (time_index > 1000) then
                    if (s(cur) == 10) then
                        r0(1, cur) = (-3d-4) + r0(1, cur)
                    else if (s(cur) == 11) then
                        r0(1, cur) = (3d-4) + r0(1, cur)
                    end if
                end if
            else if (my_job_force_mode == 6) then
                if (time_index > 1000) then
                    if (s(cur) == 10) then
                        !r0(1, cur) = (-6d-4) + r0(1, cur)
                        r0(1, cur) = (-para_double(1)) + r0(1, cur)
                    else if (s(cur) == 11) then
                        !r0(1, cur) = (6d-4) + r0(1, cur)
                        r0(1, cur) = (para_double(1)) + r0(1, cur)
                    end if
                end if
            else if (my_job_force_mode == 7) then
                if (time_index > 1000) then
                    if (s(cur) == 10) then
                        r0(1, cur) = (-12d-4) + r0(1, cur)
                    else if (s(cur) == 11) then
                        r0(1, cur) = (12d-4) + r0(1, cur)
                    end if
                end if                                    
            else if (my_job_force_mode == 8) then
                if (time_index > 1000) then
                    if (s(cur) == 10) then
                        r0(1, cur) = (-36d-4) + r0(1, cur)
                    else if (s(cur) == 11) then
                        r0(1, cur) = (36d-4) + r0(1, cur)
                    end if
                end if
                
            else if (my_job_force_mode == 9) then            
                if (time_index > 1000) then
                    if (s(cur) .ne. 3) then
                        r0(1, cur) = r0(2,cur)*2d-6 + r0(1, cur)
                    end if
                end if  
            else if (my_job_force_mode == 10) then            
                if (time_index > 1000) then
                    if (s(cur) < 200) then
                        r0(1, cur) = r0(2,cur)*4d-6 + r0(1, cur)
                    end if
                end if                  
                
            else if (my_job_force_mode == 11) then            
                if (mod(time_index, 500) == 0 .and. time_index .ne. 0) then
                    T_wanted = T_wanted/2d0
                end if          
            else if (my_job_force_mode == 12) then
                
                if (time_index > 1) then

                    dr_cen(:) = r0(:, cur) - f_cen(:)
                    dr_norm_sq = dr_cen(1)**2 + dr_cen(2)**2 + dr_cen(3)**2
                    dr_norm_inv = 1d0/sqrt(dr_norm_sq)
                    dr2_ex(:, cur) = (para_double(1)*(dr_cen(:)*dr_norm_inv))!(5d-5*time_exp)*exp(-dr_norm_sq/200d0)+
                    dr2_ex(2, cur) =  dr2_ex(2, cur) + sign(para_double(2), r0(2, cur))
                else                         
                    dr2_ex(:, cur) = 0d0*vect_one(:)
                end if   
                
            else if (my_job_force_mode == 13) then
                 if (time_index > 1) then

                    dr_cen(:) = r0(:, cur) - f_cen(:)
                    dr_norm_sq = dr_cen(1)**2 + dr_cen(2)**2 + dr_cen(3)**2
                    dr_norm_inv = 1d0/sqrt(dr_norm_sq)
                    dr2_ex(:, cur) = (para_double(1)*(dr_cen(:)*dr_norm_inv))!(5d-5*time_exp)*exp(-dr_norm_sq/200d0)+
                    if ((r0(2, cur) < para_double(3)) .or. (r0(2, cur) > para_double(4))) then
                        dr2_ex(2, cur) =  dr2_ex(2, cur) + sign(para_double(2), r0(2, cur))
                    end if
                else                         
                    dr2_ex(:, cur) = 0d0*vect_one(:)
                end if  
            else if (my_job_force_mode == 14) then
                 if (time_index > 1) then

                    dr_cen(:) = r0(:, cur) - f_cen(:)
                    dr_norm_sq = dr_cen(1)**2 + dr_cen(2)**2 + dr_cen(3)**2
                    dr_norm_inv = 1d0/sqrt(dr_norm_sq)
                    dr2_ex(:, cur) = (para_double(1)*(dr_cen(:)*dr_norm_inv))!(5d-5*time_exp)*exp(-dr_norm_sq/200d0)+
                    !if ((r0(2, cur) < para_double(3)) .or. (r0(2, cur) > para_double(4))) then
                    dr2_ex(2, cur) =  dr2_ex(2, cur) + para_double(2)*r0(2, cur)/srw(2,2)
                    !end if
                else                         
                    dr2_ex(:, cur) = 0d0*vect_one(:)
                end if       
            else if (my_job_force_mode == 15) then
                 if (time_index > 1) then

                    dr_cen(:) = r0(:, cur) - f_cen(:)
                    dr_norm_sq = dr_cen(1)**2 + dr_cen(2)**2 + dr_cen(3)**2
                    dr_norm_inv = 1d0/sqrt(dr_norm_sq)
                    dr2_ex(:, cur) = (para_double(1)*(dr_cen(:)*dr_norm_inv))!(5d-5*time_exp)*exp(-dr_norm_sq/200d0)+
                    !if ((r0(2, cur) < para_double(3)) .or. (r0(2, cur) > para_double(4))) then
                    dr2_ex(2, cur) =  dr2_ex(2, cur) + para_double(2)*r0(2, cur)/srw(2,2)
                    r0(1, cur) = r0(2,cur)* para_double(3) + r0(1, cur)
                    !end if
                else                         
                    dr2_ex(:, cur) = 0d0*vect_one(:)
                end if       
            else if (my_job_force_mode == 16) then
                !Shear
                 if (time_index > 1000) then
                    if (r0(2, cur) > para_double(1)) then
                        sign_x = r0(1, cur) - para_double(2)
                        mag = -para_double(4)*(sign_x) + sign(para_double(5),sign_x)
                        dr2_ex(:, cur) = mag * para_vect(:)
                        
                    else 
                        dr2_ex (:, cur) = 0d0*vect_one(:)
                    end if
                else                         
                    dr2_ex(:, cur) = 0d0*vect_one(:)
                end if                  
            else if (my_job_force_mode == 17) then
                if (time_index > 1000) then
                    if (s(cur) == 10) then
                        r0(1, cur) = (-para_double(1)) + r0(1, cur)
                    else if (s(cur) == 11) then
                        r0(1, cur) = (para_double(1)) + r0(1, cur)
                    end if
                end if                                          
             else if (my_job_force_mode == 18) then
                ! 0 2 3e-4 3e-4 3e-4
                if (time_index > 1000) then
                    if ((s(cur) == 101) .or. (s(cur) == 10)) then
                        r0(1, cur) = (-para_double(1)) + r0(1, cur)
                        r0(2, cur) = (-para_double(2)) + r0(2, cur)
                    else if ((s(cur) == 102) .or. (s(cur) == 11)) then
                        r0(1, cur) = (para_double(1)) + r0(1, cur)
                        r0(2, cur) = (para_double(2)) + r0(2, cur)
                    end if
                end if     
             else if (my_job_force_mode == 19) then
                ! 0 2 3e-4 3e-4 3e-4
                if (time_index > 1000) then
                    if ((s(cur) == 101) .or. (s(cur) == 10)) then
                        r0(1, cur) = (para_double(1)) + r0(1, cur)
                        r0(2, cur) = (para_double(2)) + r0(2, cur)
                        dr2_ex(:, cur) = mag * para_vect(:)
                    else if ((s(cur) == 102) .or. (s(cur) == 11)) then
                        r0(1, cur) = (para_double(3)) + r0(1, cur)
                        r0(2, cur) = (para_double(4)) + r0(2, cur)
                    else if ((s(cur) == 103) .or. (s(cur) == 12)) then
                        r0(1, cur) = (para_double(5)) + r0(1, cur)
                        r0(2, cur) = (para_double(6)) + r0(2, cur)    
                        
                    end if
                end if
            else if (my_job_force_mode == 20) then
                ! 0 2 3e-4 3e-4 3e-4
                if ((mod(time_index,10) == 0) .and. (time_index .ne. 0)) then
                    if ((s(cur) == 101) .or. (s(cur) == 10)) then
                        r0(1, cur) = (-para_double(1)) + r0(1, cur)
                        r0(2, cur) = (-para_double(2)) + r0(2, cur)
                    else if ((s(cur) == 102) .or. (s(cur) == 11)) then
                        r0(1, cur) =  (para_double(1)) + r0(1, cur)
                        r0(2, cur) =  (para_double(2)) + r0(2, cur)
                    end if      
                end if
            else if (my_job_force_mode == 21) then
                ! 0 2 1.001 1.0
                if ((mod(time_index,50) == 0) .and. (time_index .ne. 0)) then
                    
                    r0(1, cur) = para_double(1)*r0(1, cur)
                    r0(2, cur) = para_double(2)*r0(2, cur)
                    
                end if                                   
            
            end if
            
                
    end subroutine
    
    !Force Mode Initiation
    subroutine MPI_MD_force_mode_init()

        f_cen(:) = (sr(2,:) + sr(1,:))/2d0

        if (myid .ne. 0) then
            dr2_ex(:, :) = 0d0*dr2_ex(:, :)
            
            if (my_job_force_mode == 17) then
                para_vect(1) = cos(para_double(3)*PI/180d0)
                para_vect(2) = sin(para_double(3)*PI/180d0)
                para_vect(3) = 0d0
            end if
        end if

    
    end subroutine
    
    !Temperature Control
    subroutine MD_temperature_control()
        double precision start_time1, stop_time1, elasped_time1
        !Add special T_control instruction here
        !start_time1 = MPI_Wtime()
        
        !if (time_index >= 500) then
        if (mod(time_index + 250, 500) == 0) then
            T_ratio_sq = sqrt(T_wanted/T_current)
        else
            T_ratio_sq = 1
        End if
        !sqrt((T_wanted - T_current)*dt/(T_current)/1d0 + 1d0)
        !end if
        
        !stop_time1 = MPI_Wtime()
        !elasped_time1 = stop_time1 - start_time1
        !time_temperature_control = time_temperature_control + elasped_time1
        !operation_time = operation_time + elasped_time1
    end subroutine
        
    !Fold Z shadow    
    subroutine MD_atom_folding_XP()
        integer i, cur
        double precision r0_temp(3)
        double precision cache_zero(3) /0.0, 0.0, 0.0/
        
        
        
        cur = ch_start

        do i = 1, ch_len
            
            ! Fold Z and Y atoms
            if (s(cur) > 0) then

                if (r0(3, cur) < sr(1, 3)) then
                    r0(3,cur) = r0(3,cur) + srw(3,3)
                else if (r0(3, cur) > sr(2, 3)) then
                    r0(3,cur) = r0(3,cur) - srw(3,3)
                end if

                if (r0(2, cur) < sr(1, 2)) then
                    r0(2,cur) = r0(2,cur) + srw(2,2)
                else if (r0(2, cur) > sr(2, 2)) then
                    r0(2,cur) = r0(2,cur) - srw(2,2)
                end if
                
            end if 
        
            ! Fold Front Z face Shadow
            if ((r0(3, cur) < (sr(1, 3) + f_cut_2)) .and. (r0(3, cur) >= sr(1,3))) then
                if (r0(3, cur) <= (sr(1, 3) + f_cut)) then
                    call add_atom(r0(:,cur) + srw(:,3),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, t(cur), 0, 0)  
                else
                    call add_atom(r0(:,cur) + srw(:,3),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, t(cur), -1, 0) 
                end if  
  
            end if
            
            ! Fold Back Z face Shadow
            if ((r0(3, cur) > (sr(2, 3) - f_cut_2)) .and. (r0(3, cur) <= sr(2,3))) then
                if (r0(3, cur) >= (sr(2, 3) - f_cut)) then
                    call add_atom(r0(:,cur) - srw(:,3),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, t(cur), 0, 0)  
                else
                    call add_atom(r0(:,cur) - srw(:,3),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, t(cur), -1, 0) 
                end if    
  
            end if
#ifdef Perodic_Y            
            ! Fold Front Z face Shadow
            if ((r0(2, cur) < (sr(1, 2) + f_cut_2)) .and. (r0(2, cur) >= sr(1,2))) then
                if (r0(2, cur) <= (sr(1, 2) + f_cut)) then
                    call add_atom(r0(:,cur) + srw(:,2),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, t(cur), 0, 0)  
                else
                    call add_atom(r0(:,cur) + srw(:,2),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, t(cur), -1, 0) 
                end if  
  
            end if
            
            ! Fold Back Z face Shadow
            if ((r0(2, cur) > (sr(2, 2) - f_cut_2)) .and. (r0(2, cur) <= sr(2,2))) then
                if (r0(2, cur) >= (sr(2, 2) - f_cut)) then
                    call add_atom(r0(:,cur) - srw(:,2),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, t(cur), 0, 0)  
                else
                    call add_atom(r0(:,cur) - srw(:,2),cache_zero,cache_zero,&
                                 &cache_zero, cache_zero, cache_zero, t(cur), -1, 0) 
                end if    
  
            end if
#endif

        cur = next(cur)
        end do    
    end subroutine    

    !Gear Prediction
    subroutine MD_predict()
        integer i, cur
        double precision start_time1, stop_time1, elasped_time1
        
        start_time1 = MPI_wtime()
        cur = ch_start
        
        time_exp = (2d0**int(time_index/2000))
        
        



        do i = 1, ch_len
        call apply_strain_stress(cur)

        if ((s(cur) >= 1) .and. (s(cur) <= 20)) then
        r1(:,cur) = T_ratio_sq * r1(:, cur)
        
        r0(:,cur) = r0(:,cur) + dt*r1(:,cur) + dt2d2*r2(:,cur) + dt3d6*r3(:,cur) + dt4d24*r4(:,cur) + dt5d120*r5(:,cur)
        r1(:,cur) = r1(:,cur) + dt*r2(:,cur) + dt2d2*r3(:,cur) + dt3d6*r4(:,cur) + dt4d24*r5(:,cur)                 
        r2(:,cur) = r2(:,cur) + dt*r3(:,cur) + dt2d2*r4(:,cur) + dt3d6*r5(:,cur)
        r3(:,cur) = r3(:,cur) + dt*r4(:,cur) + dt2d2*r5(:,cur)
        r4(:,cur) = r4(:,cur) + dt*r5(:,cur)
        
        end if
        cur = next(cur)
        end do
        
        stop_time1 = MPI_Wtime()
        elasped_time1 = stop_time1 - start_time1
        time_predict = time_predict + elasped_time1
        operation_time = operation_time + elasped_time1
    
    end subroutine

    !Calculate Force    
    subroutine MD_EAM_force_XP()
        integer xp, xn, yp, yn, zp,zn, i, j, k, n1, n2
        integer cur1, cur2
        integer ix(3)
        double precision dx(3), dx_inv(3), dr0, dr0_inv, dr0_sq
        double precision dFj_cache(3, 400), dF_cache
        double precision dFi_sum(3), dVi_sum(3), drho(3), Vi_sum
        double precision rho_sum, KE_sum
        double precision start_time1, stop_time1, elasped_time1
        integer n_count, n_list(400)       

     
        
        start_time1 = MPI_wtime()
        
        !Outter Loop
        Force(:, :) = 0d0*Force(:, :)
        cur1 = ch_start
        do n1 =1, ch_len

            dFi_sum(:) = 0d0*vect_one(:)
            dVi_sum(:) = dFi_sum(:)
            rho_sum = 0d0
            Vi_sum = 0d0
            n_count = 0

      
            
            !cur1 = s1_list(n1)
            if ((s(cur1) >= 0) .and. (s(cur1) <= 20)) then
                
                ix(:) = to_int3( (r0(:, cur1) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))
                    !Set the upper bound and lower bound of the index
                if (ix(1) > 0) then
                    xn = ix(1)-1
                else 
                    xn = 0
                end if
                    
                if (ix(1) < NOD(1)) then
                    xp = ix(1)+1
                else
                    xp = NOD(1)
                end if
                    
                if (ix(2) > 0) then
                    yn = ix(2)-1
                else 
                    yn = 0
                end if
                    
                if (ix(2) < NOD(2)) then
                    yp = ix(2)+1
                else
                    yp = NOD(2)
                end if

                if (ix(3) > 0) then
                    zn = ix(3)-1
                else 
                    zn = 0
                end if
                    
                if (ix(3) < NOD(3)) then
                    zp = ix(3)+1
                else
                    zp = NOD(3)
                end if
                
                KE_sum = KE_sum + m_d3k(t(cur1))*(r1(1, cur1)**2 + r1(2, cur1)**2 + r1(3, cur1)**2)

                do k = zn, zp                    
                do j = yn, yp
                do i = xn, xp

                if (nbox_time(i,j,k) == time_index) then
                            
                        
                    cur2 = nbox_start(i,j,k)
                                                    
                    !Inner Loop
                    do n2 = 1, nbox_count(i,j,k)
                    
                    if (cur1 .ne. cur2) then 
                            
                        dx(:) = r0(:,cur1)-r0(:,cur2)

                        dr0_sq = (dx(1)**2 + dx(2)**2 + dx(3)**2)
                                    !if (myid == 3) write (*, fmt1) "dr 1",dr0_sq, r0(1,cur1)- br(myid, 1, 1), r0(1,cur2)- br(myid, 1, 1),&
                                    !                                       & cur1,cur2,s(cur2),i,j,j,&
                                    !                                       & nbox_start(i,j,k), nbox_next(cur2)
                        if ((dr0_sq < f_cut_sq_max(t(cur1)+t(cur2))) .and. (cur1 .ne. cur2))then
                                    
                            dr0 = sqrt(dr0_sq)
                            dr0_inv = 1d0/dr0
                            dx_inv(:) = dx(:)*dr0_inv
                            rho_sum = rho_sum + get_rho(t(cur2), dr0)
                            
                            
                            dFi_sum(:) = dFi_sum(:) + get_drho(t(cur2), dr0)*dx_inv(:)
                            dVi_sum(:) = dVi_sum(:) + get_dv(t(cur1), t(cur2), dr0)*dx_inv(:)
                            Vi_sum = Vi_sum + get_v(t(cur1), t(cur2), dr0)
                            
                                  
                            n_count = n_count + 1
                            n_list(n_count) = cur2
                            dFj_cache(:, n_count) = get_drho(t(cur1), dr0)*dx_inv(:)

                                        
                        end if
                    end if
                    
                    cur2 = nbox_next(cur2)
                    end do
                        
                end if
                                 
                end do
                end do
                end do
                
                dF_cache = get_dF(t(cur1), rho_sum)
                dFi_sum(:) = (-1d0)*dF_cache*dFi_sum(:)
                dVi_sum(:) = (-0.5d0)*dVi_sum(:)
                Force(:, cur1) = Force(:, cur1) +  dFi_sum(:) + dVi_sum(:)
                Energy(cur1) = (0.5d0)*Vi_sum + get_F(t(cur1), rho_sum)

               
                do j = 1, n_count
                    Force(:, n_list(j)) = Force(:, n_list(j)) + dF_cache*dFj_cache(:,j)
                end do
  
           
            end if   
            
        cur1 = next(cur1)
 
        end do
 
        if (NOAA .ne. 0) then
            T_current = KE_sum / dble(NOAA)
        end if
 
        stop_time1 = MPI_Wtime()
        elasped_time1 = stop_time1 - start_time1
        time_EAM = time_EAM + elasped_time1
        operation_time = operation_time + elasped_time1
    end subroutine
    
    ! Force with Stress 
    subroutine MD_EAM_force_stress_XP()
        integer xp, xn, yp, yn, zp,zn, i, j, k, n1, n2, m, n
        integer cur1, cur2
        integer ix(3)
        integer n_count, n_list(400)            
        double precision dx(3), dx_inv(3), dr0, dr0_inv, dr0_sq
        double precision dFj_cache(3, 400), dF_cache, dx_cache(3, 400)
        double precision dFi_sum(3), dVi_sum(3), drho(3), Vi_sum
        double precision rho_sum, KE_sum
        double precision stress_tmp(3,3,3), dv_tmp(3), drho_tmp(3)
        double precision start_time1, stop_time1, elasped_time1
        
        
        
        start_time1 = MPI_wtime()
        
        !Outter Loop
        Force(:, :) = 0d0*Force(:, :)
        Stress(:, :, :) = 0d0 * Stress(:, :, :)
                
        cur1 = ch_start
        do n1 =1, ch_len

            dFi_sum(:) = 0d0*vect_one(:)
            dVi_sum(:) = dFi_sum(:)
            rho_sum = 0d0
            Vi_sum = 0d0
            n_count = 0

      
            
            !cur1 = s1_list(n1)
            !if ((s(cur1) >= 0) .and. (s(cur1) <= 20)) then
                
            stress_tmp(:, :, :) = 0d0*stress_tmp(:, :, :)
            
            ix(:) = to_int3( (r0(:, cur1) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))
                !Set the upper bound and lower bound of the index
            if (ix(1) > 0) then
                xn = ix(1)-1
            else 
                xn = 0
            end if
                
            if (ix(1) < NOD(1)) then
                xp = ix(1)+1
            else
                xp = NOD(1)
            end if
                
            if (ix(2) > 0) then
                yn = ix(2)-1
            else 
                yn = 0
            end if
                
            if (ix(2) < NOD(2)) then
                yp = ix(2)+1
            else
                yp = NOD(2)
            end if

            if (ix(3) > 0) then
                zn = ix(3)-1
            else 
                zn = 0
            end if
                
            if (ix(3) < NOD(3)) then
                zp = ix(3)+1
            else
                zp = NOD(3)
            end if
            
            KE_sum = KE_sum + m_d3k(t(cur1))*(r1(1, cur1)**2 + r1(2, cur1)**2 + r1(3, cur1)**2)

            do k = zn, zp                
            do j = yn, yp
            do i = xn, xp

            if (nbox_time(i,j,k) == time_index) then
                        
                    
                cur2 = nbox_start(i,j,k)
                                                
                !Inner Loop
                do n2 = 1, nbox_count(i,j,k)
                
                if (cur1 .ne. cur2) then 
                        
                    dx(:) = r0(:,cur1)-r0(:,cur2)

                    dr0_sq = (dx(1)**2 + dx(2)**2 + dx(3)**2)
                                !if (myid == 3) write (*, fmt1) "dr 1",dr0_sq, r0(1,cur1)- br(myid, 1, 1), r0(1,cur2)- br(myid, 1, 1),&
                                !                                       & cur1,cur2,s(cur2),i,j,j,&
                                !                                       & nbox_start(i,j,k), nbox_next(cur2)
                    if ((dr0_sq < f_cut_sq_max(t(cur1)+t(cur2))) .and. (cur1 .ne. cur2))then
                                
                        dr0 = sqrt(dr0_sq)
                        dr0_inv = 1d0/dr0
                        dx_inv(:) = dx(:)*dr0_inv
                        rho_sum = rho_sum + get_rho(t(cur2), dr0)
                        
                        dv_tmp(:) = get_dv(t(cur1), t(cur2), dr0)*dx_inv(:)
                        !drho_tmp(:) = get_drho(t(cur2), dr0)*dx_inv(:)
                        
                        dVi_sum(:) = dVi_sum(:) + dv_tmp(:)
                        dFi_sum(:) = dFi_sum(:) + get_drho(t(cur2), dr0)*dx_inv(:)
                        Vi_sum = Vi_sum + get_v(t(cur1), t(cur2), dr0)
                        
                              
                        n_count = n_count + 1
                        n_list(n_count) = cur2
                        dx_cache(:, n_count) = dx(:)
                        
                        dFj_cache(:, n_count) = get_drho(t(cur1), dr0)*dx_inv(:)

                        do m = 1, 3
                          stress_tmp(m, :, 1) = stress_tmp(m, :, 1) + dv_tmp(m)*dx(:)
                          stress_tmp(m, :, 2) = stress_tmp(m, :, 2) + get_drho(t(cur1), dr0)*dx_inv(m)*dx(:)
                        end do
                                    
                    end if
                end if
                
                cur2 = nbox_next(cur2)
                end do
                    
            end if
                             
            end do
            end do
            end do
            
            dF_cache = get_dF(t(cur1), rho_sum)
            dFi_sum(:) = (-1d0)*dF_cache*dFi_sum(:)
            dVi_sum(:) = (-0.5d0)*dVi_sum(:)
            Force(:, cur1) = Force(:, cur1) +  dFi_sum(:) + dVi_sum(:)
            Energy(cur1) = (0.5d0)*Vi_sum + get_F(t(cur1), rho_sum)
            Stress(:, :, cur1) = Stress(:, :, cur1) + (-0.5d0)*stress_tmp(:,:,1) + (-1d0)*dF_cache*stress_tmp(:,:,2)
           
            do j = 1, n_count
                Force(:, n_list(j)) = Force(:, n_list(j)) + dF_cache * dFj_cache(:,j)
                do m = 1, 3    
                    Stress(m, :, n_list(j)) = Stress(m, :, n_list(j)) + (-1d0)*dF_cache * dFj_cache(m,j)*dx_cache(:, j)
                end do
                
            end do

           
            !end if   
            
        cur1 = next(cur1)
 
        end do
 
        if (NOAA .ne. 0) then
            T_current = KE_sum / dble(NOAA)
        end if
 
        stop_time1 = MPI_Wtime()
        elasped_time1 = stop_time1 - start_time1
        time_EAM = time_EAM + elasped_time1
        operation_time = operation_time + elasped_time1
    end subroutine
 
    ! Average out all Stresses
    subroutine MD_Stress_Average()
        integer xp, xn, yp, yn, zp,zn, i, j, k, n1, n2, m, n
        integer cur1, cur2
        integer ix(3), ncount
        double precision dx(3), dx_inv(3), dr0, dr0_inv, dr0_sq
        double precision sa(3,3)
        !double precision Fi 
        !double precision Vi 
        double precision dFi_sum(3), dVi_sum(3), Vi_sum
        double precision rho_sum, KE_sum
        double precision start_time1, stop_time1, elasped_time1
        double precision Iv1_tmp, Iv2_tmp, Iv3_tmp
        double precision stress_tmp(3,3)
        
        start_time1 = MPI_wtime()
        cur1 = ch_start
        
        !Stress_ave(:, :, :) =  Stress(:, :, :)
        do n1 =1, ch_len    

            ncount = 1
            
            stress_tmp(:, :) = stress(:, :, cur1) 

            !if ((s(cur1) == 1) .or. ((s(cur1) >= 2) .and. free_pin)) then
                                 
                
                ix(:) = to_int3( (r0(:, cur1) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))
                    !Set the upper bound and lower bound of the index
                if (ix(1) > 0) then
                    xn = ix(1)-1
                else 
                    xn = 0
                end if
                    
                if (ix(1) < NOD(1)) then
                    xp = ix(1)+1
                else
                    xp = NOD(1)
                end if
                    
                if (ix(2) > 0) then
                    yn = ix(2)-1
                else 
                    yn = 0
                end if
                    
                if (ix(2) < NOD(2)) then
                    yp = ix(2)+1
                else
                    yp = NOD(2)
                end if

                if (ix(3) > 0) then
                    zn = ix(3)-1
                else 
                    zn = 0
                end if
                    
                if (ix(3) < NOD(3)) then
                    zp = ix(3)+1
                else
                    zp = NOD(3)
                end if
                
                !KE_sum = KE_sum + m_d3k(t(cur1))*(r1(1, cur1)**2 + r1(2, cur1)**2 + r1(3, cur1)**2)
                    
                    !if (myid == 4) print *, "p2", time_index, ix(1), ix(2), ix(3)
                do i = xn, xp
                do j = yn, yp
                do k = zn, zp
                        !if (myid == 4) print *, "p2b", time_index, ix(1), ix(2), ix(3),i, j,k, cur2
                if (nbox_time(i,j,k) == time_index) then
                            
                        
                    cur2 = nbox_start(i,j,k)
                                                    
                            !if (myid == 3) print *, "p3", time_index, ix(1), ix(2), ix(3),i, j,k, cur2

                    do n2 = 1, nbox_count(i,j,k)
                        if (cur1 .ne. cur2) then
                        !if ((s(cur2) == 1) .or. ((s(cur2) == 2) .and. free_pin)) then
                            dx(:) = r0(:,cur1)-r0(:,cur2)

                            dr0_sq = (dx(1)**2 + dx(2)**2 + dx(3)**2)
                                        !if (myid == 3) write (*, fmt1) "dr 1",dr0_sq, r0(1,cur1)- br(myid, 1, 1), r0(1,cur2)- br(myid, 1, 1),&
                                        !                                       & cur1,cur2,s(cur2),i,j,j,&
                                        !                                       & nbox_start(i,j,k), nbox_next(cur2)
                            if ((dr0_sq < f_cut_sq) .and. (cur1 .ne. cur2))then
                                        
                                stress_tmp(:,:) = stress_tmp(:,:) + Stress(:,:,cur2)
                                ncount = ncount + 1
                                            
                            end if
                            
                        end if
                        cur2 = nbox_next(cur2)
                    end do

                        
                end if
                                 
                end do
                end do
                end do
                
                Stress_ave(:,:,cur1) = Stress_tmp(:,:)/dble(ncount)/dble(ncount)

            !end if   
   

            cur1 = next(cur1)
 
        end do
 
        !print *, "A6b", myid, NOAA, NOPA, "A6b"
        cur1 = ch_start

        do n1 =1, ch_len
            sa(:,:) = Stress_ave(:, :, cur1)
            Iv1(cur1) = sa(1,1) + sa(2,2) + sa(3,3)
            Iv2(cur1) = sa(1,1)*sa(2,2) + sa(2,2)*sa(3,3) + sa(3,3)*sa(1,1)
            Iv2(cur1) = Iv2(cur1) - sa(1,2)**2 - sa(2,3)**2 - sa(3,1)**2
            Iv3(cur1) = sa(1,1)*sa(2,2)*sa(3,3) + 2d0*sa(1,2)*sa(2,3)*sa(3,1)
            Iv3(cur1) = Iv3(cur1) - (sa(1,2)**2)*sa(3,3)-(sa(2,3)**2)*sa(1,1)-(sa(1,3)**2)*sa(2,2)
            Jv2(cur1) = (Iv1(cur1)**2)/(3d0) - Iv2(cur1)
            Jv3(cur1) = (2d0/27d0)*Iv1(cur1)**3 - Iv1(cur1)*Iv2(cur1)/3d0 + Iv3(cur1)
            cur1 = next(cur1)
        end do
 
        stop_time1 = MPI_Wtime()
        elasped_time1 = stop_time1 - start_time1
        time_EAM = time_EAM + elasped_time1
        operation_time = operation_time + elasped_time1
    end subroutine
    
    !Gear Correct
    subroutine MD_correct()
        integer i, cur
        double precision start_time1, stop_time1, elasped_time1
        
        
#ifdef HARD_PARTICLE_MODE
        integer HPM_i
        double precision dr0_rot(3)
        
        do i = 1, HPM_NOP
            HPM_dr0(:, i) = dt2d2*(HPM_force(:, i)/HPM_mass(i))
        end do
#ifdef HARD_PARTICLE_torque
        HPM_dtheta(:) = dt2d2*(HPM_torque(:)/HPM_mass(:))
#endif
#endif

        
        start_time1 = MPI_wtime()
        cur = ch_start
        !do i = 1, s1_count
        do i = 1, ch_len
        !cur = s1_list(i)
            if ((s(cur) >= 1) .and. (s(cur) <= 20)) then
            dr2(:, cur) = q_dm(t(cur))*Force(:, cur) - r2(:, cur)
            r0(:,cur) = r0(:,cur) + a0*(dr2(:,cur) + dr2_ex(:,cur))
            r1(:,cur) = r1(:,cur) + a1*(dr2(:,cur) + dr2_ex(:,cur))
            r2(:,cur) = r2(:,cur) +    (dr2(:,cur) + dr2_ex(:,cur))
            r3(:,cur) = r3(:,cur) + a3*(dr2(:,cur) + dr2_ex(:,cur))
            r4(:,cur) = r4(:,cur) + a4*(dr2(:,cur) + dr2_ex(:,cur))
            r5(:,cur) = r5(:,cur) + a5*(dr2(:,cur) + dr2_ex(:,cur))
            r1(:,cur) = r1(:,cur) - r1_ave(:)
            
#ifdef HARD_PARTICLE_MODE
            else if ((s(cur) > 200) .and. (s(cur) < 210)) then

            HPM_i = s(cur)-200
            
#ifdef HARD_PARTICLE_torque
            dr0_rot(:) = HPM_dtheta *(r0(:, cur) - HPM_cm(:, HPM_i))
            r0(1, cur) = r0(1, cur) - HPM_dtheta(HPM_i)*dr0_rot(2)
            r0(2, cur) = r0(2, cur) + HPM_dtheta(HPM_i)*dr0_rot(1)
#endif                    
            r0(:, cur) = r0(:, cur) + HPM_dr0(:, HPM_i)
#endif 
            
            
            
            end if        
        cur = next(cur)
        end do
        r1_ave(:) = 0d0*vect_one(:)
        
        stop_time1 = MPI_Wtime()
        elasped_time1 = stop_time1 - start_time1
        time_correct = time_correct + elasped_time1
        operation_time = operation_time + elasped_time1
        
    end subroutine
    
    !Remove Shadow()
    subroutine MD_remove_shadow()
        double precision start_time, stop_time, elasped_time
        integer i, n, my_next
        start_time = MPI_wtime()
        n = ch_start
        do i = 1, ch_len
            my_next = next(n)
            if (s(n) <= 0)then
                !if (.not. in_domain(r0(:,n)))then
                call remove_atom(n)
                !end if
            end if
            n = my_next
        end do
        
        stop_time = MPI_wtime()
        elasped_time = stop_time - start_time
        time_remove_shadow = time_remove_shadow + elasped_time
        operation_time = operation_time + elasped_time
        
    end subroutine
    
    !Remove Escaped Atoms
    subroutine MD_remove_atom()
        double precision start_time, stop_time, elasped_time
        integer i, si
        
        start_time = MPI_Wtime()
        
        do i = 1, satom_count
        si = satom_list(i) + 1
        call remove_atom(si)
        end do
    
        stop_time = MPI_wtime()
        elasped_time = stop_time - start_time
        time_remove_atom = time_remove_atom + elasped_time
        operation_time = operation_time + elasped_time
        
    end subroutine
    
    !Shadow Boundary folding
    subroutine MD_shadow_boundary_folding(rt1, n, rt2)
        double precision rt1(3), rt2(3,8)
        integer n, i
        logical in(6)
        n = 1
        rt2(:,1) = rt1(:)
        
        in(1) = (rt1(1) < sr(1,1) + f_cut)
        in(2) = (rt1(1) >= sr(2,1) - f_cut)
        in(3) = (rt1(2) < sr(1,2) + f_cut)
        in(4) = (rt1(2) >= sr(2,2) - f_cut)
        in(5) = (rt1(3) < sr(1,3) + f_cut)
        in(6) = (rt1(3) >= sr(2,3) - f_cut)
        
        if(in(1)) then
            select case (b_type(1,1))
            case (0)
            
            case (1)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i) + srw(1,:)
                end do
                n = 2*n
            case (2)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i)
                    rt2(1, n+i) = sr2(1,1) - rt2(1, n+i) 
                end do
                n = 2*n
            case default
                print *, "ERROR, UNRECONIZED BOUNDARY"
            end select
        else if (in(2)) then
            select case (b_type(1,2))
            case (0)
            
            case (1)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i) - srw(1,:)
                end do
                n = 2*n
            case (2)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i)
                    rt2(1, n+i) = sr2(2,1) - rt2(1, n+i) 
                end do
                n = 2*n
            case default
                print *, "ERROR, UNRECONIZED BOUNDARY"
            end select
        end if

        if(in(3)) then
            select case (b_type(2,1))
            case (0)
            
            case (1)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i) + srw(2,:)  
                end do
                n = 2*n
            case (2)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i)
                    rt2(2, n+i) = sr2(1,2) - rt2(2, n+i) 
                end do
                n = 2*n
            case default
                print *, "ERROR, UNRECONIZED BOUNDARY"
            end select
        else if (in(4)) then
            select case (b_type(2,2))
            case (0)
            
            case (1)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i) - srw(2,:)
                end do
                n = 2*n
            case (2)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i)
                    rt2(2, n+i) = sr2(2,2) - rt2(2, n+i) 
                end do
                n = 2*n
            case default
                print *, "ERROR, UNRECONIZED BOUNDARY"
            end select
        end if
        
       if(in(5)) then
            select case (b_type(3,1))
            case (0)
            
            case (1)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i) + srw(3,:)
                end do
                n = 2*n
            case (2)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i)
                    rt2(3, n+i) = sr2(1,3) - rt2(3, n+i) 
                end do
                n = 2*n
            case default
                print *, "ERROR, UNRECONIZED BOUNDARY"
            end select
        else if (in(6)) then
            select case (b_type(3,2))
            case (0)
            
            case (1)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i) - srw(3,:)
                end do
                n = 2*n
            case (2)
                do i = 1, n
                    rt2(:, n+i) = rt2(:, i)
                    rt2(3, n+i) = sr2(2,3) - rt2(3, n+i) 
                end do
                n = 2*n
            case default
                print *, "ERROR, UNRECONIZED BOUNDARY"
            end select
        end if
        
    end subroutine
    
    !Boundary folding
    subroutine MD_boundary_folding()
        double precision start_time, stop_time, elasped_time
        integer i
                
        ! Perform Boundary folding to all active atoms (s = 1)
        ! Loop unrolled to improve efficiency
        start_time = MPI_Wtime()
        
        do i = 1, NOS
            if (s(i) > 0) then
               
                !X Axis boundary
                if (r0(1, i) < sr(1,1)) then
                    select case(b_type(1,1))
                    case (0) !No Boundary
                    
                    case (1) !Periodic Boundary
                        r0(1, i) = r0(1, i) + sr(2,1) - sr(1,1)
                    case (2) !Mirror Boundary
                        r0(1, i) = sr2(1,1) - r0(1, i)                    
                    case default
                        print *, "ERROR, UNRECONIZED BOUNDARY"
                    end select
                else if (sr(2,1) <= r0(1, i))then
                    select case(b_type(1,2))
                    case (0) !No Boundary
                    
                    case (1) !Periodic Boundary
                        r0(1, i) = r0(1, i) - sr(2,1) + sr(1,1)
                    case (2) !Mirror Boundary
                        r0(1, i) = sr2(2,1) - r0(1, i)                    
                    case default
                        print *, "ERROR, UNRECONIZED BOUNDARY"
                    end select
                end if
                
                !Y Axis boundary
                if (r0(2, i) < sr(1, 2)) then
                    select case(b_type(2,1))
                    case (0) !No Boundary
                    
                    case (1) !Periodic Boundary
                        r0(2, i) = r0(2, i) + sr(2,2) - sr(1,2)
                    case (2) !Mirror Boundary
                        r0(2, i) = sr2(1,2) - r0(2, i)                    
                    case default
                        print *, "ERROR, UNRECONIZED BOUNDARY"
                    end select                    
                else if (sr(2,2) <= r0(2, i)) then
                    select case(b_type(2,2))
                    case (0) !No Boundary
                    
                    case (1) !Periodic Boundary
                        r0(2, i) = r0(2, i) - sr(2,2) + sr(1,2)
                    case (2) !Mirror Boundary
                        r0(2, i) = sr2(2,2) - r0(2, i)                    
                    case default
                        print *, "ERROR, UNRECONIZED BOUNDARY"
                    end select
                end if
                
                !Z Axis boundary
                if (r0(3, i) < sr(1, 3)) then
                    select case(b_type(3,1))
                    case (0) !No Boundary
                    
                    case (1) !Periodic Boundary
                        r0(3, i) = r0(3, i) + sr(2,3) - sr(1,3)
                    case (2) !Mirror Boundary
                        r0(3, i) = sr2(1,3) - r0(3, i)                    
                    case default
                        print *, "ERROR, UNRECONIZED BOUNDARY"
                    end select
                else if (sr(2,3) <= r0(3, i)) then
                    select case(b_type(3,2))
                    case (0) !No Boundary
                    
                    case (1) !Periodic Boundary
                        r0(3, i) = r0(3, i) - sr(2,3) + sr(1,3)
                    case (2) !Mirror Boundary
                        r0(3, i) = sr2(2,3) - r0(3, i)                    
                    case default
                        print *, "ERROR, UNRECONIZED BOUNDARY"
                    end select
                end if
               
            
            end if
        end do
        stop_time = MPI_wtime()
        elasped_time = stop_time - start_time
        time_boundary_folding = time_boundary_folding + elasped_time
        operation_time = operation_time + elasped_time     
        
    end subroutine
         
    !Generate Shadow List to be send for Linear Mode
    subroutine MD_generate_send_list_Linear_XP(sign1)
        integer i, n, sign1
        double precision edge1, edge2, edge3, edge4
        
        edge1 = br(myid, 1, 1)
        edge2 = br(myid, 1, 1) + f_cut_2
        edge3 = br(myid, 2, 1) - f_cut_2
        edge4 = br(myid, 2, 1)
        
        sshadow_count = 0
        satom_count = 0
      
        
        n = ch_start
        if (sign1  == 1) then
            do i = 1, ch_len
                !debug            

                if (s(n) > 0) then
                
                    if ((r0(1,n) >= edge3) .and. (r0(1, n) < edge4)) then
                        sshadow_count = sshadow_count + 1
                        sshadow_list(sshadow_count) = n-1
                    end if
                    
                    if (r0(1,n) >= edge4) then
                        satom_count = satom_count + 1
                        satom_list(satom_count) = n-1
                        image_count = image_count + 1
                        image_list(image_count) = n
                    end if
                                        
                end if
                
                n = next(n)
            end do
        else
            do i = 1, ch_len
                !debug            

                if (s(n) > 0) then
                    if ((r0(1,n) < edge2) .and. (r0(1, n) >= edge1)) then
                        sshadow_count = sshadow_count + 1
                        sshadow_list(sshadow_count) = n-1
                    end if
                    

                    if (r0(1,n) < edge1) then
                        satom_count = satom_count + 1
                        satom_list(satom_count) = n-1
                        image_count = image_count + 1
                        image_list(image_count) = n
                    end if
                    
                end if
                
                n = next(n)
            end do
        
        end if
        
    end subroutine

    !Special Functions
    subroutine MPI_MD_MEASURE_elastic()
        double precision M_F_Ave(3,2), M_double(3,2), M_F_dot(2), M_Energy
        integer M_counter(3), M_int(3), i1, myfile, cur
        character*500 mystring
        
        M_Energy = 0d0
        M_counter(:) = 0 * M_counter(:)
        M_int(:) = 0 * M_int(:)
        M_F_ave(:,:) = 0d0 * M_F_ave(:,:)
        M_double(:,:) = 0d0 * M_double(:,:)
        
        if (myid == 0) then
        
            myfile = 6000 + my_job_id
        
            if (time_index == 0) then
                open (UNIT = (myfile), FILE=root_path(1:root_len)//"out/measurements.txt", STATUS = "replace")
                print *, "File Created", root_path(1:root_len)//"out/measurements.txt"
            end if
            
            call MPI_Reduce(M_double, M_F_ave, 6, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(M_double, M_Energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(M_int, M_counter, 3, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            
            M_F_dot(:) = M_F_ave(1,:)*0.52223071303152020074280170549084 + M_F_ave(2,:)*0.85280424621749509507957435287515
            
            write (mystring, *) time_index, M_Energy, M_F_ave(:,1), M_F_dot(1),  M_F_ave(:, 2), M_F_dot(2), M_counter(:)
            write (myfile, '(A)') TRIM(ADJUSTL(mystring))
            if (time_index == end_time) then
                close (myfile)
            end if
        else 
        
            !Sum all Forces
            
            cur = ch_start
            do i1 = 1, ch_len
                if ((r0(2, cur)  < sr(2, 2) - 3*f_cut) .and. ((r0(2, cur)  > sr(1, 2) + 3*f_cut))) then
                    if (s(cur) == 101) then
                        M_F_ave(:,1) = M_F_ave(:,1) + Force(:, cur)
                        M_counter(1) = M_counter(1) + 1
                    else if (s(cur) == 102) then
                        M_F_ave(:,2) = M_F_ave(:,2) + Force(:, cur)
                        M_counter(2) = M_counter(2) + 1
                    else if (s(cur) == 1) then
                        M_counter(3) = M_counter(3) + 1
                        M_Energy = M_Energy + Energy(cur)
                    end if
                end if
                cur = next(cur)
            end do         
            
            call MPI_Reduce(M_F_ave, M_double, 6, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(M_Energy, M_double,1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(M_counter, M_int, 3, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
               
        end if
        
        
        
    end subroutine 

    !Atom Sorting
    subroutine MD_Atom_Chain_Sorting()
        use qsort_library
        integer i
        
        if (mod(time_index, 200) .eq. 0) then
        
            !call MPI_Dump_Atoms(1, NOS, time_index, myid*100)
            
            !call debug_find_my_atom()
            
            call Qsort_Atoms(1, NOS, 0)
        
            call Qsort_Rebuild_Chain()
            
            !call debug_find_my_atom()
            
            !call MPI_Dump_Atoms(1, NOS, time_index, myid*100 + 1)
        end if
    end subroutine
   
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MPI Implemented Debug Subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  
#ifdef COMPILE_DEBUG
 


    ! Dump parameters onto screen   
    subroutine MPI_MD_debug_print_parameter()
        print *, myid, f_cut, dt, end_time, movie_time
    end subroutine
    
   ! Dump all variables into file
    subroutine MPI_MD_debug_dump_data()
        integer myfile, i, cur
        character*20 str_node
        
        
        myfile = 500 + myid
        
        write (str_node, "(I2)") myid
        str_node = "debug_atoms"//str_node(1:2)//".txt"
        open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/"//str_node(1:18), STATUS="replace")
        write (myfile, *) ch_start, ch_end, ch_len, ch_empty_count
        write (myfile, *)
        
        cur = ch_start
        
        do i = 1, ch_len
        
        
            write (myfile, *) "i",id(i),t(i),s(i),prev(i),next(i),ch_empty(i)
            write (myfile, *) "r",r(1, i),r(2, i),r(3, i)
            write (myfile, *) "v",v(1, i),v(2, i),v(3, i)
            write (myfile, *) "a",a(1, i),a(2, i),a(3, i)
            write (myfile, *)
            cur = next(cur)
            
        end do
        close (myfile)
    end subroutine

   ! Track Atom 1
    subroutine MPI_MD_debug_track_atom(track)
        integer n, track
        double precision rt1(3), rt2(3)        

        if (myid == 0) then 
            if (time_index == 0) then
                open (UNIT = (600), FILE=root_path(1:root_len)//"out/track.txt", STATUS = "replace")
            end if
            !call MPI_Barrier(MPI_COMM_WORLD, ierr)
            !print *, "receiving", time_index
            call MPI_RECV(rt1, 3, MPI_DOUBLE_PRECISION , MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(rt2, 3, MPI_DOUBLE_PRECISION , MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, ierr)
            write (600, *) time_index, rt1(1), rt1(2), rt1(3), rt2(1), rt2(2), rt2(3)
            if (time_index == end_time)then
                close (600)
            end if
        else 
            !call MPI_Barrier(MPI_COMM_WORLD, ierr)
            !n = ch_start
            do n = 1, NOS
                
                if ((id(n) .eq. track) .and. (s(n) .eq. 1)) then
                    !print *, "sending", time_index
                    call MPI_SEND(r0(1, n), 3, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(r1(1, n), 3, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
                end if
                !n = next(n)
            end do
        endif
    end subroutine
    
    !Track pinned atoms
    subroutine MPI_MD_debug_pinned_atoms()
        integer i, cur

        cur = ch_start
        do i = 1, ch_len
        if (s(cur) == 2) then
            print *, "pinned at", myid, id(cur)
        end if
        cur = next(cur)
        end do
    
    end subroutine
    
    !Dump Rshadow Data
    subroutine MPI_MD_debug_receive_data(from)
        integer myfile, i, from
        character*5 str_time
        character str_from
        character*30 str_node
        
        myfile = 1000 + mygpid + from*100
        
        write (str_node, "(I2)") mygpid
        write (str_time, "(I3)") time_index/sync_time
        write (str_from, "(I1)") from
        str_node = "fout_rshawdow."//str_time(1:3)//"."//str_from//"."//str_node(1:2)//".txt"
        open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/"//str_node(1:26), STATUS="replace")
        write (myfile, *) "Shadow List", myid, rshadow_count
        do i = 1, rshadow_count
            write (myfile, *) rshadow_t(i), rshadow_r0(1,i),  rshadow_r0(2,i),  rshadow_r0(3,i)
        end do
        write (myfile, *)
        write (myfile, *) "Atom List", ratom_count
        do i = 1, ratom_count
            write (myfile, *)ratom_t(i), ratom_s(i), ratom_id(i)
            write (myfile, *)ratom_r0(1, i), ratom_r0(2, i), ratom_r0(3, i)
            write (myfile, *)ratom_r1(1, i), ratom_r1(2, i), ratom_r1(3, i)
            write (myfile, *)ratom_r2(1, i), ratom_r2(2, i), ratom_r2(3, i)
            write (myfile, *)ratom_r3(1, i), ratom_r3(2, i), ratom_r3(3, i)
            write (myfile, *)ratom_r4(1, i), ratom_r4(2, i), ratom_r4(3, i)
            write (myfile, *)ratom_r5(1, i), ratom_r5(2, i), ratom_r5(3, i)

            write (myfile, *)
        end do
        close (myfile)
    end subroutine
    
    !Dump EAM parameters
    subroutine MPI_MD_debug_EAM()
        integer myfile, i1, i2
        character*30 str_node
        
        myfile = 50000 + myid
        
        write (str_node, "(I2)") myid
        str_node = "fout_EAM"//str_node(1:2)//".csv"
        open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/"//str_node(1:14), STATUS="replace")
        
        write (myfile, *) NOP_Spline
        write (myfile, *)
        
        do i1 = 1, 3
            write (myfile, *) "ID=", i1
            write (myfile, *)
            write (myfile, *) "Max, Min", v_max(i1),",", v_min(i1)
            write (myfile, *)
            write (myfile, *) "Pair Potential"
            do i2 = 1, NOP_Spline
                write (myfile, *) v_inc(i1)*(i2-1) + v_min(i1), ",",pot_V(i1, i2, 1),",", pot_V(i1, i2, 2)
            end do
            write (myfile, *)
            write (myfile, *) "Diff Pair Force"
            do i2 = 1, NOP_Spline
                write (myfile, *) v_inc(i1)*(i2-1) + v_min(i1),",", for_V(i1, i2, 1),",", for_V(i1, i2, 2)
            end do
        end do
        
        do i1 = 1, 2
            write (myfile, *) "ID=", i1
            write (myfile, *)
            write (myfile, *) "Max", f_max(i1)
            write (myfile, *)
            write (myfile, *) "Embedding Energy"
            do i2 = 1, NOP_Spline
                write (myfile, *) f_inc(i1)*(i2-1),",", pot_F(i1, i2, 1),",", pot_F(i1, i2, 2)
            end do
            write (myfile, *)
            write (myfile, *) "Diff Embedding Energy"
            do i2 = 1, NOP_Spline
                write (myfile, *) f_inc(i1)*(i2-1),",", for_F(i1, i2, 1),",", for_F(i1, i2, 2)
            end do
        end do
        
        do i1 = 1, 2
            write (myfile, *) "ID=", i1
            write (myfile, *)
            write (myfile, *) "Max, min", rho_max(i1), rho_min(i1)
            write (myfile, *)
            write (myfile, *) "Rho"
            do i2 = 1, NOP_Spline
                write (myfile, *) rho_inc(i1)*(i2-1)+rho_min(i1),",", rho(i1, i2, 1),",", rho(i1, i2, 2)
            end do
  
        end do
        
        
    end subroutine
        
    !Dump Sshadow Data
    subroutine MPI_MD_debug_send_data(from)
        integer myfile, i, from, si
        character*5 str_time
        character str_from
        character*30 str_node
        
        myfile = 20000 + mygpid + from*10 + time_index*100
        
        write (str_node, "(I2)") mygpid
        write (str_time, "(I3)") time_index/sync_time
        write (str_from, "(I1)") from
        str_node = "fout_sshawdow."//str_time(1:3)//"."//str_from//"."//str_node(1:2)//".txt"
        open (UNIT=(myfile),FILE=root_path(1:root_len)//str_node(1:26), STATUS="replace")
        write (myfile, *) "Shadow List", sshadow_count
        do i = 1, sshadow_count
            si = sshadow_list(i) + 1
            write (myfile, *) "", id(si), r0(1, si), r0(2, si), r0(3, si)
        end do
        write (myfile, *)
        write (myfile, *) "Atom List", satom_count
        do i = 1, satom_count
            si = satom_list(i) + 1
            write (myfile, *) "", id(si), r0(1, si), r0(2, si), r0(3, si)
        end do
        close (myfile)
    end subroutine
 
    !Dump all Atoms
    subroutine MPI_MD_debug_atoms()
        integer count_list((NON - 1) *2), count_sum
        integer send_count(2), count_index, buffer_index
        integer shadow_list(ch_len), atom_list(ch_len)        
        integer shadow_sum, atom_sum
        integer shadow_count, atom_count, send_size
        integer indexed_real_type, indexed_int_type
        integer i, j, cur
        double precision atom_r0(3,NOTA), atom_id(NOTA), shadow_r0(3,2*NOTA), shadow_id(2*NOTA)
        character, allocatable, dimension (:) :: buffer
        
        
        shadow_count = 0
        atom_count = 0
    
        if (myid == 0) then
            open (30000, FILE=root_path(1:root_len)//"out/atoms.csv")
            open (40000, FILE=root_path(1:root_len)//"out/shadow.csv")
            allocate (buffer(NOTA*3*(size_of_real*3+size_of_int)))
            
            count_sum = 1
            count_index = 1
            buffer_index = 1
            
            
            do i = 1, NON - 1
                call MPI_RECV(send_count, 2, MPI_INTEGER, i, i, MPI_COMM_WORLD, status, ierr)
                
                send_size = send_count(1) *(3*size_of_real + size_of_int) 
                
                call MPI_RECV(buffer(count_sum), send_size, MPI_PACKED, i, i, MPI_COMM_WORLD, status, ierr)
                
                
                
call MPI_Unpack(buffer,send_size,buffer_index,shadow_r0(1,shadow_sum+1),send_count(1)*3,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
call MPI_Unpack(buffer,send_size,buffer_index,shadow_id(shadow_sum),send_count(1),MPI_INTEGER,MPI_COMM_WORLD,ierr)
call MPI_Unpack(buffer,send_size,buffer_index,atom_r0(1,atom_sum+1), send_count(2)*3,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
call MPI_Unpack(buffer,send_size,buffer_index,atom_id(atom_sum),send_count(2),MPI_INTEGER,MPI_COMM_WORLD,ierr)
                count_list(count_index) = send_count(1)
                count_list(count_index+1) = send_count(2)
                count_sum = count_sum + (send_count(1)+send_count(2)) * (3*size_of_real + size_of_int)
                shadow_sum = shadow_sum + send_count(1)
                atom_sum = atom_sum + send_count(2)
                count_index = count_index + 2
            end do
            
            count_index = 1
            atom_sum = 0
            shadow_sum = 0
            do i = 1, NON - 1
                do j = 1, count_list(count_index)
                    atom_sum = atom_sum+1
                    write (30000, *) i, atom_id(atom_sum), atom_r0(1,atom_sum),atom_r0(2,atom_sum),atom_r0(3,atom_sum)
                end do
                do j = 1, count_list(count_index+1)
                    shadow_sum = shadow_sum + 1
                    write (40000,*)i, shadow_id(shadow_sum),shadow_r0(1,shadow_sum),shadow_r0(2,shadow_sum),shadow_r0(3,shadow_sum)
                    
                end do
            
            end do
            
        else
            
            cur = ch_start
            do i = 1, ch_len
                if (s(cur) == 0) then
                    shadow_count = shadow_count + 1
                    shadow_list(shadow_count) = cur - 1
                    
                else
                    atom_count = atom_count + 1 
                    atom_list(atom_count) = cur - 1
                    
                end if
                
                cur = next(cur)
            
            end do
        
            buffer_index = 1
            send_count(1) = atom_count
            send_count(2) = shadow_count   
            send_size = (send_count(1) + send_count(2))*(3*size_of_real + size_of_int) 
            
            call MPI_SEND(send_count, 2, MPI_INTEGER, 0, myid, MPI_COMM_WORLD, ierr)
            allocate (buffer(send_size))
        
            call MPI_TYPE_indexed (shadow_count, block_size3, 3*shadow_list, MPI_DOUBLE_PRECISION, indexed_real_type, ierr)
            call MPI_TYPE_COMMIT(indexed_real_type, ierr)
                        
            call MPI_TYPE_indexed (shadow_count, block_size1, shadow_list, MPI_INTEGER, indexed_int_type, ierr)
            call MPI_TYPE_COMMIT(indexed_int_type, ierr)
            
            call MPI_Pack(r0(1,1), 1 , indexed_real_type, buffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
            call MPI_Pack(id(1),   1, indexed_int_type,  buffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
            
            call MPI_TYPE_FREE(indexed_real_type, ierr)
            call MPI_TYPE_FREE(indexed_int_type, ierr)
            
            call MPI_TYPE_indexed (satom_count, block_size3, 3*satom_list, MPI_DOUBLE_PRECISION, indexed_real_type, ierr)
            call MPI_TYPE_COMMIT(indexed_real_type, ierr)
                    
            call MPI_TYPE_indexed (satom_count, block_size1, satom_list, MPI_INTEGER, indexed_int_type, ierr)
            call MPI_TYPE_COMMIT(indexed_int_type, ierr)
                    
            call MPI_Pack(r0(1,1), 1, indexed_real_type,  buffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
            call MPI_Pack(id(1),  1, indexed_int_type,   buffer, send_size, buffer_index, MPI_COMM_WORLD, ierr)
            
            call MPI_TYPE_FREE(indexed_real_type, ierr)
            call MPI_TYPE_FREE(indexed_int_type, ierr)
            
            call MPI_SEND(buffer, send_size, MPI_PACKED, 0, myid, MPI_COMM_WORLD, ierr)
            
            deallocate(buffer)
            
  
        end if 
        
    end subroutine
    
    !Debug Dump atom to screen
    subroutine MPI_MD_debug_atoms_to_screen()
        integer i, cur
        
        if (myid .ne. 0) then
            cur = ch_start
            do i = 1, ch_len
                print *, "",myid, id(cur), s(cur), r0(1,cur),r0(2,cur),r0(3,cur)
                cur = next(cur)
            end do
    
        end if 
        
    end subroutine
 
    !Debug Nlist to screen
    subroutine MPI_MD_debug_nlist_to_screen()
        integer xp, xn, yp, yn, zp,zn, i, j, k, n1, n2
        integer cur1, cur2
        integer ix(3)
        double precision dx(3)
        double precision start_time1, stop_time1, elasped_time1
        double precision start_time2, stop_time2, elasped_time2
        !logical debug1(10,10) /100*.false./
        !logical debug2(10,10) /100*.false./
        
        if (myid .ne. 0) then
        
        call MPI_MD_nlist_renew()
        
        print *, "============================="
        cur1 = ch_start
        
        start_time1 = MPI_wtime()
        do n1 =1, ch_len
            ix(:) = to_int3( (r0(:, cur1) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))
            
            !Set the upper bound and lower bound of the index
            if (ix(1) > 0) then
                xn = ix(1)-1
            else 
                xn = ix(1)
            end if
            
            if (ix(1) < NOD(1)) then
                xp = ix(1)+1
            else
                xp = ix(1)
            end if
            
            if (ix(2) > 0) then
                yn = ix(2)-1
            else 
                yn = ix(2)
            end if
            
            if (ix(2) < NOD(2)) then
                yp = ix(2)+1
            else
                yp = ix(2)
            end if

            if (ix(3) > 0) then
                zn = ix(3)-1
            else 
                zn = ix(3)
            end if
            
            if (ix(3) < NOD(3)) then
                zp = ix(3)+1
            else
                zp = ix(3)
            end if

            do i = xn, xp
            do j = yn, yp
            do k = zn, zp
                if (nbox_time(i,j,k) == time_index) then
                    cur2 = nbox_start(i,j,k)
                    do n2 = 1, nbox_count(i,j,k)
                        dx(:) = r0(:,cur1)-r0(:,cur2)
                        
                        if ((dx(1)**2 + dx(2)**2 + dx(3)**2 < f_cut_sq) .and. (cur1 .ne. cur2))then
                            print *, time_index, myid, " nlist", cur1, cur2
                        end if
                        cur2 = nbox_next(cur2)
                    end do
                
                end if
                        
            
            end do
            end do
            end do
            cur1 = next(cur1)
        end do
        
        stop_time1 = MPI_Wtime()
        elasped_time1 = stop_time1 - start_time1
        
        print *, time_index, myid, " nlist Ends"
        
        start_time1 = MPI_wtime()
        
        cur1 = ch_start
        do n1 = 1, ch_len
        cur2 = ch_start
        do n2 = 1, ch_len
           dx(:) = r0(:,cur1)-r0(:,cur2)             
            if ((dx(1)**2 + dx(2)**2 + dx(3)**2 < f_cut_sq) .and. (cur1 .ne. cur2))then
                print *, time_index, myid, " nlist", cur1, cur2
            end if
            cur2 = next(cur2)
        end do
        cur1 = next(cur1)
        
        end do
        
        stop_time2 = MPI_Wtime()
        elasped_time2 = stop_time2 - start_time2

        
        
        print *, time_index, myid, "nlist Test Ends", stop_time1, stop_time2
        print *, "============================="
        
        end if
    end subroutine
 
    !Print EAM Spline
    subroutine MPI_MD_debug_EAM_spline()
        integer myfile, div1, i1
        double precision drt1, end1, inc1
        
        end1 = 6.690d0
        div1 = 10000
        inc1 = end1/dble(div1)
        myfile = 70000
        
        if (myid .eq. 1) then
            open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/EAM_spline.csv", STATUS="replace")
            
            write (myfile, *) 1
            
            do i1 = 1, div1
                drt1 = inc1*(i1 -1)
                write (myfile, *) drt1,",", get_f(1, drt1),",", get_df(1, drt1)
            end do
        end if
    end subroutine

    !Dump atoms with specific ID    
    subroutine MPI_MD_debug_dump_data2()
        integer i1, cur
        integer myfile
        character*20 str_node
        
        myfile = 1500 + myid
        
        write (str_node, "(I2)") myid
        str_node = "cut_atoms"//str_node(1:2)//".txt"
        open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/"//str_node(1:15), STATUS="replace")
                
        cur = ch_start
        
        do i1 = 1, ch_len
            if ((r0(1, cur) > 349) .and. (r0(1, cur) < 430)) then 
            if ((r0(2, cur) > 288) .and. (r0(2, cur) < 323)) then
                write (myfile, *) r0(1, cur), r0(2, cur), r0(3, cur), s(cur), t(cur)
            end if
            end if          

            cur = next(cur)
        end do
        
    end subroutine

    !Dump all atoms at different stage of frame 2
    subroutine MPI_MD_debug_dump_data3(stage)
        integer i1, cur
        integer myfile, stage
        character*30 str_node
        character*200 mystring
        double precision r_p(3)
        double precision mybox_width_ex(3)
        
        myfile = 1500 + myid
       
        write (str_node, "(I2.2,I3.3)") stage, myid
        str_node = "dump_all_atoms"//str_node(1:5)//".cfg"
        open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/"//str_node(1:23), STATUS="replace")
                
        mybox_width_ex(:) = sr(2,:) - sr(1,:) + 5d0*f_cut_vect(:)


        write (mystring, *) "Number of particles = ", ch_len
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "A = 1 Angstorm"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(1,1) = ", mybox_width_ex(1), " A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(1,2) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(1,3) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(2,1) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(2,2) = ", mybox_width_ex(2), " A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(2,3) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(3,1) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(3,2) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(3,3) = ", mybox_width_ex(3), " A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) ".NO_VELOCITY."
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "entry_count = 7"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[0] = PE [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[1] = ID [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[2] = Type [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[3] = Status [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))    
        write (mystring, *) "27.0"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "Al"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))


        cur = ch_start
           
        do i1 = 1, ch_len
            r_p(:) = (r0(:,i1)-sr(1,:)+1.5d0*f_cut_vect(:))/mybox_width_ex(:)

            write (mystring, "(4ES14.5,I8,2I3)") r_p(:), Energy(i1), id(i1), t(i1), s(i1)
            write (myfile, '(A)') TRIM(ADJUSTL(mystring))         

            cur = next(cur)
        end do


        
    end subroutine
    
    
#endif
    
    !Dump Sorted Atoms CFG
    subroutine MPI_Dump_Sorted_Atoms(stage)
        integer i1, cur
        integer myfile, stage, ix(3)
        character*30 str_node
        character*200 mystring
        double precision r_p(3)
        double precision mybox_width_ex(3)
        
        myfile = 1500 + myid
       
        write (str_node, "(I5.5,I3.3)") stage, myid
        str_node = "dump_sorted_atoms"//str_node(1:8)//".cfg"
        open (UNIT=(myfile),FILE=root_path(1:root_len)//"out/"//str_node(1:29), STATUS="replace")
                
        mybox_width_ex(:) = sr(2,:) - sr(1,:) + 4d0*f_cut_vect(:)


        write (mystring, *) "Number of particles = ", ch_len
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "A = 1 Angstorm"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(1,1) = ", mybox_width_ex(1), " A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(1,2) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(1,3) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(2,1) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(2,2) = ", mybox_width_ex(2), " A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(2,3) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(3,1) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(3,2) = 0 A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "H0(3,3) = ", mybox_width_ex(3), " A"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) ".NO_VELOCITY."
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "entry_count = 11"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[0] = ID [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[1] = Type [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[2] = Status [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[3] = Index [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[4] = Current Index [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[5] = I [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[6] = J [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "auxiliary[7] = K [reduced unit]"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "27.0"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))
        write (mystring, *) "Al"
        write (myfile, '(A)') TRIM(ADJUSTL(mystring))


        cur = ch_start
           
        do i1 = 1, ch_len
            r_p(:) = (r0(:,cur)-sr(1,:)+ f_cut_vect_2(:))/mybox_width_ex(:)
            ix(:) = to_int3( (r0(:, cur) - br(myid, 1, :) + f_cut_vect_2(:))/nbox_width(:))

            write (mystring, "(3ES14.5,I8,I3, I5, 2I8, 3I6)") r_p(:), id(cur), t(cur), s(cur), i1, cur, ix(:)
            write (myfile, '(A)') TRIM(ADJUSTL(mystring))         

            cur = next(cur)
        end do

    end subroutine 
    
    !General Stage Debuging
    subroutine debug_print_stage(index1)
        integer index1, temp_value
        character*100 mystring
        temp_value = myid + 100*index1 + 10000*time_index
        
        write (mystring, *) "Stage:" , index1 , "Time:", time_index  , "ID:" , myid
        print *, mystring
        stage_index = index1 + 1
    end subroutine
    
    subroutine debug_fprint_stage(index1)
        integer index1, temp_value, myfile1
        character*100 mystring
        character*100 str_node
      
        temp_value = myid + 100*index1 + 10000*time_index
                        
        if (myid .ne. 0) then
            myfile1 = 17000 +  myid + time_index + stage_index
            write (mystring, "(A6,I3.3,A8,I7.7,A5,I3.3)" ) "Stage:" , index1 , "  Time:", time_index  , "  ID:" , myid            
            
            print *, mystring
            write (mystring, "(I7.7)") temp_value
            str_node = "Debug_"//mystring(1:32)//".cfg"
            open (UNIT=myfile1,FILE=root_path(1:root_len)//"out/"//mystring(1:7), STATUS="replace") !root_path(1:root_len)//"out/"//
            close (myfile1)
       end if

       stage_index = index1 + 1
    end subroutine
   
    !Find Atom
    subroutine debug_find_my_atom()
	integer i1, cur
	    cur = ch_start
	    do i1 = 1, ch_len
		    if ((id(cur) == 390231) .or. (id(cur) == 400431) .or. (id(cur) == 405150)) then
			    print *, "My ID", myid, time_index, id(cur), cur
		    end if
		    cur = next(cur)
	    end do
	end subroutine

    !Dump all Atoms
    subroutine MPI_MD_write_debug_cfg()
        integer i, j, cur
        integer Stripe(0:NON-1)
        integer TNOAA(0:NON-1)
        integer buffer_index
        integer buffer_size, buffer_size_sum
        integer rbuffer_size(0:NON-1)
        integer myfile1
        integer wid, wt, ws        
        double precision wr0(3), wr1(3), wr2(3), wr3(3), wr4(3), wr5(3)
        double precision T_ave, wPE
        double precision r_p(3)
        double precision mybox_width_ex(3)
        character*30 str_node
        character*500 mystring
        character, allocatable, dimension (:) :: buffer
        
                
        if (myid .ne. 0) then
            myfile1 = 25000 +  myid + time_index
            
            write (str_node, "(I2)") myid
            str_node = "Debug_"//str_node(1:2)//".cfg"
            open (UNIT=myfile1,FILE=root_path(1:root_len)//"out/"//str_node(1:16), STATUS="replace")
                      
            
            mybox_width_ex(:) = sr(2,:) - sr(1,:) + 2d0 * f_cut_vect_2(:)

            write (mystring, *) "Number of particles = ", ch_len
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "A = 1 Angstorm"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,1) = ", mybox_width_ex(1), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,2) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(1,3) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,1) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,2) = ", mybox_width_ex(2), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(2,3) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,1) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,2) = 0 A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "H0(3,3) = ", mybox_width_ex(3), " A"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) ".NO_VELOCITY."
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "entry_count = 19"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[0] = R1_X [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[1] = R1-Y [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))  
            write (mystring, *) "auxiliary[2] = R1-Z [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[3] = R2-X [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[4] = R2-Y [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))  
            write (mystring, *) "auxiliary[5] = R2-Z [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[6] = dR2-X [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[7] = dR2-Y [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))  
            write (mystring, *) "auxiliary[8] = dR2-Z [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[9] = dR2_ex-X [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[10] = dR2_ex-Y [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))  
            write (mystring, *) "auxiliary[11] = dR2_ex-Z [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
                      
            write (mystring, *) "auxiliary[12] = PE [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[13] = ID [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[14] = Type [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "auxiliary[15] = Status [reduced unit]"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))    
            write (mystring, *) "27.0"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
            write (mystring, *) "Al"
            write (myfile1, '(A)') TRIM(ADJUSTL(mystring))            

            cur = ch_start
            
            do j = 1, ch_len
                r_p(:) = (r0(:,cur)- sr(1,:) +f_cut_vect_2(:)) / mybox_width_ex(:)
                       
write (mystring, "(16ES14.5,I8,2I4)") r_p(:),r1(:,cur), r2(:,cur), dr2(:,cur), dr2_ex(:,cur), Energy(cur), id(cur), t(cur), s(cur)
                write (myfile1, '(A)') TRIM(ADJUSTL(mystring))
                cur = next(cur)
            end do
           
            close (myfile1)
        end if
        
  
        
    end subroutine      
    
    

    
end module

