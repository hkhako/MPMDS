!Version 1.0.0 Sept 20, 2010

program md_mpi
#include "PRESET.INI"
    use function_library
    implicit none

   
    call MPI_MD_init()
   
    do while (MPI_MD_designate_instruction())

        call MPI_MD_distribute_parameters()

        call MPI_MD_data_init()

        call MPI_MD_distribute_border()

        call MPI_MD_distribute_atoms()

        call MPI_MD_force_mode_init()

        call MPI_MD_distribute_EAM()

        call MPI_MD_nlist_init_XP()

        init_time = MPI_Wtime()

        do time_index = 0, end_time
            if (myid .ne. 0) then

              call MD_remove_shadow()

              call MD_predict()

              call MPI_MD_redist_atoms_Linear_XP()

              call MD_atom_folding_XP()

              call MD_atom_chain_sorting()
              
              call MPI_MD_nlist_renew()

            if (write_stress .and. (mod(time_index, movie_time*sync_time) == 0)) then
                  
                  call MD_EAM_force_stress_XP()
                  
                  call MD_Stress_Average()
                  
            else
                  
                  call MD_EAM_force_XP()

            end if         
              
              call MD_correct()

              call MD_temperature_control()

            end if
            
#ifdef Stress_Strain                  
            call MPI_MD_MEASURE_elastic()
#endif    
            
            !Recording
            if (mod(time_index, sync_time) == 0) then
            
                !call MPI_Barrier(MPI_COMM_WORLD, ierr)
           
                call MPI_MD_collect_atoms_binary()
           
                if (mod(time_index, movie_time*sync_time) == 0) then
            
                    if (write_stress) then
                        
                        call MPI_MD_write_stress_cfg()
                        
                    else
                    
                        call MPI_MD_write_cfg()
                    
                    end if
            
                end if
            end if
                 
        end do
    
        call MPI_MD_write_input()
        
        call MPI_MD_record_time()
        
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
        call MPI_MD_finalize()
    
    end do
    
    !call MPI_Finalize()
    
    
  
end program
    

