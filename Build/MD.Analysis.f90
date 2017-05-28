module common_variable
    implicit none
    include 'mpif.h'
    
    double precision, allocatable, dimension (:, :, :, :), save :: r
    double precision, allocatable, dimension (:), save :: time, TEMP
    double precision, allocatable, dimension (:, :), save :: PE
    double precision, allocatable, dimension (:), save :: E_bin_max, E_bin_min
    integer, allocatable, dimension (:), save :: t, s
    integer end_time, sync_time, movie_time, NOTA, myid, size
    integer time_index, n, time_span, time1
    integer NOF  !number of frame
    integer frame(100000)
    double precision dt
    double precision :: dwidth(3)
    data dwidth /90.560753088741482704571533583617, 90.560753088741482704571533583617, 81/
    double precision :: a
    data a /4.05/
    data time1 /10000/
    integer, parameter :: NON = 9
    integer, parameter :: bin_size = 1001
    double precision :: E_select1, E_select2
    logical write_formatted
end module


program md_analysis
    use common_variable
    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) 
    call MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr )
    call read_data_binary()
    call write_last_frame()
    call write_format_xyz_traj_series()
    !call generate_POV_RAY_script_binary()
    !call energy_spectrum()
    !call write_format_xyz()
    !call write_formatted_frame()
    
    !call write_frame()
    !call generate_pair_correlation()
end program

subroutine energy_bound(f1)
    use common_variable
    integer myfile, i, f1,index, cf  
    double precision E_min, E_max
    
        cf = frame(f1)
    
        E_min = PE(1, cf)
        E_max = PE(1, cf)
             

        do i = 1, NOTA
            if (E_min > PE(i, cf)) then
                E_min = PE(i, cf)
            end if
            
            if (E_max < PE(i, cf)) then
                E_max = PE(i, cf)
            end if
        
        end do
        
        E_max = E_max*0.99d0
        E_bin_max(f1) = E_max
        E_bin_min(f1) = E_min
end subroutine

subroutine energy_spectrum()
    use common_variable
    integer myfile, i, f1,index, cf
    character*21 str_node
    
    double precision E_min, E_max
    integer*8 E_bin(bin_size)
    double precision bin_inc
    
    print *, "Printing Energy Spectrum"
    
    
    do f1 = 1, NOF
        cf = frame(f1)    !current frame
        E_min = PE(1, cf)
        E_max = PE(1, cf)
        
        do i = 1, bin_size
            E_bin(i) = 0
        end do
        

        do i = 1, NOTA
            if (E_min > PE(i, cf)) then
                E_min = PE(i, cf)
            end if
            
            if (E_max < PE(i, cf)) then
                E_max = PE(i, cf)
            end if
        
        end do
        
        E_max = E_max*0.99d0
        E_bin_max(f1) = E_max
        E_bin_min(f1) = E_min
        bin_inc = (E_max - E_min)/bin_size
        
        do i = 1, NOTA
            index = ceiling((PE(i, cf)-E_min)/bin_inc)
            E_bin(index) = E_bin(index) + 1
        end do
    
    
    
        write (str_node, "(I4.4)") cf
        
        str_node = "output_energy"//str_node(1:4)//".csv"
        myfile = i+100
        open (UNIT=(myfile),FILE="out/"//str_node(1:21), STATUS="replace")
        
        do i = 1, bin_size
            write (myfile, *) E_min+(i-1)*bin_inc, ",", E_bin(i)
        end do
        
        close (myfile)
    end do
    
    print *, "Finished Printing Energy Spectrum"

end subroutine

subroutine energy_spectrum_current(f1)
    use common_variable
    integer myfile, i, f1,index, cf
    character*21 str_node
    
    double precision E_min, E_max
    integer*8 E_bin(bin_size)
    double precision bin_inc
    
    
        cf = frame(f1)    !current frame
        E_min = PE(1, cf)
        E_max = PE(1, cf)
        
        do i = 1, bin_size
            E_bin(i) = 0
        end do
        

        do i = 1, NOTA
            if (E_min > PE(i, cf)) then
                E_min = PE(i, cf)
            end if
            
            if (E_max < PE(i, cf)) then
                E_max = PE(i, cf)
            end if
        
        end do
        
        E_max = E_max*0.99d0
        E_bin_max(f1) = E_max
        E_bin_min(f1) = E_min
        bin_inc = (E_max - E_min)/bin_size
        
        do i = 1, NOTA
            index = ceiling((PE(i, cf)-E_min)/bin_inc)
            E_bin(index) = E_bin(index) + 1
        end do
    
    
    
        write (str_node, "(I4.4)") cf
        
        str_node = "output_energy"//str_node(1:4)//".csv"
        myfile = i+100
        open (UNIT=(myfile),FILE="out/"//str_node(1:21), STATUS="replace")
        
        do i = 1, bin_size
            write (myfile, *) E_min+(i-1)*bin_inc, ",", E_bin(i)
        end do
        
        close (myfile)
    
    
    print *, "Finished Printing Energy Spectrum", cf

end subroutine

subroutine write_format_xyz()
    use common_variable
    integer myfile, i, j, cf, NOTA_selected
    character*20 str_node
    double precision E_limit1, E_limit2
    
    do i = 1, NOF
        cf = frame(i)    !current frame
        NOTA_selected = 0
        
        write (str_node, "(I4.4)") cf
        str_node = "output_atoms"//str_node(1:4)//".xyz"
        myfile = i+100
        open (UNIT=(myfile),FILE="out/"//str_node(1:20), STATUS="replace")
        
        
        E_limit1 = (E_bin_max(i) - E_bin_min(i))*(E_select1) + E_bin_min(i)
        E_limit2 = (E_bin_max(i) - E_bin_min(i))*(E_select2) + E_bin_min(i)
        
        
        do j = 1, NOTA
            if ((PE(j, cf) > E_limit1) .and. (PE(j, cf) < E_limit2)) then
                NOTA_selected = NOTA_selected + 1
            end if
        end do
        
        print *, cf, NOTA_selected
        
        write (myfile, *) NOTA_selected
        write (myfile, "(A3,ES14.5,A3,ES14.5)") "E1 ", E_select1, "E2 ", E_select2
        

        
        do j = 1, NOTA
            if ((PE(j, cf) > E_limit1) .and. (PE(j, cf) < E_limit2)) then
                if (t(j) == 1) then
                    write (myfile, "(A4,4ES14.5)") "Al", r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf), PE(j, cf)
                else if (t(j) == 2) then
                    write (myfile, "(A4,4ES14.5)") "Mg", r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf), PE(j, cf)
                end if    
            end if    
        end do
        
                

        close (myfile)
    end do


end subroutine

subroutine write_format_xyz_current(i)
    use common_variable
    integer myfile, i, j, cf, NOTA_selected
    character*20 str_node
    double precision E_limit1, E_limit2
    
    
        cf = frame(i)    !current frame
        NOTA_selected = 0
        print *, "a1", cf
        write (str_node, "(I4.4)") cf
        str_node = "output_atoms"//str_node(1:4)//".xyz"
        myfile = i+500
        open (UNIT=(myfile),FILE="out/"//str_node(1:20), STATUS="replace")
        print *, "a2"
        
        E_limit1 = (E_bin_max(i) - E_bin_min(i))*(E_select1) + E_bin_min(i)
        E_limit2 = (E_bin_max(i) - E_bin_min(i))*(E_select2) + E_bin_min(i)
        
        
        do j = 1, NOTA
            if ((PE(j, cf) > E_limit1) .and. (PE(j, cf) < E_limit2)) then
                NOTA_selected = NOTA_selected + 1
            end if
        end do
        print *, "a3"
        
        print *, cf, NOTA_selected
        
        write (myfile, *) NOTA_selected
        write (myfile, "(A3,ES14.5,A3,ES14.5)") "E1 ", E_select1, "E2 ", E_select2
        

        
        do j = 1, NOTA
            if ((PE(j, cf) > E_limit1) .and. (PE(j, cf) < E_limit2)) then
                if (t(j) == 1) then
                    write (myfile, "(A4,4ES14.5)") "Al", r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf), PE(j, cf)
                else if (t(j) == 2) then
                    write (myfile, "(A4,4ES14.5)") "Mg", r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf), PE(j, cf)
                end if    
            end if    
        end do
                

        close (myfile)
    


end subroutine

subroutine write_format_xyz_trajectory(i)
    use common_variable
    integer myfile, i, j, cf, NOTA_selected
    character*21 str_node
    double precision E_limit1, E_limit2
    
    
        cf = frame(i)    !current frame
        NOTA_selected = 0
        myfile = 1000
        !write (str_node, "(I4.4)") cf
        if (i .eq. 1) then
            str_node = "output_atoms_traj.xyz"
            
            open (UNIT=(myfile),FILE="out/"//str_node(1:21), STATUS="replace")
        end if
        
        
        E_limit1 = (E_bin_max(i) - E_bin_min(i))*(E_select1) + E_bin_min(i)
        E_limit2 = (E_bin_max(i) - E_bin_min(i))*(E_select2) + E_bin_min(i)
        
        
        do j = 1, NOTA
            if ((PE(j, cf) > E_limit1) .and. (PE(j, cf) < E_limit2)) then
                NOTA_selected = NOTA_selected + 1
            end if
        end do
        
        
        print *, cf, NOTA_selected
        
        write (myfile, *) NOTA_selected
        write (myfile, "(A3,ES14.5,A3,ES14.5)") "E1 ", E_select1, "E2 ", E_select2
        

        
        do j = 1, NOTA
            if ((PE(j, cf) > E_limit1) .and. (PE(j, cf) < E_limit2)) then
                if (t(j) == 1) then
                    write (myfile, "(A4,4ES14.5)") "Al", r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf), PE(j, cf)
                else if (t(j) == 2) then
                    write (myfile, "(A4,4ES14.5)") "Mg", r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf), PE(j, cf)
                end if    
            end if    
        end do
                
        if (i .eq. NOF) then
            close (myfile)
        end if
    


end subroutine

subroutine write_format_xyz_traj_series()
    use common_variable
    integer myfile, i, j, cf, NOTA_selected
    character*21 str_node
    double precision E_limit1, E_limit2
    
    
    do i = 1, NOF
        call energy_bound(i)
        cf = frame(i)    !current frame
        NOTA_selected = 0
        myfile = 1000
        !write (str_node, "(I4.4)") cf
        if (i .eq. 1) then
            str_node = "output_atoms_traj.xyz"
            
            open (UNIT=(myfile),FILE="out/"//str_node(1:21), STATUS="replace")
        end if
        
        
        E_limit1 = (E_bin_max(i) - E_bin_min(i))*(E_select1) + E_bin_min(i)
        E_limit2 = (E_bin_max(i) - E_bin_min(i))*(E_select2) + E_bin_min(i)
        
        
        do j = 1, NOTA
            if ((PE(j, cf) > E_limit1) .and. (PE(j, cf) < E_limit2)) then
                NOTA_selected = NOTA_selected + 1
            end if
        end do
        
        
        print *, cf, NOTA_selected
        
        write (myfile, *) NOTA_selected
        write (myfile, "(A3,ES14.5,A3,ES14.5)") "E1 ", E_select1, "E2 ", E_select2
        

        
        do j = 1, NOTA
            if ((PE(j, cf) > E_limit1) .and. (PE(j, cf) < E_limit2)) then
                if (t(j) == 1) then
                    write (myfile, "(A4,4ES14.5)") "Al", r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf), PE(j, cf)
                else if (t(j) == 2) then
                    write (myfile, "(A4,4ES14.5)") "Mg", r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf), PE(j, cf)
                end if    
            end if    
        end do
                
        if (i .eq. NOF) then
            close (myfile)
        end if
    
    end do

end subroutine

subroutine write_last_frame()
    use common_variable
    integer myfile, i, j, cf, NOTA_selected
    character*21 str_node
    double precision E_limit1, E_limit2
    
    
    
    i = NOF

    cf = frame(i)    !current frame
    NOTA_selected = 0
    myfile = 1000
        !write (str_node, "(I4.4)") cf
    str_node = "input_atom.txt"
            
    open (UNIT=(myfile),FILE="out/"//str_node(1:14), STATUS="replace")
    
    write (myfile, *) NOTA
    
    do j = 1, NOTA
        write (myfile, "(I9, 2I4,18ES14.5)") j, t(j), s(j), r(1,0,j,cf), r(2,0,j,cf), r(3,0,j,cf), &
              &r(1,1,j,cf), r(2,1,j,cf), r(3,1,j,cf), r(1,2,j,cf), r(2,2,j,cf), r(3,2,j,cf), &
              &r(1,3,j,cf), r(2,3,j,cf), r(3,3,j,cf), r(1,4,j,cf), r(2,4,j,cf), r(3,4,j,cf), &
              &r(1,5,j,cf), r(2,5,j,cf), r(3,5,j,cf)
    end do
    
    close (myfile)
        
  
end subroutine

subroutine generate_POV_RAY_script()
    use common_variable
    integer :: infile1
    integer :: outfile1
    double precision PEi
    double precision r0(3), r1(3), r2(3), r3(3), r4(3), r5(3)
    integer  id1, t1, s1, ts1
    integer i
    character*10 fname1
    character*14 fname2
    
    open (UNIT=infile1,FILE="out/atoms.out", STATUS="old")
    open (UNIT=outfile1,FILE="out/animate.pov", STATUS="replace")
    
    infile1 = 100
    outfile1 = 101
    
    
    read  (infile1, "(1I10)")       end_time
    read  (infile1, "(1I10)")       sync_time
    read  (infile1, "(1E20.14)")    dt
    read  (infile1, "(1I10)")       NOTA
    read  (infile1, *) 
    print *, end_time, sync_time, dt, NOTA
    
    time_span = end_time/sync_time

    allocate (r(3, 0:5, NOTA, 0:time_span), time(0:time_span), TEMP(0:time_span), t(NOTA), s(NOTA), PE(NOTA, 0:time_span))
    
    do time_index = 0, time_span
        read (infile1,"(1I10)") ts1 !, time(time_index)
        
        do n = 1, NOTA
            read (infile1, *) id1, PEi, r0(1), r0(2), r0(3), r1(1), r1(2), r1(3), r2(1), r2(2), r2(3), &
                                & r3(1), r3(2), r3(3), r4(1), r4(2), r4(3), r5(1), r5(2), r5(3), t1, s1
            !print *, time_index, n
            t(id1) = t1
            s(id1) = s1
            r(:, 0, id1, time_index) = r0(:)
            r(:, 1, id1, time_index) = r1(:)
            r(:, 2, id1, time_index) = r2(:)
            r(:, 3, id1, time_index) = r3(:)
            r(:, 4, id1, time_index) = r4(:)
            r(:, 5, id1, time_index) = r5(:)
            PE(id1, time_index) = PEi
            
        end do
        read (infile1,"(1A20, 1E20.10)")  TEMP(time_index)
        read (infile1, *)

    end do



    write (outfile1,*) "#version 3.6;"
    write (outfile1,*) "global_settings {  assumed_gamma 1.0 }"

    write (outfile1,*) "#include ", '"', "colors.inc", '"'

    write (outfile1,*) "#include ", '"', "textures.inc", '"'

    write (outfile1,*) "#include ", '"', "math.inc", '"'
    write (outfile1,*) "#include ", '"', "transforms.inc", '"'

    write (outfile1,*) "#declare Camera_1 = camera { "
    write (outfile1,*) "                             angle 15"
    write (outfile1,*) "                             location  <0.0 , 50.0 ,330.0>"
    write (outfile1,*) "                             right     x*image_width/image_height"
    write (outfile1,*) "                             look_at   <0.0 , 25.0 , 0.0>"
    write (outfile1,*) "                           }"
    write (outfile1,*) "camera{Camera_1}"

    write (outfile1,*) "light_source{<1000,2500,2500> color White}"
    write (outfile1,*) "light_source{<-1000,-2500,-2500> color White}"


    do n = 1, NOTA
        write (fname1, "(1I0.1)") n
        fname2 = "spl_"//fname1


        write (outfile1,"(1A9,1A14, 1A3)") "#declare ", fname2, " = "
        write (outfile1,*) "  spline {"
        write (outfile1,*) "    natural_spline"

        do time_index = 0, time_span-1
        
        write (outfile1, "(1I9, 1A3, 1e15.6,1A2, 1e15.6,1A2, 1e15.6, 1A2)") time_index, ", <", &
        & r(1,0, n,time_index), ", ",r(2,0, n,time_index), ", ", r(3, 0, n,time_index),">,"
        
        end do
        
        write (outfile1, "(1I9, 1A3, 1e15.6,1A2, 1e15.6,1A2, 1e15.6, 1A2)" ) time_index, ", <", &
        & r(1, 0, n,time_span), ", ",r(2,0, n,time_span), ", ", r(3,0, n,time_span),"> "
        write (outfile1, *) "  }"        
    end do
 
 
    write (outfile1, *)
    write (outfile1, *) "#declare Atom1 ="
    write (outfile1, *) "sphere { <0,0,0>, 0.25 "
    write (outfile1, *) "         texture { pigment{ color rgb<1,0,0>}"
    write (outfile1, *) "                   finish { ambient 0.1 diffuse 0.85  phong 1}"
    write (outfile1, *) "                 }"
    write (outfile1, *) "          scale<2.86,2.86,2.86>"
    write (outfile1, *) "       }"
    
    write (outfile1, *)
    write (outfile1, *) "#declare Atom2 ="
    write (outfile1, *) "sphere { <0,0,0>, 0.25 "
    write (outfile1, *) "         texture { pigment{ color rgb<0,1,0>}"
    write (outfile1, *) "                   finish { ambient 0.1 diffuse 0.85  phong 1}"
    write (outfile1, *) "                 }"
    write (outfile1, *) "          scale<3,3,3>"
    write (outfile1, *) "       }"
    
    write (outfile1, *)
    write (outfile1, *) "#declare Atom3 ="
    write (outfile1, *) "sphere { <0,0,0>, 0.25 "
    write (outfile1, *) "         texture { pigment{ color rgb<0,0,1>}"
    write (outfile1, *) "                   finish { ambient 0.1 diffuse 0.85  phong 1}"
    write (outfile1, *) "                 }"
    write (outfile1, *) "          scale<3,3,3>"
    write (outfile1, *) "       }"
    
    do n = 1, NOTA
        write (fname1, "(1I0.1)") n
        fname2 = "spl_"//fname1
        
        if (t(n) == 1) then 
            if (s(n) == 2) then
                write (outfile1,*) "object { Atom3"
            else 
                write (outfile1,*) "object { Atom1"
            end if
        else if (t(n) == 2) then
            if (s(n) == 2) then
                write (outfile1,*) "object { Atom3"
            else 
                write (outfile1,*) "object { Atom2"
            end if
            
        end if
        
        write (outfile1,*) "translate ", fname2, "(clock)"
        write (outfile1,*) "}"       
    end do

    
end subroutine

subroutine read_data_binary()
    use common_variable
    integer :: infile1
    integer :: outfile1
    integer :: infile2 
    double precision PEi
    double precision r0(3), r1(3), r2(3), r3(3), r4(3), r5(3)
    double precision time_index_dt
    double precision T_current, T_ave
    integer time_index_data
    integer TNOAA(0:NON-1)
    integer  id1, t1, s1, ts1
    integer i, j, i1, ierr
    integer Stripe(0:NON-1)
    integer buffer_index, buffer_size_sum
    integer rbuffer_size(0:NON-1)
    character*10 fname1
    character*14 fname2
    integer, allocatable, dimension (:) :: id_track
    character, allocatable, dimension (:) :: buffer
    CHARACTER(LEN=600) :: FMT1 = "(1I15, 19E25.15, 2I5)"
    
    
     
    
    infile1 = 100
    infile2 = 101
    outfile1 = 102
    
    i1 = 1
    
    open (UNIT=infile1,FILE="out/atoms.out", STATUS="old", FORM="unformatted")
    open (UNIT=infile2,FILE="in/analyze.in", STATUS="old", FORM="formatted")
    !if (write_formatted) open (UNIT=outfile1,FILE="out/atoms_formated.out", STATUS="replace")
    
    read (infile2, *)      E_select1, E_select2
    read (infile2, *)      NOF
    do i = 1, NOF
        read (infile2, *) frame(i)
    end do
    
    allocate (E_bin_max(NOF), E_bin_min(NOF)) 
    
    read  (infile1)       end_time
    read  (infile1)       sync_time
    read  (infile1)       movie_time
    read  (infile1)       dt
    read  (infile1)       NOTA


    !if (write_formatted) write (outfile2, "(1I10)") end_time
    !if (write_formatted) write (outfile2, "(1I10)") sync_time
    !if (write_formatted) write (outfile2, "(1I10)") movie_time
    !if (write_formatted) write (outfile2, "(1D20.14)") dt
    !if (write_formatted) write (outfile2, "(1I10)") NOTA
    !if (write_formatted) write (outfile2, *)
    
    time_span = nint(end_time/dble(sync_time*movie_time))

    print *, end_time, sync_time, movie_time, time_span, dt, NOTA

    allocate (r(3, 0:5, NOTA, 0:time_span), time(0:time_span), TEMP(0:time_span), t(NOTA), s(NOTA), &
        & PE(NOTA, 0:time_span), id_track(NOTA))

    do i = 1, NOTA
        id_track(i) = -1
    end do
    
    do time_index = 0, time_span
    
        print *, "reading ", time_index
        read (infile1) TNOAA
               
        read (infile1) Stripe, buffer_size_sum
                       
        Stripe(0) = 0
        
        T_ave = 0
        
        allocate (buffer(buffer_size_sum)) 
        
        read (infile1) time_index_data, time_index_dt
                
        !if (write_formatted) write (outfile2, "(1I10, 1E10.5)") time_index, time_index_dt
        
        read (infile1) buffer
        
        buffer_index = 1
        
        print *, "Buffer Starts"
        !print *, buffer
        !print *, "Buffer Ends"
        
        do i = 1, NON-1
            call MPI_Unpack(buffer, buffer_size_sum, buffer_index, T_current, 1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD, ierr)
            T_ave = T_ave + T_current * TNOAA(i)/ dble(NOTA)
            do j = 1, TNOAA(i)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, id1, 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, PEi, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, r0, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, r1, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, r2, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, r3, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, r4, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, r5, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, t1, 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                
                call MPI_Unpack(buffer, buffer_size_sum, buffer_index, s1, 1, MPI_INTEGER , MPI_COMM_WORLD, ierr)
                
                id_track(id1) = id_track(id1) + 1
                
                t(id1) = t1
                
                s(id1) = s1
                
                r(:, 0, id1, time_index) = r0(:)
                
                r(:, 1, id1, time_index) = r1(:)
                
                r(:, 2, id1, time_index) = r2(:)
                
                r(:, 3, id1, time_index) = r3(:)
                
                r(:, 4, id1, time_index) = r4(:)
                
                r(:, 5, id1, time_index) = r5(:)
                
                PE(id1, time_index) = PEi
                
                
                !if (write_formatted) write (outfile2, FMT1)     id1, PEi, r0(1), r0(2), r0(3), r1(1), r1(2), r1(3), &
                !                             & r2(1), r2(2), r2(3), r3(1), r3(2), r3(3), &
                !                             & r4(1), r4(2), r4(3), r5(1), r5(2), r5(3), t1, s1
                
                !print *, id1, PEi, r0(1), r0(2), r0(3), r1(1), r1(2), r1(3), &
                                             !& r2(1), r2(2), r2(3), r3(1), r3(2), r3(3), &
                                             !& r4(1), r4(2), r4(3), r5(1), r5(2), r5(3), t1, s1  
            
            end do
            
            buffer_index = buffer_index + 2
        end do
        print *, "time_index_data", time_index_data, "Buffer Size ", buffer_size_sum, "Temperature ", T_ave
        !if (write_formatted) write (outfile2, "(1E20.10)") T_ave
        !if (write_formatted) write (outfile2, *)
    
        deallocate(buffer)
        
        do i = 1, NOTA
            if (id_track(i) == (time_index - 1)) then
                print *, "Missing ", i, " at ", time_index
            end if
        end do
        
        
       ! if (time_index .eq. frame(i1)) then
       !     print *, "Writting ", frame(i1)
       !     call energy_bound(i1)
            !call energy_spectrum_current(i1)
            !call write_format_xyz_current(i1)
            !call write_format_xyz_trajectory(i1)
       !     i1 = i1 + 1
       ! end if
        

    end do
end subroutine

! not finished
subroutine write_formatted_frame()
    use common_variable
    integer myfile, i, j, cf
    character*20 str_node    


    do i = 1, NOF
        cf = frame(i)    !current frame
        
        write (str_node, "(I4.4)") cf
        str_node = "input_atom"//str_node(1:4)//".txt"
        myfile = i+100
        open (UNIT=(myfile),FILE="out/"//str_node(1:19), STATUS="replace")
        
        write (myfile, *) NOTA
        
        do j = 1, NOTA
            write (myfile, "(I10,1x,I2,1x,I2,1x)", advance="no") j, t(j), s(j)
            write (myfile, "(ES14.5, 1x,ES14.5, 1x,ES14.5, 1x)", advance="no") r(1, 0, j, cf), r(2, 0, j, cf), r(3, 0, j, cf)
            write (myfile, "(ES14.5, 1x,ES14.5, 1x,ES14.5, 1x)", advance="no") r(1, 1, j, cf), r(2, 1, j, cf), r(3, 1, j, cf)
            write (myfile, "(ES14.5, 1x,ES14.5, 1x,ES14.5, 1x)", advance="no") r(1, 2, j, cf), r(2, 2, j, cf), r(3, 2, j, cf)
            write (myfile, "(ES14.5, 1x,ES14.5, 1x,ES14.5, 1x)", advance="no") r(1, 3, j, cf), r(2, 3, j, cf), r(3, 3, j, cf)
            write (myfile, "(ES14.5, 1x,ES14.5, 1x,ES14.5, 1x)", advance="no") r(1, 4, j, cf), r(2, 4, j, cf), r(3, 4, j, cf)
            write (myfile, "(ES14.5, 1x,ES14.5, 1x,ES14.5, 1x)", advance="yes") r(1, 5, j, cf), r(2, 5, j, cf), r(3, 5, j, cf)

        end do

        close (myfile)
    end do



end subroutine

subroutine generate_POV_RAY_script_binary()
    use common_variable
    integer :: outfile1
    integer :: outfile2 
    double precision PEi
    double precision r0(3), r1(3), r2(3), r3(3), r4(3), r5(3)
    double precision time_index_dt
    double precision T_current, T_ave
    integer time_index_data
    integer TNOAA(0:NON-1)
    integer  id1, t1, s1, ts1
    integer i, j, ierr
    integer Stripe(0:NON-1)
    integer buffer_index, buffer_size_sum
    integer rbuffer_size(0:NON-1)
    character*10 fname1
    character*14 fname2
    integer, allocatable, dimension (:) :: id_track
    character, allocatable, dimension (:) :: buffer
    CHARACTER(LEN=600) :: FMT1 = "(1I15, 19E25.15, 2I5)"
    
    
     
    
    outfile1 = 101
        
    open (UNIT=outfile1,FILE="out/animate.pov", STATUS="replace")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Generate POVRay Here
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    write (outfile1,*) "#version 3.6;"
    write (outfile1,*) "global_settings {  assumed_gamma 1.0 }"

    write (outfile1,*) "#include ", '"', "colors.inc", '"'

    write (outfile1,*) "#include ", '"', "textures.inc", '"'

    write (outfile1,*) "#include ", '"', "math.inc", '"'
    write (outfile1,*) "#include ", '"', "transforms.inc", '"'

    write (outfile1,*) "#declare Camera_1 = camera { "
    write (outfile1,*) "                             angle 15"
    write (outfile1,*) "                             location  <0.0 , ", dwidth(2)/2d0, " ,330.0>"
    write (outfile1,*) "                             right     x*image_width/image_height"
    write (outfile1,*) "                             look_at   <0.0 , ", dwidth(2)/2d0, ", 0.0>"
    write (outfile1,*) "                           }"
    write (outfile1,*) "camera{Camera_1}"

    write (outfile1,*) "light_source{<1000,2500,2500> color White}"
    write (outfile1,*) "light_source{<-1000,-2500,-2500> color White}"


    do n = 1, NOTA
        write (fname1, "(1I0.1)") n
        fname2 = "spl_"//fname1


        write (outfile1,"(1A9,1A14, 1A3)") "#declare ", fname2, " = "
        write (outfile1,*) "  spline {"
        write (outfile1,*) "    natural_spline"

        do time_index = 0, time_span-1
        
        write (outfile1, "(1I9, 1A3, 1e15.6,1A2, 1e15.6,1A2, 1e15.6, 1A2)") time_index, ", <", &
        & r(1,0, n,time_index), ", ",r(2,0, n,time_index), ", ", r(3, 0, n,time_index),">,"
        
        end do
        
        write (outfile1, "(1I9, 1A3, 1e15.6,1A2, 1e15.6,1A2, 1e15.6, 1A2)" ) time_index, ", <", &
        & r(1, 0, n,time_span), ", ",r(2,0, n,time_span), ", ", r(3,0, n,time_span),"> "
        write (outfile1, *) "  }"        
    end do
 
 
    write (outfile1, *)
    write (outfile1, *) "#declare Atom1 ="
    write (outfile1, *) "sphere { <0,0,0>, 0.25 "
    write (outfile1, *) "         texture { pigment{ color rgb<1,0,0>}"
    write (outfile1, *) "                   finish { ambient 0.1 diffuse 0.85  phong 1}"
    write (outfile1, *) "                 }"
    write (outfile1, *) "          scale<2.86,2.86,2.86>"
    write (outfile1, *) "       }"
    
    write (outfile1, *)
    write (outfile1, *) "#declare Atom2 ="
    write (outfile1, *) "sphere { <0,0,0>, 0.25 "
    write (outfile1, *) "         texture { pigment{ color rgb<0,1,0>}"
    write (outfile1, *) "                   finish { ambient 0.1 diffuse 0.85  phong 1}"
    write (outfile1, *) "                 }"
    write (outfile1, *) "          scale<3,3,3>"
    write (outfile1, *) "       }"
        
    do n = 1, NOTA
        write (fname1, "(1I0.1)") n
        fname2 = "spl_"//fname1
        
        if (t(n) == 1) then 
            write (outfile1,*) "object { Atom1"
        else if (t(n) == 2) then
            write (outfile1,*) "object { Atom2"
        end if
        
        write (outfile1,*) "translate ", fname2, "(clock)"
        write (outfile1,*) "}"       
    end do

    
end subroutine

subroutine generate_pair_correlation()
    use common_variable
    integer i, j, t_index, r_index
    integer k1, k2, k3
    integer bin(3, 0:100)
    double precision dr(3), dr_norm, k_width(3)
    integer :: outfile2 
    
    
    
    outfile2 = 103
    
    open (UNIT=outfile2,FILE="out/pair_correlation.csv", STATUS="replace")
    
    do i = 1, 3
    do j = 0, 100
       bin (i, j) = 0
    end do
    end do
    
    
    do i = 1, NOTA
        print *, "Checking ", i
        do j = i+1, NOTA
            if (t(i) == t(j)) then
                t_index = t(i)
            else 
                t_index = 3
            end if
            
            do k1 = -1, 1
            do k2 = -1, 1
            do k3 = -1, 1
            
            k_width(1) = dble(k1)*dwidth(1)
            k_width(2) = dble(k2)*dwidth(2)
            k_width(3) = dble(k3)*dwidth(3)
            
            dr(:) = r(:, 0, i, time_span) - r(:, 0, j, time_span) + k_width(:)
            dr_norm = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
            if (dr_norm <= 10d0) then
                r_index = int(dr_norm*10)
                bin(t_index, r_index) = bin(t_index, r_index) + 2
            end if
            
            end do
            end do
            end do
                        
        end do
    end do
    
    do i = 0, 100
        write (outfile2, *) i, ",",bin(1, i),",", bin(2, i),",", bin(3, i)
    end do
end subroutine