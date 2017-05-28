mpif90 -O3 -c MD_MPI_DLL.f90 mpif90
mpif90 -O3 -o MD_MPI MD.MPI.f90 MD_MPI_DLL.o
mpif90 -fast -Mr8 -Munroll=n:10 -Minline=levels:30 -Mvect -Minfo=loop -o MD_Analyze MD.Analysis.f90
