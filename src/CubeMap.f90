program cube_map

  implicit none
  integer :: nt, n_total
  double precision :: v_err, intg
  double precision :: start_time, end_time, elapsed_time
  double precision, dimension(3) :: thr 
  double precision, dimension(:,:), allocatable :: grid2, grid3 
  double precision, dimension(:,:), allocatable :: points2, points3  
  logical :: c, multic, ex
  character(len=300) :: cube, cube1, cube2
  character(len=300) :: refm1, refm2, boxfile
  character(len=300) :: mon1, mon2, exclu
  ! Get the start time 
  call cpu_time(start_time) 

  !Flags
  call call_flags(c, multic, ex, &
                cube, cube1, cube2, refm1, refm2, boxfile, &
                mon1, mon2, exclu)
  
  if (c) then
    call simple_cube(cube,mon1,mon2,ex,exclu,grid2,grid3,nt,thr)
  else if (multic) then
    call multi_cube(cube1,cube2,refm1,refm2,boxfile,mon1,mon2,ex,exclu,grid2,grid3,nt,thr)
  end if

  ! Detectar los puntos que se solapan 
  call find_overlapping_points(grid2, grid3, thr, nt, points2, points3, n_total, v_err, multic) 

  call write_grid(points2,n_total,'points2.dat')
  call write_grid(points3,n_total,'points3.dat')
                                                                                         
  ! Detectar los puntos que se solapan
  call cal_integral(points2, points3, thr, intg)  

  !Guardar el error y el tiempo de calculo 
  open(unit=77, file='results.dat')
  write(77,*) "Error: ", v_err
  call cpu_time(end_time) 
  elapsed_time = end_time - start_time 
  write(77,*) "Elapsed time: ", elapsed_time
  write(77,*) "Integral: ", intg 
  close(77) 

contains 
        subroutine simple_cube(cube,mon1,mon2,ex,exclu,grid2,grid3,nt,thr)
                implicit none
                integer :: natoms, nexclude, i
                integer, intent(out) :: nt
                double precision :: box1
                double precision, dimension(3) :: centroid1, centroid2, centroid3
                double precision, dimension(3), intent(out):: thr
                double precision, dimension(3,3) :: Rr2, Rr3 
                double precision, dimension(:,:), allocatable :: grid1, coord1, coord2, coord3 
                double precision, dimension(:,:), allocatable :: dens
                double precision, dimension(:,:), allocatable, intent(out) :: grid2, grid3 
                logical, intent(in) :: ex
                character(len=300), intent(in) :: cube, mon1, mon2, exclu

                call load_xyz(mon1, coord2, natoms)  
                call load_xyz(mon2, coord3, natoms)  
                call cubereader(cube, natoms, box1, nt, coord1, dens, grid1)
                
                do i = 1, natoms 
                  coord1(i, :) = coord1(i, :) / 1.889725886 
                enddo 
                call write_grid(grid1,nt,'grid1.dat')
                
                if (ex) then
                  call xyzrestrain(coord1, natoms, exclu, nexclude)
                  call xyzrestrain(coord2, natoms, exclu, nexclude)
                  call xyzrestrain(coord3, natoms, exclu, nexclude)
                  natoms=natoms-nexclude
                end if

                call cal_centroide(coord1, natoms, centroid1)
                call cal_centroide(coord2, natoms, centroid2)
                call cal_centroide(coord3, natoms, centroid3)

                call kabsch_rmsd(coord1, coord2, centroid1, centroid2, natoms, Rr2)
                call kabsch_rmsd(coord1, coord3, centroid1, centroid3, natoms, Rr3)

                thr = box1

                call rotate_grid(grid1, dens, centroid1, centroid2, Rr2, nt, grid2) 
                call write_grid(grid2,nt,'grid2.dat')
                call rotate_grid(grid1, dens, centroid1, centroid3, Rr3, nt, grid3)
                call write_grid(grid3,nt,'grid3.dat')

        end subroutine simple_cube

        subroutine multi_cube(cube1,cube2,refm1,refm2,boxfile,mon1,mon2,ex,exclu,grid2,grid3,nt,thr)
                implicit none
                integer :: natoms, nexclude, i
                integer, intent(out):: nt
                double precision, dimension(3) :: centroidref1, centroidref2, centroid2, centroid3 
                double precision, dimension(3), intent(out):: thr
                double precision, dimension(3,3) :: Rr2, Rr3 
                double precision, dimension(:,:), allocatable :: gridref1, gridref2 
                double precision, dimension(:,:), allocatable :: coordref1, coordref2, coord2, coord3 
                double precision, dimension(:,:), allocatable :: densref1, densref2
                double precision, dimension(:,:), allocatable, intent(out) :: grid2, grid3 
                logical, intent(in) :: ex
                character(len=300), intent(in) :: cube1, cube2, refm1, refm2, boxfile, mon1, mon2, exclu

                call load_xyz(mon1, coord2, natoms)  
                call load_xyz(mon2, coord3, natoms)  
                call load_xyz(refm1, coordref1, natoms)  
                call load_xyz(refm2, coordref2, natoms)  


                if (ex) then
                  call xyzrestrain(coord2, natoms, exclu, nexclude)
                  call xyzrestrain(coord3, natoms, exclu, nexclude)
                  call xyzrestrain(coordref1, natoms, exclu, nexclude)
                  call xyzrestrain(coordref2, natoms, exclu, nexclude)
                  natoms=natoms-nexclude
                end if

                call cal_centroide(coordref1, natoms, centroidref1)
                call cal_centroide(coordref2, natoms, centroidref2)
                call cal_centroide(coord2, natoms, centroid2)
                call cal_centroide(coord3, natoms, centroid3)

                call kabsch_rmsd(coordref1, coord2, centroidref1, centroid2, natoms, Rr2)
                call kabsch_rmsd(coordref2, coord3, centroidref2, centroid3, natoms, Rr3)

                open(unit=44,file=boxfile)
                read(44,*) thr(1) , thr(2), thr(3), nt
                close(44)

                thr=thr*1.889725886
                                                      
                allocate(densref1(nt,1))
                allocate(densref2(nt,1))
                                                      
                allocate(gridref1(nt,4))
                allocate(gridref2(nt,4))
                
                open(unit=44,file=cube1)
                do i=1,nt
                   read(44,*)gridref1(i,:)
                   densref1(i,1)=gridref1(i,4)!*1.889725886
                end do
                close(44)
                                               
                open(unit=44,file=cube2)
                do i=1,nt
                   read(44,*)gridref2(i,:)
                   densref2(i,1)=gridref2(i,4)!*1.889725886
                end do
                close(44)

                !Los grids aqui están en anstrong también

                call rotate_grid(gridref1*1.889725886 , densref1, centroidref1, centroid2, Rr2, nt, grid2)
                call write_grid(grid2/1.889725886,nt,'grid2.dat')
                call rotate_grid(gridref2*1.889725886 , densref2, centroidref2, centroid3, Rr3, nt, grid3)
                call write_grid(grid3/1.889725886,nt,'grid3.dat')

        end subroutine multi_cube

        subroutine call_flags(c, multic, ex, &
                      cube, cube1, cube2, refm1, refm2, boxfile, &
                      mon1, mon2, exclu)
          implicit none
          logical, intent(out) :: c, multic, ex
          logical :: mon
          character(len=300), intent(out) :: cube, cube1, cube2
          character(len=300), intent(out) :: refm1, refm2, boxfile
          character(len=300), intent(out) :: mon1, mon2, exclu
          integer :: i, narg
          character(len=300) :: arg
          
          c      = .false.
          multic = .false.
          mon    = .false.
          ex     = .false.
          
          cube    = ''
          cube1   = ''
          cube2   = ''
          refm1   = ''
          refm2   = ''
          boxfile = ''
          mon1    = ''
          mon2    = ''
          exclu   = ''
          
          narg = command_argument_count()
          i = 1
          
          do while (i <= narg)
              call get_command_argument(i, arg)
              select case (trim(arg))
              case ('--c','--cube')
                  if (multic) then
                      write(*,*) 'ERROR: Use either --cube or --multicube, not both'
                      stop
                  end if
                  c = .true.
                  if (i+1 > narg) then
                      write(*,*) 'ERROR: Missing cube file after --cube'
                      stop
                  end if
                  i = i + 1
                  call get_command_argument(i, cube)
              case ('--mc','--multicube')
                  if (c) then
                      write(*,*) 'ERROR: Use either --cube or --multicube, not both'
                      stop
                  end if
                  multic = .true.
                  if (i+5 > narg) then
                      write(*,*) 'ERROR: --multicube needs 5 arguments'
                      stop
                  end if
                  i = i + 1
                  call get_command_argument(i, cube1)
                  i = i + 1
                  call get_command_argument(i, cube2)
                  i = i + 1
                  call get_command_argument(i, refm1)
                  i = i + 1
                  call get_command_argument(i, refm2)
                  i = i + 1
                  call get_command_argument(i, boxfile)
              case ('--m','--monomers')
                  mon = .true.
                  if (i+2 > narg) then
                      write(*,*) 'ERROR: --monomers needs 2 arguments'
                      stop
                  end if
                  i = i + 1
                  call get_command_argument(i, mon1)
                  i = i + 1
                  call get_command_argument(i, mon2)
              case ('--ex','--exclude')
                  ex = .true.
                  if (i+1 > narg) then
                      write(*,*) 'ERROR: Missing file after --exclude'
                      stop
                  end if
                  i = i + 1
                  call get_command_argument(i, exclu)
              case ('--h','--help')
                  write(*,*) 'Usage:'
                  write(*,*) 'cubemap --c ref.cube --m mon1.xyz mon2.xyz [--ex Excludefile]'
                  write(*,*) 'cubemap --mc ref1.cube ref2.cube ref1.xyz ref2.xyz box.dat --m mon1.xyz mon2.xyz [--ex Excludefile]'
                  write(*,*) 'OPTIONS:'
                  write(*,*) '--ex or --exclude, To exclude selected atoms in the CubeMap approach'
                  write(*,*) '--m or --monomers, To provide the selected monomers'
                  write(*,*) '--mc or --multicube, The program use more than one reference provided cube'
                  write(*,*) '--c or --cube, The program use only the cube provided'
                  write(*,*)
                  stop
              case default
                  write(*,*) 'ERROR: Unknown option "', trim(arg), '", use cubemap --h or --help for more information'
                  stop
              end select
              i = i + 1
          end do
          
          if (.not. c .and. .not. multic) then
              write(*,*) 'ERROR: You must specify --cube or --multicube'
              stop
          end if
          
        end subroutine call_flags


        subroutine load_xyz(file_name, coords, natoms) 

                implicit none 
                character(len=*), intent(in) :: file_name 
                double precision, dimension(:,:), allocatable, intent(out) :: coords
                integer :: i, natoms
                character(len=2) :: atom_type

                ! Abrir archivo 
                open(unit=20, file=file_name, status='old') 
                
                read(20,*) natoms 

                !Almacenar las coordenadas 
                allocate(coords(natoms,3)) 
                
                do i = 1, natoms 
                        read(20,*) atom_type, coords(i,1), coords(i,2), coords(i,3) 
                enddo 

                close(20) 
        
        end subroutine load_xyz 

                subroutine xyzrestrain(coords,NA,filename,NEx)

                implicit none
                character(len=*), intent(in):: filename
                double precision, dimension(:,:), allocatable :: coordsnew
                double precision, dimension(:,:), allocatable, intent(inout) :: coords
                integer, dimension(:), allocatable :: exclude
                integer:: i, j, k
                integer, intent(inout):: NA,NEx
                logical:: outex
                
                open (unit=30, file=filename, status='old')
                read(30,*) NEx
                allocate(exclude(NEx))
                allocate(coordsnew(NA-NEx,3))
                read(30,*) exclude(1:NEx)
                close(30)


                j = 0
                do i = 1, NA
                        outex = .false.
                        do k = 1, NEx
                             if (i == exclude(k)) then
                                outex = .true.
                                exit
                             end if
                        end do

                        if (.not. outex) then
                                j = j + 1
                                coordsnew(j,:) = coords(i,:)
                        end if
                end do

                deallocate(exclude)
                call move_alloc(coordsnew,coords)
        end subroutine xyzrestrain

        subroutine cubereader(cube_file, natoms, box1, nt, coord, dens_column, grid1)

                 implicit none
                 character(len=100), intent(in) :: cube_file
                 integer, intent(in) :: natoms
                 double precision :: a, b, c
                 double precision, intent(out) :: box1
                 double precision, dimension(:,:), allocatable :: grid1
                 double precision, dimension(:,:), allocatable, intent(out) :: coord, dens_column
                 double precision, dimension(:), allocatable :: dens
                 integer :: i, j, k, t
                 integer :: Nx, Ny, Nz
                 integer, intent(out) :: nt
                 integer, dimension(:), allocatable :: z
                 double precision, dimension(3) :: origin, box
                 double precision, dimension(:,:), allocatable :: grid

                 ! Open cube file
                 open(unit=10, file=cube_file)

                 ! Read header information
                 read(10,*)
                 read(10,*)
                 read(10,*) c, origin(1:3)
                 read(10,*) Nx, box(1)
                 read(10,*) Ny, a, box(2)
                 read(10,*) Nz, a, b, box(3)

                 ! Calculate total number of grid points
                 nt = Nx * Ny * Nz

                 ! Allocate memory for arrays
                 allocate(z(natoms))
                 allocate(coord(natoms,3))
                 allocate(grid(nt, 3))
                 allocate(dens(nt))
                 allocate(grid1(nt, 4))
                 allocate(dens_column(nt,1)) 

                 box1 = box(1)

                 ! Read atom coordinates
                 do i = 1, natoms
                   read(10,*) z(i), a, coord(i, 1), coord(i, 2), coord(i, 3)
                 end do

                 ! Skip one line
                 read(10,*)

                 ! Read density values
                 t = 0
                 do i = 1, Nx
                   do j = 1, Ny
                     t = t + 1
                     Read(10,'(6E13.5)') dens((t-1)*Nz+1:(t-1)*Nz+Nz)
                   end do
                 end do

                 ! Populate grid array
                 t = 1
                 do i = 1, Nx
                   do j = 1, Ny
                     do k = 1, Nz
                       grid(t,1) = origin(1) + (i-1) * box(1)
                       grid(t,2) = origin(2) + (j-1) * box(2)
                       grid(t,3) = origin(3) + (k-1) * box(3)
                       t = t + 1
                     end do
                   end do
                 end do

                 ! Copy data to grid1 array
                 do i = 1, nt
                   grid1(i, 1:3) = grid(i, :)
                   grid1(i, 4) = dens(i)
                 end do

                 ! Convertir dens a una columna
                 dens_column(:,1) = dens

                 ! Clean up
                 close(10)
                 deallocate(z, grid, dens)

        end subroutine cubereader


        subroutine write_grid(grid,nt,name)
                implicit none
                integer :: i, nt
                double precision, dimension(:,:), intent(in) :: grid
                character(len=*)::name

                open(unit=66, file=name)
                do i =1 , nt
                  write(66,*) grid(i,:)
                end do
                close(66)
        end subroutine write_grid

        subroutine cal_centroide(coord, natoms, centroide)
        
                implicit none
                integer :: i
                integer, intent(in) :: natoms
                double precision, dimension(:,:), intent(in) :: coord
                double precision, dimension(3), intent(out) :: centroide
                
                centroide(1:3) = 0
                
                do i = 1, natoms
                        centroide(1:3) = centroide(1:3) + coord(i,1:3)
                end do
                
                centroide(1:3) = centroide(1:3) / natoms
                return
        
        end subroutine cal_centroide

        subroutine kabsch_rmsd(coord1, coord2, centroid1, centroid2, natoms, Rr)
            
                implicit none
                integer :: i, j 
                double precision, dimension(:,:), intent(in)  :: coord1, coord2
                double precision, dimension(:), intent(in) :: centroid1, centroid2
                double precision, dimension(3,3), intent(out) :: Rr
                integer, intent(in) :: natoms
                double precision, dimension(natoms,3) :: P, Q
                double precision :: C(3,3), S(min(3,3)), U(3,3), Vt(3,3)
                integer :: info, lwork
                double precision, dimension(:), allocatable :: work
                integer, dimension(:), allocatable :: iwork
                double precision :: det_U, det_Vt ,d
    
                ! Calcular las matrices P y Q
                do i = 1, natoms
                    do j = 1, 3
                        P(i,j) = coord1(i,j) - centroid1(j)
                        Q(i,j) = coord2(i,j) - centroid2(j)
                    end do
                end do
            
                ! Calcular la matriz de covarianza C
                C = matmul(transpose(P), Q)
            
                ! Espacio de trabajo para la descomposición SVD
                lwork = -1
                allocate(work(1))
                allocate(iwork(8*min(3,3)))
                call dgesdd('A', 3, 3, C, 3, S, U, 3, Vt, 3, work, lwork, iwork, info)
                lwork = max(1, int(work(1)))

                deallocate(work)
                allocate(work(lwork))
            
                ! Realizar la descomposición SVD
                call dgesdd('A', 3, 3, C, 3, S, U, 3, Vt, 3, work, lwork, iwork, info)
    
                ! Verifica si se necesita corrección de reflexión
                call deter(U, det_U) 
                call deter(Vt, det_Vt)
    
                !d = sign(1.00d0,det_U) * sign(1.00d0,det_Vt)
                d = (det_U)/abs(det_U) * det_Vt/abs(det_Vt)
                if (d < 0.0) then
                    ! Cambia el signo de la última columna de V y W
                    U(:, 3) = -U(:, 3)
                    !S(:, 3) = -S(:, 3)
                endif 
    
                ! Crear la matriz de rotación U
                Rr = matmul(U, Vt)
                
                deallocate(work)
                deallocate(iwork)
                
                return
        end subroutine kabsch_rmsd
        
        subroutine deter(A, det_A) 
            
            implicit none 
            double precision, dimension(3,3), intent(in)  :: A 
            double precision, intent(out) :: det_A 
            double precision :: T1, T2, T3 

            T1  = A(2,2)*A(3,3) - A(3,2)*A(2,3)
            T2  =-A(2,1)*A(3,3) + A(3,1)*A(2,3)
            T3  = A(2,1)*A(3,2) - A(3,1)*A(2,2)
            det_A  = A(1,1)*T1 + A(1,2)*T2 + A(1,3)*T3

            return 
        end subroutine deter 

        subroutine rotate_grid(grid1, dens, centroid1, centroid2, Rr, nt, grid2) 

                implicit none 
                integer :: i 
                integer, intent(in) :: nt 
                double precision, dimension (:,:), intent(in) :: grid1, dens 
                double precision, dimension (:,:), intent(in) :: Rr 
                double precision, dimension (:), intent(in) :: centroid1, centroid2
                double precision, dimension (3) :: cent1, cent2
                double precision, dimension (:,:), allocatable :: grid, grid_rot
                double precision, dimension (:,:), allocatable, intent(out) :: grid2 
                
                allocate(grid2(nt,4))
                allocate(grid(nt,4))
                allocate(grid_rot(nt,4))
                
                !Move the grid to the mass center of molecule 1
                cent1 = centroid1 * 1.8897259886 
                do i = 1, nt 
                        grid(i,1:3) = grid1(i,1:3) - cent1(1:3)
                enddo         

                ! Rotate the grid
                grid_rot = matmul(grid(:,1:3), Rr) 
                
                ! Add the mass center of the molecule 2 to the grid 
                cent2 = centroid2 * 1.8897259886 
                do i = 1, nt 
                        grid_rot(i,1:3) = grid_rot(i,1:3) + cent2(1:3)
                enddo         

                 ! Copy data to grid1 array
                 do i = 1, nt
                   grid2(i, 1:3) = grid_rot(i, :)
                   grid2(i, 4) = dens(i,1)
                 end do
                
                 deallocate(grid_rot, grid)  
                return 
        
        end subroutine rotate_grid 

        subroutine find_overlapping_points(grid1, grid2, thr, nt, points1, points2, n_total, v_err, mc) 

                implicit none
                double precision, dimension(:,:), intent(in) :: grid1, grid2 
                double precision, dimension(3), intent(in) :: thr
                integer, intent(in) :: nt
                double precision, dimension(:,:), allocatable, intent(out) :: points1, points2 
                integer, intent(out) :: n_total
                double precision, intent(out) :: v_err
                integer :: i, j, n_k, n_v 
                double precision :: n_nodup 
                logical, intent(in):: mc
                logical, dimension(:), allocatable :: k_set, v_set 

                
                ! k_set  and v_set son matrices logicas, es decir, son matrices que utilizaremos para realizar un seguimiento de si
                ! un punto del grid1 y grid2 ha sido encontrado como superpuesto durante el proceso de búsqueda. Incialmente, todos
                ! los elementos de estas matrices se establecen en '.false.' indicando que ningún punto ha sido encontrado.
                ! Estas matrices se utilizan en el bucle doble dentro de la subroutine "find_overlapping_points" para realizar un
                ! seguimineto de los puntos que ya se han encontrado como superpuestos y eviat contarlos múltiples veces. Cuando se
                ! encuentra un punto superpuesto entre 'grid1' y 'grid2', el valor correspondiuente en 'v_set' y 'k_set' se
                ! establece en '.true.' para indicar que ese punto ya se ha contado como superpuesto. Esto garantiza que un punto no
                ! se cuente como superpuesto más de una vez en el cálculo del número total de puntos superpuestos. 

                n_total = 0

                ! Allocate memory for points1, points2, k_set and v_set arrays 
                allocate(points1(nt, 4)) 
                allocate(points2(nt, 4)) 
                points1=0.0
                points2=0.0
               
                allocate(k_set(nt)) 
                allocate(v_set(nt)) 

                k_set = .false.
                v_set = .false.  

                do i = 1, nt 
                    do j = 1, nt
                         if (abs(grid1(i,1) - grid2(j,1)) .le. thr(1)/2.0 .and. abs(grid1(i,2) - grid2(j,2)) .le. thr(2)/2.0 .and. &
                            abs(grid1(i,3) - grid2(j,3)) .le. thr(3)/2.0) then 
                            n_total = n_total + 1
                            k_set(i) = .true.
                            v_set(j) = .true.
                            points1(n_total, :) = grid1(i,:)  
                            points2(n_total, :) = grid2(j,:)
                            if (mc) then
                                exit
                            endif
                        endif
                    enddo 
                enddo 

                ! Calcular error  
                n_k = count(k_set)
                n_v = count(v_set) 
                n_nodup = (n_total - n_k) + (n_total - n_v) 
                v_err = real(n_nodup) / real(n_total) 
                 
                deallocate(k_set) 
                deallocate(v_set) 
                return

        end subroutine find_overlapping_points

        subroutine cal_integral(points1, points2, box1, intg) 

               implicit none 
               double precision, dimension(:,:), intent(in) :: points1, points2
               double precision, dimension(3), intent(in) :: box1 
               double precision, intent(out) :: intg
               double precision :: vol
               integer :: i 

               vol = box1(1)*box1(2)*box1(3)
               intg = 0.0 
               do i = 1, size(points1, 1) 
                    intg = intg + points1(i,4) * points2(i,4) * vol 
               enddo    
               return

        end subroutine cal_integral

end program cube_map
