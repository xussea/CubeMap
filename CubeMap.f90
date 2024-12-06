program cube_map

  implicit none
  character(len=100) :: cube_file, file2, file3
  real*8 :: box1, thr, v_err, intg
  real*8 :: start_time, end_time, elapsed_time
  integer :: natoms, nt, n_total
  integer :: i
  real*8, dimension(:,:), allocatable :: grid1, grid2, grid3 
  real*8, dimension(:,:), allocatable :: coord1, coord2, coord3 
  real*8, dimension(:,:), allocatable :: points2, points3  
  real*8, dimension(:,:), allocatable :: dens
  real*8, dimension(3) :: centroid1, centroid2, centroid3 
  real*8, dimension(3,3) :: Rr2, Rr3 

  ! Get the start time 
  call cpu_time(start_time) 

  ! Read cube file name from command line argument
  call getarg(1, cube_file)
  call getarg(2, file2)
  call getarg(3, file3)

  ! Call read xyz file subroutine 
  call load_xyz(file2, coord2, natoms)  
  call load_xyz(file3, coord3, natoms)  
  
  ! Call cubereader subroutine
  call cubereader(cube_file, natoms, box1, nt, coord1, dens, grid1)
  do i = 1, natoms 
        coord1(i, :) = coord1(i, :) / 1.889725886 
  enddo 

  ! Write grid1 array to standard output
  open(unit=22, file='grid1.dat') 
  do i = 1, nt
    write(22,*) grid1(i, :)
  end do
  close(22) 

  ! Calcular centroides 
  call cal_centroide(coord1, natoms, centroid1)
  call cal_centroide(coord2, natoms, centroid2)
  call cal_centroide(coord3, natoms, centroid3)

  ! Calcular Matriz de Rotacion 
  call kabsch_rmsd(coord1, coord2, centroid1, centroid2, natoms, Rr2)
  call kabsch_rmsd(coord1, coord3, centroid1, centroid3, natoms, Rr3)

  ! Rotate the grid2 and grid3 using grid1 
  call rotate_grid(grid1, dens, centroid1, centroid2, Rr2, nt, grid2) 
  open(unit=33, file='grid2.dat') 
  do i = 1, nt
    write(33,*) grid2(i, :)
  end do
  close(33) 
       
  call rotate_grid(grid1, dens, centroid1, centroid3, Rr3, nt, grid3) 
  open(unit=44, file='grid3.dat') 
  do i = 1, nt
    write(44,*) grid3(i, :)
  end do
  close(44) 

  ! Definir threshold 
  thr = box1 / 2.0 

  ! Detectar los puntos que se solapan 
  call find_overlapping_points(grid2, grid3, thr, nt, points2, points3, n_total, v_err) 

  ! Guardar los puntos solapados en el un archivo 
  open(unit=55, file='points2.dat') 
  do i = 1, n_total 
          write(55,*) points2(i,:)
  enddo
  close(55) 
  
  open(unit=66, file='points3.dat') 
  do i = 1, n_total 
          write(66,*) points3(i,:)
  enddo
  close(66) 
  
  ! Detectar los puntos que se solapan
  call cal_integral(points2, points3, box1, intg)  
  
  !Guardar el error y el tiempo de calculo 
  open(unit=77, file='results.dat')
  write(77,*) "Error: ", v_err
  call cpu_time(end_time) 
  elapsed_time = end_time - start_time 
  write(77,*) "Elapsed time: ", elapsed_time
  write(77,*) "Integral: ", intg 
  close(77) 

contains 

        subroutine load_xyz(file_name, coords, natoms) 

                implicit none 
                character(len=*), intent(in) :: file_name 
                real*8, dimension(:,:), allocatable, intent(out) :: coords
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

        subroutine cubereader(cube_file, natoms, box1, nt, coord, dens_column, grid1)

                 implicit none
                 character(len=100), intent(in) :: cube_file
                 integer, intent(in) :: natoms
                 real*8 :: a, b, c
                 real*8, intent(out) :: box1
                 real*8, dimension(:,:), allocatable :: grid1
                 real*8, dimension(:,:), allocatable, intent(out) :: coord, dens_column
                 real*8, dimension(:), allocatable :: dens
                 integer :: i, j, k, t, dimensiones
                 integer :: Nx, Ny, Nz
                 integer, intent(out) :: nt
                 integer, dimension(:), allocatable :: z
                 real*8, dimension(3) :: origin, box
                 real*8, dimension(:,:), allocatable :: grid

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

        subroutine cal_centroide(coord, natoms, centroide)
        
                implicit none
                integer*8 :: i
                integer, intent(in) :: natoms
                real*8, dimension(:,:), intent(in) :: coord
                real*8, dimension(3), intent(out) :: centroide
                
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
                integer :: info, lwork, d
                double precision, dimension(:), allocatable :: work
                integer, dimension(:), allocatable :: iwork
                double precision :: det_U, det_Vt 
    
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
            
                ! Realizar la descomposición SVD
                call dgesdd('A', 3, 3, C, 3, S, U, 3, Vt, 3, work, lwork, iwork, info)
    
                ! Verifica si se necesita corrección de reflexión
                call deter(U, det_U) 
                call deter(Vt, det_Vt)
    
                d = sign(1.00d0,det_U) * sign(1.00d0,det_Vt)
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
                real*8, dimension (:,:), intent(in) :: grid1, dens 
                real*8, dimension (:,:), intent(in) :: Rr 
                real*8, dimension (:), intent(in) :: centroid1, centroid2
                real*8, dimension (3) :: cent1, cent2
                real*8, dimension (:,:), allocatable :: grid, grid_rot
                real*8, dimension (:,:), allocatable, intent(out) :: grid2 

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

        subroutine find_overlapping_points(grid1, grid2, thr, nt, points1, points2, n_total, v_err) 

                implicit none
                real*8, dimension(:,:), intent(in) :: grid1, grid2 
                real*8, intent(in) :: thr 
                integer, intent(in) :: nt
                real*8, dimension(:,:), allocatable, intent(out) :: points1, points2 
                integer, intent(out) :: n_total
                real*8, intent(out) :: v_err
                integer :: i, j, n_k, n_v 
                real*8 :: n_nodup 
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
               
                allocate(k_set(nt)) 
                allocate(v_set(nt)) 
               

                k_set = .false.
                v_set = .false.  

                do i = 1, nt 
                    do j = 1, nt
                         if (abs(grid1(i,1) - grid2(j,1)) <= thr .and. abs(grid1(i,2) - grid2(j,2)) <= thr .and. &
                            abs(grid1(i,3) - grid2(j,3)) <= thr) then 
                            n_total = n_total + 1
                            k_set(i) = .true.
                            v_set(j) = .true.
                            points1(n_total, :) = grid1(i,:)  
                            points2(n_total, :) = grid2(j,:)
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
               real*8, dimension(:,:), intent(in) :: points1, points2
               real*8, intent(in) :: box1 
               real*8, intent(out) :: intg
               real*8 :: vol
               integer :: i 

               vol = box1 * box1 * box1 
               intg = 0.0 

               do i = 1, size(points1, 1) 
                    intg = intg + points1(i,4) * points2(i,4) * vol 
               enddo    
               return

        end subroutine cal_integral

end program cube_map
