!**********************************************************************
! This revised program aims to utilize the sub_bs_dc subroutine 
! from the fk_openmp module to generate a Green's function bank.
!**********************************************************************
module green_function_utilities
    implicit none
    contains

    subroutine output_gf_bank_to_file(gf_data, gf_filename)
        complex(kind=8), intent(in) :: gf_data(:,:,:,:)
        character(len=*), intent(in) :: gf_filename
        integer :: i, j, k, l

        open(unit=10, file=gf_filename, form='unformatted')
        do i = 1, size(gf_data, 1)
            do j = 1, size(gf_data, 2)
                do k = 1, size(gf_data, 3)
                    do l = 1, size(gf_data, 4)
                        write(10) gf_data(i, j, k, l)
                    end do
                end do
            end do
        end do
        close(10)
    end subroutine output_gf_bank_to_file

end module green_function_utilities

program generate_gf_bank_using_fk_openmp
   use green_function_utilities
   use constants, only : nlay, ndis, nt
   use retrieve_gf, only : block_gg, dt, lnpt, dist_max, dist_min, d_step, dep_max, dep_min, dep_step, get_gf_data, t_cor
   use vel_model_data, only : read_vel_model, update_model 
   use wave_travel, only : trav_fk
   use bessel2, only : load_bessel, deallocate_bessel
   use fk_openmp, only : sub_bs_dc
   implicit none

   integer :: k, iz, nx
   real(kind=8) :: depth
   real(kind=8), allocatable :: distances(:), t0(:)
   real(kind=8) :: tmin(ndis)  ! Assuming ndis is defined and represents the number of distance points
   complex(kind=8), allocatable :: green_function_bank(:,:,:,:)
   character(len=500) :: gf_file, vel_model_path
   character(len=100) :: gf_bank = 'green_function_bank.bin'
   logical :: disp
   
   ! Initialize parameters
   integer, parameter :: nz = 5, nComponents = 8

   call get_command_argument(1, vel_model_path)
   print *, "Velocity Model Path:", trim(vel_model_path)
   if (trim(vel_model_path) == '') then
      print *, "No velocity model path provided. Exiting."
      stop
   endif
   disp = .true.

   ! gf_file = 'Green_strong.txt'
   ! if (disp) gf_file = 'Green_cgps.txt'
   ! call get_gf_data(gf_file, vel_model_path, gf_bank)
   dist_max = 100.0
   dist_min = 0.0
   d_step = 10.0


   ! Calculate nx and distances based on input parameters
   nx = int((dist_max - dist_min) / d_step) + 1
   print *, "Calculated nx: ", nx
   allocate(distances(nx))
   allocate(t0(nx))
   t0 = 0.0

   print *, "Calculating 'dist' and 'ndis' values."
   print *, "dist_min:", dist_min, "dist_max:", dist_max, "d_step:", d_step
   print *, "Calculated nx (number of distances):", nx
   print *, "'ndis' value from module:", ndis
   flush(0)

   call read_vel_model(vel_model_path)

   ! Mimic the distance array calculation as in the working example
   print *, "Distances Array:"
   do k = 1, nx
      distances(k) = dist_min + (k-1) * d_step
      print *, "distances(", k, "): ", distances(k)
   end do

   ! Before loop to generate Green's functions
   print *, "Before entering loop to generate Green's functions."
   
   print *, "Initial Parameters:"
   print *, "dist_max: ", dist_max, " dist_min: ", dist_min, " d_step: ", d_step
   call load_bessel(distances, nx)  ! Make sure load_bessel is modified to handle the dynamic size of `distances`

   ! Generate Green's function bank
   allocate(green_function_bank(nz, nx, nt, nComponents))

   do iz = 1, nz
      depth = dep_min + (iz - 1) * dep_step
      call update_model(depth, 0.0d0)
      call trav_fk(distances, tmin, nx)
      t0 = tmin - t_cor  ! Adjust 't_cor' as needed

      print *, "Before calling sub_bs_dc for depth:", depth
      print *, "nx:", nx, "First distance:", distances(1), "First t0:", t0(1)
      flush(0)

      call sub_bs_dc(nx, distances, t0, green_function_bank(iz, :, :, :), disp, nt, ndis)
      
      print *, "After calling sub_bs_dc for depth:", depth
      flush(0)
   end do

   call output_gf_bank_to_file(green_function_bank, gf_bank)

   call deallocate_bessel()
   deallocate(distances, t0, green_function_bank)

end program generate_gf_bank_using_fk_openmp

