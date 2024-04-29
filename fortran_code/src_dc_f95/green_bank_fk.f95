!       Make the Green Function bank 
!       Input dist_max dist_min d_step
!       Input dep_max dep_min dep_step
!       Input Lnpt dt
!       Input Green_name
!       output Green
program green_bank_fk


   use constants, only : nlay, ndis, nt
   use retrieve_gf, only : block_gg, dt, lnpt, dist_max, dist_min, d_step, dep_max, dep_min, dep_step, get_gf_data, t_cor
   use vel_model_data, only : read_vel_model, update_model, src, rcv 
   use wave_travel, only : trav_fk
   use bessel2, only : load_bessel
   use fk, only : sub_bs_dc
   implicit none
   real(kind=8) dist(ndis)
   real(kind=8) t0(ndis)
   real(kind=8) tmin(ndis)
   real(kind=8) depth
   real(kind=8) :: zero_real8 = 0.0d0
   real :: green(nt, 8, ndis)
   integer iz, k, ll, n_com, nd_max, npt, ntc, nx, nz, current_ntc, current_n_com, j
   logical :: disp

   character(len=100) gf_file, vel_model, gf_bank, input
!
   call getarg(1, input)
   disp = (input.eq.'cgps')

   gf_file = 'Green.in'
   if (disp) gf_file = 'Green_cgps.in'
   call get_gf_data(gf_file, vel_model, gf_bank)
   npt=2**lnpt
!
   if(abs(dist_min-int(dist_min/d_step+0.001)*d_step).gt.1.0e-4)then
      write(*,*) "To improve the speed of this calculation, "
      write(*,*) "The minimum distance dist_min = m *d_step "
      write(*,*)dist_min,d_step,abs(dist_min-int(dist_min/d_step+0.001)*d_step)
!    pause
   endif
   if(abs(dist_max-int(dist_max/d_step+0.001)*d_step).gt.1.0e-4)then
      write(*,*) "To improve the speed of this calculation, "
      write(*,*) "The maximum distance dist_min = m *d_step "
      write(*,*)dist_max,d_step,abs(dist_max-int(dist_max/d_step+0.001)*d_step)
!	pause
   endif
   nd_max=int(dist_max/d_step+0.1)+1

   if(npt.gt.nt)then
      write(0,*)"The length of green's function must be shorten"
      stop
   endif
   if(dep_max.lt.dep_min)then
      write(*,*)"The depth region is wrong"
      stop
   endif
   if(dist_max.lt.dist_min)then
      write(*,*)"distance region is wrong"
      stop
   endif

   nx=int((dist_max-dist_min)/d_step)+1
   nz=int((dep_max-dep_min)/dep_step)+1

   write(*,*)'total =',nx,nz
!	open(100,file='stored_bessel')
!	open(101,file='computed_bessel')
   if(nx.gt.ndis) then
      write(*,*)"please reduce the number of Green function"
      stop
   endif
   !open(11,file=gf_bank,status='unknown',access='direct',recl=block_gg) ! Crear archivo binar
   open(11, file=gf_bank, form='unformatted', status='replace', access='direct', recl=10000) ! Crear archivo binario

   call read_vel_model(vel_model)

   do k=1,nx
      dist(k)=dist_min+(k-1)*d_step
   enddo
   write(*,*) '<green_bank_fk.f95>: Distance grid initialization:', dist
   write(*,*) '<green_bank_fk.f95>: Depth grid initialization:', dep_min, dep_step, nz

   call load_bessel(dist, nd_max)
   ll=0

   ! Before the loop, print the dimensions for debugging
   write(*,*) 'Writing dimensions to file:', nx, nz

   open(unit=11, file='green_function_bank.bin', form='unformatted', status='replace')

   do iz=1,nz,5  ! Skipping some depth indices for faster execution
   !do iz=1,nz  ! Skipping some depth indices for faster execution
      depth=dep_min + (iz-1)*dep_step
      call update_model(depth, zero_real8)
      call trav_fk(dist, tmin, nx)
      do k=1,nx,5  ! Skipping some distance indices for faster execution
      !do k=1,nx  ! Skipping some distance indices for faster execution
         write(*,*) 'Writing record for depth index', iz, 'and distance index', k
         write(*,*) 'Depth:', depth, 'Distance:', dist(k), 't0:', t0(k)
         write(*,*) 'tmin values before setting t0:', tmin(k)
         write(*,*) 't_cor value used:', t_cor
         ! t0(k) = tmin(k) - t_cor
         t0(k) = tmin(k)
         if (t0(k) < 0.0d0) then
            t0(k) = 0.0d0
            write(*, '(A, I3, A, F10.5)') 'Adjusted negative t0 at index ', k, ' to zero.'
         endif
         write(*, '(A, I3, A, F10.5, A)', advance='no') &
         'Debug at green_bank_fk.f95 - t0 set: Depth Index ', iz, ' Value: ', t0(k), ' (after correction)'
         print *, 'Updated t0 at ix:', k, 'with t0:', t0(k)
         write(*, '(A, I3, A, F10.5, A)', advance='no') &
         'Debug at green_bank_fk.f95 - t0 set: Depth Index ', iz, ' Value: ', t0(k), ' (after correction)'
         call sub_bs_dc(nx, dist, t0, green, disp)  ! Pass updated t0 values to sub_bs_dc
         ll = (iz - 1) * nx + k
         write(*, '(A, I3, A, I3, A, F10.5, A)', advance='no') &
         'Debug at green_bank_fk.f95 - Writing t0 to file: Depth Index ', iz, ' Distance Index ', k, ' t0 Value: ', t0(k)
         ! Writing only essential components to reduce file size and increase speed
         write(*, '(A, I3, A, F10.5, A, F10.5)') &
         'Depth Index:', iz, ' Distance:', dist(k), ' t0 Updated:', t0(k)

         write(*, '(A, I3)') 'Intensities at depth index ', iz
         do j = 1, nt  ! Assuming 'nt' is the length of the first dimension of 'green'
            write(*, '(10E10.3)') green(j, 1, k)
         end do

         write(11) iz, k, dist(k), t0(k), depth, dt, green(:, :, k)
      enddo
    enddo
   close(11)


end program green_bank_fk

