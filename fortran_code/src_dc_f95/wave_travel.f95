!      program travel
!============================================================
! Calculate travel time for horizontal layered model.
! It outputs both times for first arrival and direct arrival.
! Revision history:
!	Nov. 17, 2008	Now it can take receiver at any depth.
!			Lupei Zhu
!		
!	Feb. 20, 2012  Change it into subroutine, Chen Ji
!============================================================
!
!     be careful: src_lay and rcv_lay might be changed after calling
!     this subroutine
!
module wave_travel

    use constants, only: nlay
    use vel_model_data, only: n_layers_new, rcv, src, new_thick, slowness
    implicit none

contains

subroutine trav_fk(dis, tmin, nx)
   implicit none
   integer, intent(in) :: nx
   real(kind=8), intent(in) :: dis(nx)
   real(kind=8), intent(out) :: tmin(nx)
   real(kind=8) :: ray_len(2,nlay), t, td, x, aa, t0, p0
   complex(kind=8) :: pd, p, pc
   integer :: topp, bttm, i, k, src_lay, rcv_lay

   src_lay = src - 1
   rcv_lay = rcv - 1
   if (src_lay < rcv_lay) then
      i = src_lay
      src_lay = rcv_lay
      rcv_lay = i
   endif
   if (n_layers_new > nlay) then
      write(0,*) 'Number of layers exceeds max. allowed'
      return
   endif

   do i = 1, nx
      x = dis(i)
      topp = rcv_lay + 1
      bttm = src_lay
      aa = 1.d+20
      do k = topp, bttm
         ray_len(1,k) = new_thick(k)
         ray_len(2,k) = 0.d0
         aa = min(aa, slowness(k))
      enddo
      pd = cmplx(sqrt(aa), 1.d-20, kind=8)
      call findp0(x, pd, ray_len, topp, bttm)
      td = taup(pd, x, ray_len, topp, bttm)
      t0 = td
      p0 = pd

      do bttm = src_lay + 1, n_layers_new - 1
         ray_len(1, bttm) = 2.d0 * new_thick(bttm)
         ray_len(2, bttm) = 0.d0
         aa = min(aa, slowness(bttm))
         p = cmplx(sqrt(aa), 1.d-20, kind=8)
         call findp0(x, p, ray_len, topp, bttm)
         aa = min(aa, slowness(bttm + 1))
         pc = cmplx(sqrt(aa), 1.d-20, kind=8)
         if (real(p) > real(pc)) p = pc
         t = taup(p, x, ray_len, topp, bttm)
         if (t < t0) then
            t0 = t
            p0 = p
         endif
      enddo

      bttm = src_lay
      do topp = rcv_lay, 1, -1
         if (topp == 0) exit
         ray_len(1, topp) = 2.d0 * new_thick(topp)
         ray_len(2, topp) = 0.d0
         aa = min(aa, slowness(topp))
         p = cmplx(sqrt(aa), 1.d-20, kind=8)
         call findp0(x, p, ray_len, topp, bttm)
         if (topp > 1) aa = min(aa, slowness(topp - 1))
         pc = cmplx(sqrt(aa), 1.d-20, kind=8)
         if (real(p) > real(pc)) p = pc
         t = taup(p, x, ray_len, topp, bttm)
         if (t < t0) then
            t0 = t
            p0 = p
         endif
      enddo
      tmin(i) = t0
   enddo
   end subroutine trav_fk

   function taup(p, x, ray_len, topp, bttm) result(taup1)
      implicit none
      complex(kind=8), intent(in) :: p
      real(kind=8), intent(in) :: x, ray_len(2,*)
      integer, intent(in) :: topp, bttm
      complex(kind=8) :: taup1, pp
      integer i
      taup1 = p * x
      pp = p * p
      do i = topp, bttm
         taup1 = taup1 + sqrt(slowness(i) - pp) * ray_len(1, i) + sqrt(slowness(i) - pp) * ray_len(2, i)
      enddo
   end function taup

   FUNCTION dtdp(x, p, ray_len, topp, bttm) result(dtdp1)
      ! define d(tau)/dp
      IMPLICIT NONE
      real(kind=8), intent(in) :: x, ray_len(2,*)
      complex(kind=8), intent(in) :: p
      integer, intent(in) :: topp, bttm
      INTEGER j
      complex(kind=8) :: dtdp1, pp
      pp = p*p
      dtdp1 = cmplx(0.d0, 0.d0, kind=8) ! Ensure kind matches complex declaration
      DO j = topp, bttm
         dtdp1 = dtdp1 - ray_len(1,j) / SQRT(slowness(j) - pp) &
                  - ray_len(2,j) / SQRT(slowness(j) - pp)
      ENDDO
      dtdp1 = x + p * dtdp1
   END function dtdp



   SUBROUTINE findp0(x,p0,ray_len,topp,bttm)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     solve d(tau)/dp = 0 for p0, tau=p*x+eta*z
!     input:
!        x --- distance
!        p0 -- the largest possible p
!     output: p0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   IMPLICIT NONE
   real*8 :: ray_len(2,nlay)
   integer topp,bttm
   REAL*8 :: ZERO_0, dtdp0, x
   COMPLEX*16 :: p0, p1, p2
   ZERO_0 = 1.d-7
   p1 = cmplx(0.d0,aimag(p0),kind(1d0))
   DO WHILE ( p1.NE.p0 )
      p2 = p0
      p0 = 0.5d0*(p1+p2)
      dtdp0 = dtdp(x,p0,ray_len,topp,bttm)
      IF ( abs(dtdp0).LT.ZERO_0 .OR. p0.EQ.p1 .OR. p0.EQ.p2 ) RETURN
      IF( dtdp0 .GT. 0. ) THEN
         p1 = p0
         p0 = p2
      END IF
   ENDDO
   END subroutine findp0


end module wave_travel
