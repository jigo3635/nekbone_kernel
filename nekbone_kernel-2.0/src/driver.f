c-----------------------------------------------------------------------
      program kernel

c     Solve Au=w where A is SPD and is invoked by ax()
c
c      u - vector of length n
c      g - geometric factors for SEM operator 
c
c      Work arrays:  wk  
 
#ifdef XSMM
      USE :: LIBXSMM
      USE :: STREAM_UPDATE_KERNELS
#endif

      include 'SIZE'
      include 'INPUT'

      include 'mpif.h'

      parameter (lt=lx1*ly1*lz1*lelt)
      real w(lt),u(lt),gxyz(6*lt)          
      real ur(lt),us(lt),ut(lt),wk(lt)
      real flop_a, flop_t, dur_a, dur_t

      integer mpierror, myrank, procs
      integer iters

      iters=100

      call mpi_init(mpierror)
      call mpi_comm_size(mpi_comm_world, procs, mpierror)
      call mpi_comm_rank(mpi_comm_world, myrank, mpierror)

      call semhat(ah,wxm1,ch,dxm1,zgm1,bh,lx1-1) ! find GLL weights and pts
      call transpose(dxtm1,lx1,dxm1,lx1)         ! transpose D matrix
      call setup_g(gxyz)                         ! geometric factors

      n = lx1*ly1*lz1*lelt    
      call setup_u(u,n)                          ! fill u-vector

      call cpu_time(t0)

#ifdef XSMM
      call ax_xsmm(w,u,gxyz,ur,us,ut,wk,n,iters, dur_a) !return w=A*u
#else
      do i=1,iters
         call ax(w,u,gxyz,ur,us,ut,wk,n) !return w=A*u
      enddo
#endif

      call cpu_time(t1)
      time1 = t1-t0

!      flop_a = (15*n + 12*lx1*n)*1.e-6/time1
      flop_a = 1e-9*(12.*lx1-4)*n*iters/time1

      call mpi_barrier(mpi_comm_world,mpierror)

      call mpi_reduce(flop_a,flop_t,1,mpi_double,mpi_sum, 
     &                0,mpi_comm_world,mpierror) 
      

      call mpi_reduce(dur_a,dur_t,1,mpi_double,mpi_max, 
     &                0,mpi_comm_world,mpierror) 


      if (myrank .eq. 0) then

#ifndef XSMM
         write(6,1) "Initialization time: ",t0
         write(6,1) "Time in ax(): ",  time1
         write(6,1) "Averages Flops ", flop_t 
 1       format(a30,1e14.5) 

#else
!     Print Performance Summary and check results
         call performance(dur_t, iters*procs, lx1, ly1, lz1, lelt)
#endif
      endif

      call mpi_finalize(mpierror)

      stop
      end

c-----------------------------------------------------------------------
      subroutine transpose(a,lda,b,ldb)
c     Returns the transpone of vector b
      real a(lda,1),b(ldb,1)
 
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine setup_g(g)
c     fill g with the GLL weights calculated in semhat routine

      include 'SIZE'
      include 'INPUT'
      real g(6,lx1,ly1,lz1,lelt)
      integer e

      n = lx1*ly1*lz1*lelt


      do e=1,lelt
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         call rzero(g(1,i,j,k,e),6)
         g(1,i,j,k,e) = wxm1(i)*wxm1(j)*wxm1(k)
         g(4,i,j,k,e) = wxm1(i)*wxm1(j)*wxm1(k)
         g(6,i,j,k,e) = wxm1(i)*wxm1(j)*wxm1(k)
         g(6,i,j,k,e) = wxm1(i)*wxm1(j)*wxm1(k)
      enddo
      enddo
      enddo
      enddo

      return
      end
c-------------------------------------------------------------------------
      subroutine setup_u(u,n)
c     fill vector u  - "random"
      real u(n)

      do i=1,n
         arg  = 1.e9*(i*i)
         arg  = 1.e9*cos(arg)
         u(i) = sin(arg)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ax(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u

      include 'SIZE'
      include 'INPUT'

      parameter (lxyz=lx1*ly1*lz1)
      real w(lxyz,lelt),u(lxyz,lelt),gxyz(6,lxyz,lelt)
      parameter (lt=lx1*ly1*lz1*lelt)
      real ur(lt),us(lt),ut(lt),wk(lt)

      integer e

      do e=1,lelt                                ! ~
         call ax_e( w(1,e),u(1,e),gxyz(1,1,e)    ! w   = A  u
     $                             ,ur,us,ut,wk) !  L     L  L
      enddo                                      ! 


      return
      end
c-------------------------------------------------------------------------
      subroutine ax1(w,u,n)
      include 'SIZE'
      real w(n),u(n)
      real h2i
  
      h2i = (n+1)*(n+1)  
      do i = 2,n-1
         w(i)=h2i*(2*u(i)-u(i-1)-u(i+1))
      enddo
      w(1)  = h2i*(2*u(1)-u(2  ))
      w(n)  = h2i*(2*u(n)-u(n-1))

      return
      end
c-------------------------------------------------------------------------
      subroutine ax_e(w,u,g,ur,us,ut,wk) ! Local matrix-vector product
      include 'SIZE'
      include 'INPUT'

      parameter (lxyz=lx1*ly1*lz1)
      real w(lxyz),u(lxyz),g(6,lxyz)

      real ur(lxyz),us(lxyz),ut(lxyz),wk(lxyz)

      nxyz = lx1*ly1*lz1
      n    = lx1-1

      call local_grad3(ur,us,ut,u,n,dxm1,dxtm1)

      do i=1,nxyz
         wr = g(1,i)*ur(i) + g(2,i)*us(i) + g(3,i)*ut(i)
         ws = g(2,i)*ur(i) + g(4,i)*us(i) + g(5,i)*ut(i)
         wt = g(3,i)*ur(i) + g(5,i)*us(i) + g(6,i)*ut(i)
         ur(i) = wr
         us(i) = ws
         ut(i) = wt
      enddo

      call local_grad3_t(w,ur,us,ut,n,dxm1,dxtm1,wk)

      return
      end
c-------------------------------------------------------------------------
      subroutine local_grad3(ur,us,ut,u,n,D,Dt)
c     Output: ur,us,ut         Input:u,n,D,Dt
      real ur(0:n,0:n,0:n),us(0:n,0:n,0:n),ut(0:n,0:n,0:n)
      real u (0:n,0:n,0:n)
      real D (0:n,0:n),Dt(0:n,0:n)
      integer e

      m1 = n+1
      m2 = m1*m1

#ifndef XSMM
#ifdef SIMD
      call mxm_ism(D ,m1,u,m1,ur,m2)
      do k=0,n
         call mxm_ism(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm_ism(u,m2,Dt,m1,ut,m1)
#else
      call mxm(D ,m1,u,m1,ur,m2)
      do k=0,n
         call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u,m2,Dt,m1,ut,m1)
#endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad3_t(u,ur,us,ut,N,D,Dt,w)
c     Output: ur,us,ut         Input:u,N,D,Dt
      real u (0:N,0:N,0:N)
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      real w (0:N,0:N,0:N)
      integer e

      m1 = N+1
      m2 = m1*m1
      m3 = m1*m1*m1

#ifndef XSMM
#ifdef SIMD
      call mxm_ism(Dt,m1,ur,m1,u,m2)

      do k=0,N
         call mxm_ism(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
      enddo
      call add2(u,w,m3)

      call mxm_ism(ut,m2,D ,m1,w,m1)
      call add2(u,w,m3)
#else
      call mxm(Dt,m1,ur,m1,u,m2)

      do k=0,N
         call mxm(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
      enddo
      call add2(u,w,m3)

      call mxm(ut,m2,D ,m1,w,m1)
      call add2(u,w,m3)
#endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt0
      include 'SIZE'
      include 'INPUT'

      write(6,*) 'Exitting....'

      call exit(0)

      return
      end
c-----------------------------------------------------------------------


#ifdef XSMM
c-----------------------------------------------------------------------
      subroutine ax_xsmm(w,u,gxyz,ur2,us2,ut2,wk2,n2,iters, duration) ! Matrix-vector product: w=A*u

#ifdef XSMM
      USE :: LIBXSMM
      USE :: STREAM_UPDATE_KERNELS
#endif

      include 'SIZE'
      include 'INPUT'


      real w(lx1,ly1,lz1,lelt),u(lx1,ly1,lz1,lelt)
      real ur2(1), us2(1), ut2(1),wk2(1)

      real gxyz(2*ldim,lx1,ly1,lz1,lelt)

      real duration

      parameter (lt=lx1*ly1*lz1*lelt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)
      
      integer e
      
      INTEGER, PARAMETER :: T = KIND(0D0)
      REAL, PARAMETER :: alpha = 1, beta0 = 0, beta1 = 1
      
      REAL, allocatable, dimension(:,:,:), target :: ur, us, ut
      REAL, allocatable, target :: dx(:,:), dxt(:,:)
      REAL, ALLOCATABLE,TARGET,SAVE :: tm1(:,:,:), tm2(:,:,:),
     $     tm3(:,:,:), wk(:,:,:)

      TYPE(LIBXSMM_DMMFUNCTION) :: xmm1, xmm2, xmm3, xmm4, xmm5

      DOUBLE PRECISION :: max_diff
      INTEGER :: argc, m, n, k, routine, check
      INTEGER(8) :: i, j, ix, iy, iz, r, s, size0, size1, size, 
     $     start, it
      CHARACTER(32) :: argv
      s = lelt
      size = s
      ALLOCATE(ur(lx1,ly1,lz1), us(lx1,ly1,lz1), ut(lx1,ly1,lz1))
      ALLOCATE(wk(lx1,ly1,lz1))
      ALLOCATE(dx(lx1,lx1), dxt(ly1,ly1))

! Initialize LIBXSMM
      CALL libxsmm_init()

!  Initialize 
      do j = 1,ly1
         do i = 1, lx1
            dx(i,j)  = dxm1(i,j)
            dxt(i,j) = dxtm1(i,j)
         enddo
      enddo

c      call cpu_time(t1)

! streamed 
      
!      WRITE(*, "(A)") " Streamed... (specialized)"
      CALL libxsmm_dispatch(xmm1,lx1,ly1*lz1,lx1,alpha=alpha,beta=beta0)

      CALL libxsmm_dispatch(xmm2,lx1,ly1,ly1,alpha=alpha,beta=beta0)
      CALL libxsmm_dispatch(xmm3,lx1*ly1,lz1,lz1,alpha=alpha,beta=beta0)

      CALL libxsmm_dispatch(xmm4,lx1,ly1,ly1,alpha=alpha,beta=beta1)
      CALL libxsmm_dispatch(xmm5,lx1*ly1,lz1,lz1,alpha=alpha,beta=beta1)


      IF (libxsmm_available(xmm1).AND.libxsmm_available(xmm2) 
     $     .AND.libxsmm_available(xmm3).AND.libxsmm_available(xmm4)
     $     .AND.libxsmm_available(xmm5)) THEN         

         ALLOCATE(tm1(lx1,ly1,lz1), tm2(lx1,ly1,lz1), tm3(lx1,ly1,lz1))
         tm1 = 0; tm2 = 0; tm3 = 0
         start = libxsmm_timer_tick()

         DO it = 1,iters

         DO i = 1, lelt
C local_grad3
            CALL libxsmm_call(xmm1,  C_LOC(dx), C_LOC(u(1,1,1,i)),
     $           C_LOC(ur(1,1,1)))
            DO j = 1, ly1
               CALL libxsmm_call(xmm2, C_LOC(u(1,1,j,i)), 
     $              C_LOC(dxt), C_LOC(us(1,1,j)))
            END DO
            CALL libxsmm_call(xmm3, C_LOC(u(1,1,1,i)), C_LOC(dxt),
     $           C_LOC(ut(1,1,1)))

C local_grad3_t
            CALL libxsmm_call(xmm1,  C_LOC(dxt), C_LOC(ur(1,1,1)),
     $           C_LOC(wk(1,1,1)))
            DO j = 1, ly1
               CALL libxsmm_call(xmm4, C_LOC(us(1,1,j)), 
     $              C_LOC(dx), C_LOC(wk(1,1,j)))
            END DO
            
            CALL libxsmm_call(xmm5, C_LOC(ut(1,1,1)), C_LOC(dx),
     $           C_LOC(wk(1,1,1)))

            CALL stream_vector_copy(wk(1,1,1),w(1,1,1,i),lxyz)
         END DO
         END DO

         duration = libxsmm_timer_duration(start, libxsmm_timer_tick())

         DEALLOCATE(tm1, tm2, tm3, wk)
         DEALLOCATE(ur, us, ut)

C         IF (check.NE.0) max_diff = MAX(max_diff, 
C     $        validate(rx, ry, rz, cx, cy, cz))
      ELSE
         WRITE(*,*) "Could not build specialized function(s)!"
      END IF


! finalize LIBXSMM
      CALL libxsmm_finalize()

      return
      end


c-----------------------------------------------------------------------
      SUBROUTINE performance(duration, iters, m, n, k, size)
      DOUBLE PRECISION, INTENT(IN) :: duration
      INTEGER, INTENT(IN)    :: m, n, k
C      INTEGER(8), INTENT(IN) :: size
      integer size
      real T
      T = 8.0
      IF (0.LT.duration) THEN
         WRITE(*, 2) CHAR(9), "performance:", 
     $         (1D-9 * iters * size * m * n * k * (4*(m+n+k) - 4) / 
     $        duration),     " GFLOPS/s"
         WRITE(*, 2) CHAR(9), "bandwidth:  ", 
     $        (size*m*n*k*(2)*T*iters / (duration * ISHFT(1_8, 30)))
     $        , " GB/s"
      END IF
      WRITE(*, 2) CHAR(9), "duration:   ", (1D3 * duration), " ms"

 2    format(1A,A,F10.5,A)

      return 
      END

#endif
