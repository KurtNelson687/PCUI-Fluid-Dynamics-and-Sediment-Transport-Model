!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include 'mkl_vsl.f90'

      subroutine random_normal(r, n, seed)

      USE MKL_VSL_TYPE
      USE MKL_VSL
 
      implicit none

      TYPE (VSL_STREAM_STATE) :: stream
 
      integer :: errcode
      integer :: i,j
      integer :: brng,method,seed,n
 
      double precision, dimension(n*3) :: rall ! all random numbers
      double precision, dimension(n,3) :: r ! reshaped random numbers
      double precision :: mean, std ! parameters of normal distribution

      mean = 0.D0
      std  = 1.D0
      brng = VSL_BRNG_MT19937
      method = VSL_RNG_METHOD_GAUSSIAN_ICDF
 
!     ***** Initializing *****
      errcode=vslnewstream( stream, brng,  seed )
 
!     ***** Generating *****
      errcode=vdrnggaussian( method, stream, n*3, rall, mean, std )
      r(:,1) = rall(1:n)
      r(:,2) = rall(n+1:2*n)
      r(:,3) = rall(2*n+1:3*n)
 
!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
 
!     ***** Printing results *****
!     print *,"Sample mean of normal distribution = ", s

      return
      end
