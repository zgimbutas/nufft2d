cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a FreeBSD license
cc (see COPYING in this directory). 
c
c  NUFFT 1.3 release notes:
c
c  These codes are asymptotically fast (O(N log N)), but not optimized.
c
c  1) We initialize the FFT on every call.
c
c  2) We do not precompute the exponentials involved in "fast Gaussian
c  gridding".
c
c  3) We do not block structure the code so that irregularly placed points
c  are interpolated (gridded) in a cache-aware fashion.
c
c  4) We use the Netlib FFT library (www.netlib.org) 
c     rather than the state of the art FFTW package (www.fftw.org).
c
c  Different applications have different needs, and we have chosen
c  to provide the simplest code as a reasonable efficient template.
c
c**********************************************************************
      subroutine nufft2d1f90(nj,xj,yj,cj, iflag,eps, ms,mt,fk,ier)
      implicit real *8 (a-h,o-z)
      real*8 xj(nj),yj(nj),xc(-147:147),yc(-147:147)
      complex*16 cj(nj),fk(-ms/2:(ms-1)/2,-mt/2:(mt-1)/2)
c ----------------------------------------------------------------------
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cz,ccj,zz
c ----------------------------------------------------------------------
      complex *16, allocatable :: fm1(:)
      complex *16, allocatable :: fm2(:)
      real *8, allocatable :: cw1(:)
      real *8, allocatable :: cw2(:)
      complex *16, allocatable :: fw(:,:)
      complex *16, allocatable :: tmp(:)
      real *8, allocatable :: iw1(:)
      real *8, allocatable :: iw2(:)
c----------------------------------------------------------------------
c
c     if (iflag .ge. 0) then
c
c               1  nj
c     fk(k1,k2) = -- SUM cj(j) exp(+i (k1,k2) * (xj(j, yj(j)) )   
c               nj j=1                     
c                                   for -ms/2 <= k1 <= (ms-1)/2
c                                   for -mt/2 <= k1 <= (mt-1)/2
c
c    else 
c
c               1  nj
c     fk(k1,k2) = -- SUM cj(j) exp(-i (k1,k2) * (xj(j, yj(j)) )   
c               nj j=1                     
c                                   for -ms/2 <= k1 <= (ms-1)/2
c                                   for -mt/2 <= k1 <= (mt-1)/2
c
c     References:
c
c     [DR] Fast Fourier transforms for nonequispaced data,
c          A. Dutt and V. Rokhlin, SIAM J. Sci. Comput. 14, 
c          1368-1383, 1993.
c
c     [GL] Accelerating the Nonuniform Fast Fourier Transform,
c          L. Greengard and J.-Y. Lee, SIAM Review 46, 443-454 (2004).
c
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of sources   (integer)
c     xj,yj  location of sources (real *8)
c     cj     strengths of sources (complex *16)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
c     ms     number of Fourier modes computed (-ms/2 to (ms-1)/2 )
c     mt     number of Fourier modes computed (-mt/2 to (mt-1)/2 )
c
c     OUTPUT:
c
c     fk     Fourier transform values (2D complex *16 array)
c     ier    error return code
c   
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c     The type 1 NUFFT proceeds in three steps (see [GL]).
c
c     1) spread data to oversampled regular mesh using convolution with
c        a Gaussian 
c     2) compute FFT on uniform mesh
c     3) deconvolve each Fourier mode independently
c          (mutiplying by Fourier transform of Gaussian).
c
c ----------------------------------------------------------------------
c
c     The oversampled regular mesh is defined by 
c
c     nf1 = rat*ms  points, where rat is the oversampling ratio.
c     nf2 = rat*mt  points, where rat is the oversampling ratio.
c       
c     For simplicity, we set  
c
c         rat = 2 for eps > 1.0d-11
c         rat = 3 for eps <= 1.0d-11.
c
c     The Gaussian used for convolution is:
c
c        g(x) = exp(-x^2 / 4tau) 
c
c     It can be shown [DR] that the precision eps is achieved when
c
c     nspread = floor(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
c     and tau is chosen as
c
c     tau = pi*lambda/(ms**2)
c     lambda = nspread/(rat(rat-0.5)).
c
c     Note that the Fourier transform of g(x) is
c
c     G(s) = exp(-s^2 tau) = exp(-pi*lambda s^2/ms^2)
c
c
c ----------------------------------------------------------------------
c     Fast Gaussian gridding is based on the following observation.
c
c     Let hx = 2*pi/nf1.    
c
c     For source (xj), let kb*hx denote the closest grid point and
c     let  kx*hx be a regular grid point within the spreading
c     distance. We can write
c
c     (xj) - kx*hx = kb*hx + diff*hx - kx*hx = diff*hx - (kx-kb)*hx
c
c     where diff = (xj)/hx - kb.
c
c     Let t1 = hx*hx/(4 tau) = pi/(nf1*nf1)/lambda*ms*ms
c                            = pi/lambda/(rat*rat)
c
c     exp(-( (xj) -kx*hx)**2 / 4 tau)
c         = exp(-pi/lamb/rat^2 *(diff - (kx-kb))**2)
c         = exp(-t1 *(diff - (kx-kb))**2)
c         = exp(-t1*diff**2) * exp(2*t1*diff)**k * exp(-t1*k**2)
c           where k = kx-kb.
c 
c     In 2D/3D, this decomposition is carried out in each dimension.
c************************************************************************
c     -------------------------------
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c     -------------------------------
      ier = 0
      if ((eps.lt.1d-33) .or. (eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif
      nspread = floor(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0) + 1
      nf1 = rat*ms
      nf2 = rat*mt
      if (2*nspread.gt.min(nf1,nf2)) then
         ii = next235(2d0*nspread)
         nf1 = max(nf1, ii) 
         nf2 = max(nf2, ii)
      endif

c
c     lambda (described above) = nspread/(rat*(rat-0.5d0)) 
c     It is more convenient to define r2lamb = rat*rat*lambda
c
      r2lamb = rat*rat * nspread / (rat*(rat-.5d0))
      hx = 2*pi/nf1
      hy = 2*pi/nf2

c     -------------------------------
c     Compute workspace size and allocate
c     -------------------------------
      allocate( fw(0:(nf1-1), 0:(nf2-1)) )
      allocate( tmp(0:max(nf1,nf2)-1) )
        
      allocate( fm1(-ms/2 : (ms-1)/2) )
      allocate( fm2(-mt/2 : (mt-1)/2) )
      allocate( cw1(1:(nspread+1)) )
      allocate( cw2(1:(nspread+1)) )
      allocate( iw1(4*nf1+15) )
      allocate( iw2(4*nf2+15) )
c       
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one of the terms in fast Gaussian gridding.
c     ---------------------------------------------------------------
      t1 = pi/r2lamb
      t2 = pi/r2lamb
      do k1 = 1, nspread
         cw1(k1) = exp(-t1*k1**2)
         cw2(k1) = exp(-t2*k1**2)
      enddo
      call zffti(nf1,iw1)
      call zffti(nf2,iw2)

c     ---------------------------------------------------------------
c     Precompute deconvolution constants
c     ---------------------------------------------------------------
      tau1 = pi * r2lamb / dble(nf1)**2
      cross1 = 1.0/sqrt(r2lamb)
      do k1 = -ms/2, (ms-1)/2
         fm1(k1)=exp(tau1*dble(k1)**2)*cross1
      enddo

      tau2 = pi * r2lamb / dble(nf2)**2
      cross1 = 1.0/sqrt(r2lamb)
      do k2 = -mt/2, (mt-1)/2
         fm2(k2)=exp(tau2*dble(k2)**2)*cross1
      enddo
       
c     ---------------------------------------------------------------
c     Initialize fine grid data to zero.
c     ---------------------------------------------------------------
      do k2 = 0, nf2-1
         do k1 = 0, nf1-1
            fw(k1,k2) = dcmplx(0.0d0, 0.0d0)
         enddo
      enddo
c       
c     ---------------------------------------------------------------
c     Loop over sources (1,...,nj)
c
c     1. find closest mesh point (with periodic wrapping if necessary)
c     2. spread source data onto nearest nspread**2 grid points
c        using fast Gaussian gridding.
c
c     The following is a little hard to read because it takes
c     advantage of fast gridding and optimized to minimize the 
c     the number of multiplies in the inner loops.
c
c ---------------------------------------------------------------
      t1 = pi/r2lamb
      t2 = pi/r2lamb
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k1,jb1,diff1,k2,jb2,diff2)
C$OMP$PRIVATE(jb1d,jb1u,zz,cz,cross,cross1,istart,jstart,xc,yc,ccj)
C$OMP$reduction(+:fw)
      do j = 1, nj
         ccj = cj(j)/dble(nj)

         jb1 = floor((xj(j))/hx)
         diff1 = (xj(j))/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1

         jb2 = floor((yj(j))/hy)
         diff2 = (yj(j))/hy - jb2
         jb2 = mod(jb2, nf2)
         if (jb2.lt.0) jb2=jb2+nf2

         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = cw1(k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = cw1(k1)*cross
         enddo
            
         yc(0) = exp(-t2*diff2**2)
         cross = yc(0)
         cross1 = exp(2d0*t2 * diff2)
         do k1 = 1, nspread
            cross = cross * cross1
            yc(k1) = cw2(k1)*cross
         enddo
         cross = yc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            yc(-k1) = cw2(k1)*cross
         enddo

         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k2 = -nspread+1, nspread
            jstart = jb2+k2
            if(jstart.lt.0) then
               jstart = jstart + nf2
            elseif (jstart.ge.nf2) then
               jstart = jstart - nf2
            endif
            cz = yc(k2)*ccj
            do k1 = -nspread+1, -jb1d-1
               istart = (jb1+k1+nf1)
	       zz = xc(k1)*cz
               fw(istart,jstart) = fw(istart,jstart) + zz
            enddo
            do k1 = -jb1d, jb1u
               istart = (jb1+k1)
	       zz = xc(k1)*cz
               fw(istart,jstart) = fw(istart,jstart) + zz
            enddo
            do k1 = jb1u+1, nspread
               istart = (jb1+k1-nf1)
	       zz = xc(k1)*cz
               fw(istart,jstart) = fw(istart,jstart) + zz
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
c

c     ---------------------------------------------------------------
c     Compute 2D FFT and carry out deconvolution.
c 
c     ---------------------------------------------------------------
c
      do k2 = 0, nf2-1 
         if (iflag .ge. 0) then
            call zfftb(nf1,fw(0,k2),iw1)
         else
            call zfftf(nf1,fw(0,k2),iw1)
         endif
      enddo
c
      do k1 = 0, nf1-1 
         do ii = 0, nf2-1
         tmp(ii) = fw(k1,ii)
         enddo
         if (iflag .ge. 0) then
            call zfftb(nf2,tmp(0),iw2)
         else
            call zfftf(nf2,tmp(0),iw2)
         endif
         do ii = 0, nf2-1
         fw(k1,ii) = tmp(ii)
         enddo
      enddo
c
c
      do k2 = -mt/2, (mt-1)/2
         j2 = mod(k2, nf2)
         if (j2 .lt. 0) j2=j2+nf2
         do k1 = -ms/2, (ms-1)/2
            j1 = mod(k1, nf1)
            if (j1 .lt. 0) j1=j1+nf1
            fk(k1,k2) = fm1(k1)*fm2(k2)*fw(j1,j2)
         enddo
      enddo
  
      return
      end
c
c
c
c
c
************************************************************************
      subroutine nufft2d2f90(nj,xj,yj,cj, iflag,eps, ms,mt,fk,ier)
      implicit real *8 (a-h,o-z)
      real*8 xj(nj),yj(nj),xc(-147:147),yc(-147:147)
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(nj),fk(-ms/2:(ms-1)/2,-mt/2:(mt-1)/2),zz,cz
c ----------------------------------------------------------------------
      complex *16, allocatable :: fm1(:)
      complex *16, allocatable :: fm2(:)
      real *8, allocatable :: cw1(:)
      real *8, allocatable :: cw2(:)
      complex *16, allocatable :: fw(:,:)
      complex *16, allocatable :: tmp(:)
      real *8, allocatable :: iw1(:)
      real *8, allocatable :: iw2(:)
c----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c         (ms-1)/2 (mt-1)/2
c    cj(j) = SUM    SUM    fk(k1,k2) exp(+i (k1,k2) * (xj(j),yj(j)) )   
c         k1=-ms/2 k2=-mt/2             
c                                           for j = 1,...,nj
c     else
c
c         (ms-1)/2 (mt-1)/2
c    cj(j) = SUM    SUM    fk(k1,k2) exp(-i (k1,k2) * (xj(j),yj(j)) )   
c         k1=-ms/2 k2=-mt/2             
c                                           for j = 1,...,nj
c
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of output values   (integer)
c     xj,yj  location of output values (real *8 arrays)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
c     ms     number of Fourier modes given (1st index)  [-ms/2:(ms-1)/2]
c     mt     number of Fourier modes given (2nd index)  [-mt/2:(mt-1)/2]
c     fk     Fourier coefficient values (complex *16 2D array)
c
c     OUTPUT:
c
c     cj     output values (complex *16 array)
c     ier    error return code
c   
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c
c     The type 2 algorithm proceeds in three steps (see [GL]).
c
c     1) deconvolve (amplify) each Fourier mode first 
c     2) compute inverse FFT on uniform fine grid
c     3) spread data to regular mesh using Gaussian
c
c
c     See subroutine nufft2d1f90
c     for more comments on fast gridding and parameter selection.
c
************************************************************************
c     -------------------------------
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c     -------------------------------
      ier = 0
      if ((eps.lt.1d-33).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif
      nspread = floor(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0) + 1
      nf1 = rat*ms
      nf2 = rat*mt
      if (2*nspread.gt.min(nf1,nf2)) then
         ii = next235(2d0*nspread)
         nf1 = max(nf1, ii)
         nf2 = max(nf2, ii)
      endif
c
c     lambda (described above) = nspread/(rat*(rat-0.5d0)) 
c     It is more convenient to define r2lamb = rat*rat*lambda
c
c     -------------------------------

      r2lamb = rat*rat * nspread / (rat*(rat-.5d0))
      hx = 2*pi/nf1
      hy = 2*pi/nf2
c
c     -------------------------------
c     Compute workspace size and allocate
c
c     nw1 for fine grid, nw2 for second index, 
c     nw3 for max length to hold 1D FFT data.
c     -------------------------------
      allocate( fw(0:(nf1-1), 0:(nf2-1)) )
      allocate( tmp(0:max(nf1,nf2)-1) )
        
      allocate( fm1(-ms/2 : (ms-1)/2) )
      allocate( fm2(-mt/2 : (mt-1)/2) )
      allocate( cw1(1:(nspread+1)) )
      allocate( cw2(1:(nspread+1)) )
      allocate( iw1(4*nf1+15) )
      allocate( iw2(4*nf2+15) )
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one of the terms in fast Gaussian gridding.
c     ---------------------------------------------------------------
      t1 = pi/r2lamb
      t2 = pi/r2lamb
      do k1 = 1, nspread
         cw1(k1) = exp(-t1*k1**2)
         cw2(k1) = exp(-t2*k1**2)
      enddo
      call zffti(nf1,iw1)
      call zffti(nf2,iw2)

c     ---------------------------------------------------------------
c     Precompute deconvolution constants
c     ---------------------------------------------------------------
      tau1 = pi * r2lamb / dble(nf1)**2
      cross1 = 1.0/sqrt(r2lamb)
      do k1 = -ms/2, (ms-1)/2
         fm1(k1)=exp(tau1*dble(k1)**2)*cross1
      enddo

      tau2 = pi * r2lamb / dble(nf2)**2
      cross1 = 1.0/sqrt(r2lamb)
      do k2 = -mt/2, (mt-1)/2
         fm2(k2)=exp(tau2*dble(k2)**2)*cross1
      enddo
       
c     ---------------------------------------------------------------
c     Initialize fine grid data to zero.
c     ---------------------------------------------------------------
      do k2 = 0, nf2-1
         do k1 = 0, nf1-1
            fw(k1,k2) = dcmplx(0.0d0, 0.0d0)
         enddo
      enddo

c
c     ---------------------------------------------------------------
c     Deconvolve and compute inverse 2D FFT
c     ---------------------------------------------------------------
      do k2 = -mt/2, (mt-1)/2
         j2 = mod(k2, nf2)
         if (j2 .lt. 0) j2=j2+nf2
         do k1 = -ms/2, (ms-1)/2
            j1 = mod(k1, nf1)
            if (j1 .lt. 0) j1=j1+nf1
            fw(j1,j2) = fm1(k1)*fm2(k2)*fk(k1,k2)
         enddo
      enddo
c  
      do k2 = 0, nf2-1 
         if (iflag .ge. 0) then
            call zfftb(nf1,fw(0,k2),iw1)
         else
            call zfftf(nf1,fw(0,k2),iw1)
         endif
      enddo
c
      do k1 = 0, nf1-1 
         do ii = 0, nf2-1
         tmp(ii) = fw(k1,ii)
         enddo
         if (iflag .ge. 0) then
            call zfftb(nf2,tmp(0),iw2)
         else
            call zfftf(nf2,tmp(0),iw2)
         endif
         do ii = 0, nf2-1
         fw(k1,ii) = tmp(ii)
         enddo
      enddo
c
c
c
c     ---------------------------------------------------------------
c     Loop over target points (1,...,nj)
c
c     1. find closest mesh point (with periodic wrapping if needed)
c     2. spread source data onto nearest nspread grid points
c        using fast Gaussian gridding.
c
c     The following is a little hard to read because it takes
c     advantage of fast gridding and optimized to minimize the 
c     the number of multiplies in the inner loops.
c     ---------------------------------------------------------------
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k1,k2,jb1,diff1,jb2,diff2)
C$OMP$PRIVATE(jb1d,jb1u,zz,cz,cross,cross1,ii,xc,yc)
      do j = 1, nj
         cj(j) = dcmplx(0d0,0d0)
        
         jb1 = floor((xj(j))/hx)
         diff1 = (xj(j))/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1

         jb2 = floor((yj(j))/hy)
         diff2 = (yj(j))/hy - jb2
         jb2 = mod(jb2, nf2)
         if (jb2.lt.0) jb2=jb2+nf2

         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = cw1(k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = cw1(k1)*cross
         enddo
            
         yc(0) = exp(-t2*diff2**2)
         cross = yc(0)
         cross1 = exp(2d0*t2 * diff2)
         do k1 = 1, nspread
            cross = cross * cross1
            yc(k1) = cw2(k1)*cross
         enddo
         cross = yc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            yc(-k1) = cw2(k1)*cross
         enddo
         
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k2 = -nspread+1, nspread
            ii = jb2+k2
            if (ii.lt.0) then
               ii = ii + nf2
            elseif (ii.ge.nf2) then
               ii = ii - nf2
            endif
            cz = dcmplx(0d0, 0d0)
            do k1 = -nspread+1, -jb1d-1
               zz = fw(k1+nf1+jb1,ii)
               cz = cz + xc(k1)*zz
            enddo
            do k1 = -jb1d, jb1u
               zz = fw(k1+jb1,ii)
               cz = cz + xc(k1)*zz
            enddo
            do k1 = jb1u+1, nspread
               zz = fw(k1-nf1+jb1,ii)
               cz = cz + xc(k1)*zz
            enddo
            cj(j) = cj(j) + yc(k2)*cz
         enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine nufft2d3f90(nj,xj,yj,cj, iflag,eps, nk,sk,tk,fk,ier)
      implicit real *8 (a-h,o-z)
      real*8 xj(nj),yj(nj),sk(nk),tk(nk)
      real*8 xc(-147:147),yc(-147:147)
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(nj),fk(nk),zz,zz2,zc,zcs,cz,ccj
c ----------------------------------------------------------------------
      complex *16, allocatable :: fm1(:)
      complex *16, allocatable :: fm2(:)
      real *8, allocatable :: cw1(:)
      real *8, allocatable :: cw2(:)
      complex *16, allocatable :: fw(:,:)
      complex *16, allocatable :: tmp(:)
      real *8, allocatable :: iw1(:)
      real *8, allocatable :: iw2(:)
c ----------------------------------------------------------------------
c
c     if (iflag .ge. 0) then
c
c              nj
c     fk(k) = SUM cj(j) exp(+i s(k) xj(j)) exp(+/-i t(k) yj(j))
c             j=1                     
c                                           for  k = 1, ..., nk
c     else
c
c              nj
c     fk(k) = SUM cj(j) exp(-i s(k) xj(j)) exp(+/-i t(k) yj(j))
c             j=1                     
c                                           for  k = 1, ..., nk
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of sources   (integer)
c     xj,yj  location of sources (real *8 arrays)
c     cj     strengths of sources (complex *16 array)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
c     nk     number of (noninteger) Fourier modes computed
c     sk,tk  k-values (locations) of desired Fourier modes
c
c     OUTPUT:
c
c     fk     Fourier transform values (double complex array)
c     ier    error return code
c   
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c
c     References:
c
c     [DR] Fast Fourier transforms for nonequispaced data,
c          A. Dutt and V. Rokhlin, SIAM J. Sci. Comput. 14, 
c          1368-1383, 1993.
c
c     [LG] The type 3 nonuniform FFT and its applications
c          J.-Y. Lee and L. Greengard, J. Comput. Phys. 206, 1-5 (2005).
c
c     The algorithm is essentially a concatenation of the 
c     type 1 and 2 transforms.
c
c     1) Gaussian gridding of strengths cj(j) centered at xj 
c        to create f_\tau(n \Delta_x)  (in notation of [LG])
c     2) Deconvolve each regular grid mode 
c        to create f^{-\sigma}_\tau(n \Delta_x)  (in notation of [LG])
c     3) compute FFT on uniform mesh
c        to create F^{-\sigma}_\tau(m \Delta_s)  (in notation of [LG])
c     4) Gaussian gridding to irregular frequency points
c        to create F_\tau(s_k)  (in notation of [LG])
c     5) Deconvolution of  result 
c        to create F(s_k)  (in notation of [LG])
c
************************************************************************
      ier = 0
      if ((eps.lt.1d-33).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
c     -------------------------------
c     Compute the ranges of {xj}/{sk} 
c     -------------------------------
      t1 = xj(1)
      t2 = xj(1)
      do j = 2, nj
         if (xj(j).gt.t2) then
             t2=xj(j)
         else if (xj(j).lt.t1) then
             t1=xj(j)
         endif
      enddo
      xb = (t1+t2) / 2d0
      xm = max(t2-xb,-t1+xb)  ! max(abs(t2-xb),abs(t1-xb))
c
      t1 = sk(1)
      t2 = sk(1)
      do k1 = 2, nk
         if (sk(k1).gt.t2) then
             t2=sk(k1)
         else if (sk(k1).lt.t1) then
             t1=sk(k1)
         endif
      enddo
      sb = (t1+t2) / 2d0
      sm = max(t2-sb,-t1+sb)
c
      t1 = yj(1)
      t2 = yj(1)
      do j = 2, nj
         if (yj(j).gt.t2) then
             t2=yj(j)
         else if (yj(j).lt.t1) then
             t1=yj(j)
         endif
      enddo
      yb = (t1+t2) / 2d0
      ym = max(t2-yb,-t1+yb)
c
      t1 = tk(1)
      t2 = tk(1)
      do k1 = 2, nk
         if (tk(k1).gt.t2) then
             t2=tk(k1)
         else if (tk(k1).lt.t1) then
             t1=tk(k1)
         endif
      enddo
      tb = (t1+t2) / 2d0
      tm = max(t2-tb,-t1+tb)
c
c     -------------------------------
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c     -------------------------------
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif
      nspread = floor(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0) + 1
c
      t1 = 2d0/pi * xm*sm
      t2 = 2d0/pi * ym*tm
      nf1 = next235(rat*max(rat*t1+2*nspread,2*nspread/(rat-1)))
      nf2 = next235(rat*max(rat*t2+2*nspread,2*nspread/(rat-1)))
      rat1 = (sqrt(nf1*t1+nspread**2)-nspread)/t1
      rat2 = (sqrt(nf2*t2+nspread**2)-nspread)/t2
c
      r2lamb1 = rat1*rat1 * nspread / (rat1*(rat1-.5d0))
      r2lamb2 = rat2*rat2 * nspread / (rat2*(rat2-.5d0))
      hx = pi/(rat1*sm)
      hs = 2d0*pi/dble(nf1)/hx
      hy = pi/(rat2*tm)
      ht = 2d0*pi/dble(nf2)/hy
c
      kmax1 = floor(nf1*(r2lamb1-nspread)/r2lamb1+.1d0)
      kmax2 = floor(nf2*(r2lamb2-nspread)/r2lamb2+.1d0)
c
c     -------------------------------
c     Compute workspace size and allocate
c     -------------------------------
c  
      allocate( fw(0:(nf1-1), 0:(nf2-1)) )
      allocate( tmp(0:max(nf1,nf2)-1) )        
      allocate( cw1(1:(nspread+1)) )
      allocate( cw2(1:(nspread+1)) )
      allocate( iw1(4*nf1+15) )
      allocate( iw2(4*nf2+15) )
      allocate( fm1(-kmax1:kmax1) )
      allocate( fm2(-kmax2:kmax2) )
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one of the terms in fast Gaussian gridding.
c     ---------------------------------------------------------------
      t1 = pi/r2lamb1
      t2 = pi/r2lamb2
      do k1 = 1, nspread
         cw1(k1) = exp(-t1*k1**2)
         cw2(k1) = exp(-t2*k1**2)
      enddo
      call zffti(nf1,iw1)
      call zffti(nf2,iw2)

c     ---------------------------------------------------------------
c     Precompute deconvolution constants
c     ---------------------------------------------------------------
      tau1 = pi * r2lamb1 / dble(nf1)**2
      cross1 = 1.0/(r2lamb1)
      do k1 = -kmax1, kmax1
         fm1(k1)=exp(tau1*dble(k1)**2)*cross1
      enddo

      tau2 = pi * r2lamb2 / dble(nf2)**2
      cross1 = 1.0/(r2lamb2)
      do k1 = -kmax2, kmax2
         fm2(k1)=exp(tau2*dble(k1)**2)*cross1
      enddo
c
c
c     ---------------------------------------------------------------
c     Initialize fine grid data to zero.
c     ---------------------------------------------------------------
      do k2 = 0, nf2-1
         do k1 = 0, nf1-1
            fw(k1,k2) = dcmplx(0.0d0, 0.0d0)
         enddo
      enddo
c
c ---------------------------------------------------------------
c Step 1/5 - gridding as in type 1 transform.
c ---------------------------------------------------------------
c
      t1 = pi/r2lamb1
      t2 = pi/r2lamb2
      if (iflag .lt. 0) sb = -sb
      if (iflag .lt. 0) tb = -tb
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k1,jb1,diff1,k2,jb2,diff2)
C$OMP$PRIVATE(jb1d,jb1u,zz,cz,cross,cross1,istart,jstart,xc,yc,ccj,ang)
C$OMP$reduction(+:fw)
      do j = 1, nj
         ang = sb*xj(j) + tb*yj(j)
         ccj = dcmplx(cos(ang),sin(ang)) * cj(j)

         jb1 = floor((xj(j)-xb)/hx)
         diff1 = (xj(j)-xb)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1

         jb2 = floor((yj(j)-yb)/hy)
         diff2 = (yj(j)-yb)/hy - jb2
         jb2 = mod(jb2, nf2)
         if (jb2.lt.0) jb2=jb2+nf2

         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = cw1(k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = cw1(k1)*cross
         enddo
            
         yc(0) = exp(-t2*diff2**2)
         cross = yc(0)
         cross1 = exp(2d0*t2 * diff2)
         do k1 = 1, nspread
            cross = cross * cross1
            yc(k1) = cw2(k1)*cross
         enddo
         cross = yc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            yc(-k1) = cw2(k1)*cross
         enddo

         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k2 = -nspread+1, nspread
            jstart = jb2+k2
            if(jstart.lt.0) then
               jstart = jstart + nf2
            elseif (jstart.ge.nf2) then
               jstart = jstart - nf2
            endif
            cz = yc(k2)*ccj
            do k1 = -nspread+1, -jb1d-1
               istart = (jb1+k1+nf1)
	       zz = xc(k1)*cz
               fw(istart,jstart) = fw(istart,jstart) + zz
            enddo
            do k1 = -jb1d, jb1u
               istart = (jb1+k1)
	       zz = xc(k1)*cz
               fw(istart,jstart) = fw(istart,jstart) + zz
            enddo
            do k1 = jb1u+1, nspread
               istart = (jb1+k1-nf1)
	       zz = xc(k1)*cz
               fw(istart,jstart) = fw(istart,jstart) + zz
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
      if (iflag .lt. 0) sb = -sb
      if (iflag .lt. 0) tb = -tb
c
c     ---------------------------------------------------------------
c     Step 2/5: Deconvolve (amplify) as in Type 2 transform.
c     Step 3/5: Compute FFT
c     ---------------------------------------------------------------

      do k2 = -kmax2, kmax2
         j2 = mod(k2, nf2)
         if (j2 .lt. 0) j2=j2+nf2
         do k1 = -kmax1, kmax1
            j1 = mod(k1, nf1)
            if (j1 .lt. 0) j1=j1+nf1
            fw(j1,j2) = fm1(k1)*fm2(k2)*fw(j1,j2)
         enddo
      enddo
c  
      do k2 = 0, nf2-1 
         if (iflag .ge. 0) then
            call zfftb(nf1,fw(0,k2),iw1)
         else
            call zfftf(nf1,fw(0,k2),iw1)
         endif
      enddo
c
      do k1 = 0, nf1-1 
         do ii = 0, nf2-1
         tmp(ii) = fw(k1,ii)
         enddo
         if (iflag .ge. 0) then
            call zfftb(nf2,tmp(0),iw2)
         else
            call zfftf(nf2,tmp(0),iw2)
         endif
         do ii = 0, nf2-1
         fw(k1,ii) = tmp(ii)
         enddo
      enddo

c 
c     ---------------------------------------------------------------
c     Step 4/5: Gaussian gridding to irrregular points.
c     Step 5/5: Final deconvolution
c     ---------------------------------------------------------------
c
      t1 = pi/r2lamb1
      t2 = pi/r2lamb2
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k1,k2,jb1,diff1,jb2,diff2)
C$OMP$PRIVATE(jb1d,jb1u,zz,cz,cross,cross1,ii,xc,yc)
      do j = 1, nk
         fk(j) = dcmplx(0d0,0d0)

         jb1 = floor((sk(j)-sb)/hs)
         diff1 = (sk(j)-sb)/hs - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1

         jb2 = floor((tk(j)-tb)/ht)
         diff2 = (tk(j)-tb)/ht - jb2
         jb2 = mod(jb2, nf2)
         if (jb2.lt.0) jb2=jb2+nf2

         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = cw1(k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = cw1(k1)*cross
         enddo
            
         yc(0) = exp(-t2*diff2**2)
         cross = yc(0)
         cross1 = exp(2d0*t2 * diff2)
         do k1 = 1, nspread
            cross = cross * cross1
            yc(k1) = cw2(k1)*cross
         enddo
         cross = yc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            yc(-k1) = cw2(k1)*cross
         enddo
         
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k2 = -nspread+1, nspread
            ii = jb2+k2
            if (ii.lt.0) then
               ii = ii + nf2
            elseif (ii.ge.nf2) then
               ii = ii - nf2
            endif
            cz = dcmplx(0d0, 0d0)
            do k1 = -nspread+1, -jb1d-1
               zz = fw(k1+nf1+jb1,ii)
               cz = cz + xc(k1)*zz
            enddo
            do k1 = -jb1d, jb1u
               zz = fw(k1+jb1,ii)
               cz = cz + xc(k1)*zz
            enddo
            do k1 = jb1u+1, nspread
               zz = fw(k1-nf1+jb1,ii)
               cz = cz + xc(k1)*zz
            enddo
            fk(j) = fk(j) + yc(k2)*cz
         enddo
      enddo
C$OMP END PARALLEL DO
        
      t1 = r2lamb1/(4d0*pi) * hx**2
      t2 = r2lamb2/(4d0*pi) * hy**2
      if (iflag .lt. 0) xb = -xb
      if (iflag .lt. 0) yb = -yb
      do j = 1, nk
         fk(j) = (exp(t1*(sk(j)-sb)**2+t2*(tk(j)-tb)**2))*fk(j)
         ang = (sk(j)-sb)*xb + (tk(j)-tb)*yb
         fk(j) = dcmplx(cos(ang),sin(ang)) * fk(j)
      enddo

      return
      end
