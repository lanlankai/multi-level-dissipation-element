      double precision pi
      parameter(pi=3.1415926535897932)
      double precision jumpeps
      parameter(jumpeps=2.d-5)
      double precision pace,err
      parameter(pace=0.05,err=0.2d0)
	  double precision dimx,dimy,dimz	! domain size
	  parameter(dimx=2*pi,dimy=2*pi,dimz=2*pi)
      integer ndimx,ndimy,ndimz,winnum	! domain dimensions and maximal window size
      parameter(ndimx=80,ndimy=80,ndimz=80,winnum=15)
