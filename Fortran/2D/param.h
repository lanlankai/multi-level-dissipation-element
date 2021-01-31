      double precision pi
      parameter(pi=3.1415926535897932)
      double precision jumpeps
      parameter(jumpeps=2.d-5)
      double precision pace,err
      parameter(pace=0.05,err=0.2d0)
	  double precision dimx,dimy	! domain size
	  parameter(dimx=2*pi,dimy=2*pi)
      integer ndimx,ndimy,winnum	! domain dimensions and maximal window size
      parameter(ndimx=229,ndimy=301,winnum=15)
