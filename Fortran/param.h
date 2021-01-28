      double precision position(0:18)
      double precision pi
      parameter(pi=3.1415926535897932)
      double precision jumpeps
      parameter(jumpeps=2.d-5)
      double precision pace,err
      parameter(pace=0.05,err=0.2d0)
      integer ndimx,ndimy,ndimz,winnum
      parameter(ndimx=229,ndimy=301,winnum=15)
      integer looplength
      parameter(looplength=1500000)
