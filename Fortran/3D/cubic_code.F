c this code is used to obtain the DE infomation from cubic(3D) scalar field

      program main
      implicit none
      include "param.h"
      
      character(len=50) sc
      
      integer i,j,k,ni,nj,nk,endfind,i0,j0,k0,c,m
      integer numends,numpairs    ! number of local extremal points and 0+ level DE
      integer force0,status
      integer pairing(1000000,3)   ! pairing(*,1) - index of maximal point in *-index DE
                                                               ! pairing(*,2) - index of minimal point in *-index DE
                                                               ! pairing(*,3) - volume of *-index DE
      integer, dimension(:,:,:,:),pointer:: numbering   ! numbering(I,J,K,1) - index of maximal point in DE which grid point (I,J,K) corresponds to
                                                                                                 ! numbering(I,J,K,2) - index of minimal point in DE which grid point (I,J,K) corresponds to
      
      double precision jumpratio
      double precision tmain0,tmain1      ! timing
      double precision reed
      double precision sign(5000000)  ! to record whether each extremal point is a maximal or a minimal
      double precision switch,posx,posy,posz,endx,endy,endz
      double precision subdvx(-1:0,-1:1,-1:1),subdvy(-1:1,-1:0,-1:1),
     . subdvz(-1:1,-1:1,-1:0)     ! partial derivatives
      double precision, dimension(:,:),pointer:: points    ! points(*,:) - x, y, z coordinate of *-index extremal point
      double precision, dimension(:,:,:),pointer:: psi          ! psi(I,J,K) - scalar value of grid point (I,J,K)
      double precision, dimension(:,:,:),pointer:: coordx,coordy,coordz   ! real coordinate of grid points
      
      common/large0/ psi
      common/large1/ pairing
      common/large2/ coordx,coordy,coordz
      common/large6/ numbering
      common/large5/ points      
      
      allocate(psi(ndimx,ndimy,ndimz),STAT=status)
      allocate(points(500000,3),STAT=status)
      allocate(numbering(ndimx,ndimy,ndimz,2),STAT=status)
      allocate(coordx(ndimx,ndimy,ndimz),STAT=status)
      allocate(coordy(ndimx,ndimy,ndimz),STAT=status)
      allocate(coordz(ndimx,ndimy,ndimz),STAT=status)
      
      do m = 1,1 ! a loop for different sample field. Set [do m = 1,1] if only one sample under  consideration.
      write(sc,*) m
      write(*,*) 'we are begining sample '//Trim(AdjustL(sc))
          
      ! ======= real coordinate calculation =======
      ! set domain size based on the cases under consideration in 'param.h'
      write(*,*) '======= begin coord calculation ======='
      do k =1,ndimz
          do j=1,ndimy
              do i=1,ndimx
                  coordx(i,j,k)=(i-1)*dimx/(ndimx-1)
                  coordy(i,j,k)=(j-1)*dimy/(ndimy-1)
                  coordz(i,j,k)=(k-1)*dimz/(ndimz-1)
              enddo
          enddo
      enddo          
      write(*,*) '======= end coord calculation ======='    
      
      ! ======= data input =======
      ! assign the data to be processed to psi(:,:,:), whose dimension needs to adapt to the cases under consideration
      write(*,*) '======= begin reading data ======='
      open(20,file='ocean_summer'//Trim(AdjustL(sc))//'.dat') ! data path
      do k = 1,ndimz
          do j = 1,ndimy
              do i = 1,ndimx
                  read(20,*) psi(i,j,k)
!        psi(i,j,k)=sin(coordx(i,j,k))*sin(coordy(i,j,k))
!     .*sin(coordz(i,j,k))+0.3*sin(11*coordx(i,j,k))
!     .*sin(11*coordy(i,j,k))*sin(11*coordz(i,j,k))             
              enddo
          enddo
      enddo
      close(20)
      write(*,*) '======= end reading data ======='          
          
      ! ======= variable initialization =======
      call CPU_TIME(tmain0)   
      numends = 0
      numbering = 0
      sign = 0.d0
      points = 0.d0
      jumpratio=0.d0           
      ! ======= 0+ level algorithm =======
      ! the algorithm starts from any spatial (grid) points
      do i = 2, ndimx - 1     ! 1 and ndimx are on the boundary     
      write(*,*) i, numends
      do j = 2, ndimy - 1
      do k = 2, ndimz - 1
      force0 = 0
      endfind = 0   ! endfind=1 means the extremal point has been found
      posx = 0         ! distance between searching point and grid point,ranging 0~1
      posy = 0
      posz = 0
      ni = i
      nj = j
      nk = k
          
      ! switch=1: searching along ascending direction of gradient trajectories
      switch=1.d0
      do while(endfind.eq.0)
      if((ni.eq.1).or.(ni.eq.ndimx).or.(nj.eq.1).or.(nj.eq.ndimy)
     .  .or.(nk.eq.1).or.(nk.eq.ndimz)) then
          goto 25
      endif   
     
      ! calculation of partial derivatives
       do i0=-1,0
       do j0=-1,1
       do k0=-1,1
       subdvx(i0,j0,k0)=
     .  (psi(ni+i0+1,nj+j0,nk+k0)-psi(ni+i0,nj+j0,nk+k0))/(dimx/ndimx)
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       do k0=-1,1
       subdvy(i0,j0,k0)=
     .  (psi(ni+i0,nj+j0+1,nk+k0)-psi(ni+i0,nj+j0,nk+k0))/(dimy/ndimy)
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,0
       subdvz(i0,j0,k0)=
     .  (psi(ni+i0,nj+j0,nk+k0+1)-psi(ni+i0,nj+j0,nk+k0))/(dimz/ndimz)
       enddo
       enddo
       enddo
       
       ! subroutine [endpoint]: to search along the gradient trajectories
       call endpoint(posx,posy,posz,endx,endy,endz,switch,
     . subdvx,subdvy,subdvz,endfind,jumpratio)
	 
       ! for local eddy event
       force0=force0+1                  
       if(force0.gt.500) then
       endfind=1
       goto 25       
       endif

       ! subroutine [numberfind]: if the extreme points corresponding to the grid point (i,j,k) have been found, match them
       if(endfind.eq.1) then
       call numberfind(1,i,j,k,numends,sign,ni+endx,nj+endy,nk+endz)
       endif
       
       ! subroutine [tellnew]: to tell which grid the current searching point is in
       call tellnew(ni,nj,nk,endx,endy,endz,posx,posy,posz)
       enddo
       
25    force0=0
       endfind=0
       posx=0
       posy=0
       posz=0
       ni=i
       nj=j
       nk=k

      ! switch=-1: searching along descending direction of gradient trajectories
       switch=-1.d0
       do while(endfind.eq.0)
       if((ni.eq.1).or.(ni.eq.ndimx).or.(nj.eq.1).or.(nj.eq.ndimy)
     .  .or.(nk.eq.1).or.(nk.eq.ndimz)) then
       goto 30
       endif
     
       ! calculation of partial derivatives
       do i0=-1,0
       do j0=-1,1
       do k0=-1,1
       subdvx(i0,j0,k0)=
     .  (psi(ni+i0+1,nj+j0,nk+k0)-psi(ni+i0,nj+j0,nk+k0))/(dimx/ndimx)
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       do k0=-1,1
       subdvy(i0,j0,k0)=
     .  (psi(ni+i0,nj+j0+1,nk+k0)-psi(ni+i0,nj+j0,nk+k0))/(dimy/ndimy)
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,0
       subdvz(i0,j0,k0)=
     .  (psi(ni+i0,nj+j0,nk+k0+1)-psi(ni+i0,nj+j0,nk+k0))/(dimz/ndimz)
       enddo
       enddo
       enddo
       
       call endpoint(posx,posy,posz,endx,endy,endz,switch,
     . subdvx,subdvy,subdvz,endfind,jumpratio)

       force0=force0+1
       if(force0.gt.500) then
       endfind=1
       goto 30
       endif
       
       if(endfind.eq.1) then
       call numberfind(2,i,j,k,numends,sign,ni+endx,nj+endy,nk+endz)
       endif

       call tellnew(ni,nj,nk,endx,endy,endz,posx,posy,posz)
       enddo
      
30    enddo
      enddo
      enddo
      ! ======= 0+ level algorithm ends =======
          
      write(*,*) 'jumpration=',jumpratio/ndimx/ndimy
      write(*,*) 'numends=',numends
      ! subroutine [process]: to match the maximal and minimal points in the same DE
      call process(numends,numpairs,sc)
      write(*,*) 'numpairs=',numpairs
      call CPU_TIME(tmain1)
      write(*,*) 'main duration=',tmain1-tmain0
        
      ! output the info. of 0+ DE for multi-level algorithm program
      open(12,file='./numends'//Trim(AdjustL(sc))//'.dat')
      write(12,*) numends
      close(12)
      
      open(12,file='./numpairs'//Trim(AdjustL(sc))//'.dat')
      write(12,*) numpairs
      close(12)   
      
      open(12,file='./pairing'//Trim(AdjustL(sc))//'.dat')
      do i=1,numpairs
          write(12,*) pairing(i,1),pairing(i,2),pairing(i,3)
      enddo
      close(12)   
      
      open(12,file='./points'//Trim(AdjustL(sc))//'.dat')
      write(12,*) ' variables = "x","y","z","sign" '
      do i=1,numends          
      write(12,'(4f12.6)') points(i,1),points(i,2),points(i,3),sign(i)
      enddo
      close(12)
      
      enddo
      ! end the loop of current sample
      
      end

      
      ! subroutine [tellnew]: to tell which grid the current searching point is in
      subroutine tellnew(ni,nj,nk,endx,endy,endz,posx,posy,posz)
      implicit none
      integer ni,nj,nk
      double precision endx,endy,endz,posx,posy,posz      
      ! end(x y z) relative coordinate ranging -1~1
      ! pos(x y z) simpliest coordinate ranging -0.5~0.5
      ! e.g: for point(1.7,2.3,2.4):  (posx,posy,posz)=(-0.3,0.3,0.4) ---> (ni,nj,nk)=(2,2,2)

      if(abs(endx).le.0.5.and.abs(endy).le.0.5.and.abs(endz).le.0.5) 
     . then
      ni=ni
      nj=nj
      nk=nk
      posx=endx
      posy=endy
      posz=endz
      elseif(endx.gt.0.5.and.abs(endy).le.0.5.and.abs(endz).le.0.5)then
      ni=ni+1
      nj=nj
      nk=nk
      posx=endx-1.d0            
      posy=endy
      posz=endz
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5)
     .  .and.(abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj
      nk=nk
      posx=endx+1.d0
      posy=endy
      posz=endz
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni
      nj=nj+1
      nk=nk
      posx=endx
      posy=endy-1.d0
      posz=endz
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni
      nj=nj-1
      nk=nk
      posx=endx
      posy=endy+1.d0
      posz=endz
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz
      elseif((abs(endx).le.0.5).and.(abs(endy).le.0.5).and.
     . (endz.gt.0.5)) then
      ni=ni
      nj=nj
      nk=nk+1
      posx=endx
      posy=endy
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(abs(endy).le.0.5).and.(endz.gt.0.5))
     . then
      ni=ni+1
      nj=nj
      nk=nk+1
      posx=endx-1.d0
      posy=endy
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5).and.
     . (endz.gt.0.5))then
      ni=ni-1
      nj=nj
      nk=nk+1
      posx=endx+1.d0
      posy=endy
      posz=endz-1.d0
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.(endz.gt.0.5))
     . then
      ni=ni
      nj=nj+1
      nk=nk+1
      posx=endx
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.
     . (endz.gt.0.5)) then
      ni=ni
      nj=nj-1
      nk=nk+1
      posx=endx
      posy=endy+1.d0
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.(endz.gt.0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk+1
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.(endz.gt.0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk+1
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.(endz.gt.0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk+1
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.(endz.gt.0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk+1
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz-1.d0
      elseif((abs(endx).le.0.5).and.(abs(endy).le.0.5).and.
     . (endz.lt.-0.5)) then
      ni=ni
      nj=nj
      nk=nk-1
      posx=endx
      posy=endy
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(abs(endy).le.0.5).and.
     . (endz.lt.-0.5)) then
      ni=ni+1
      nj=nj
      nk=nk-1
      posx=endx-1.d0
      posy=endy
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni-1
      nj=nj
      nk=nk-1
      posx=endx+1.d0
      posy=endy
      posz=endz+1.d0
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni
      nj=nj+1
      nk=nk-1
      posx=endx
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni
      nj=nj-1
      nk=nk-1
      posx=endx
      posy=endy+1.d0
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk-1
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk-1
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk-1
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk-1
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz+1.d0
      endif
      
      return
      end

      
      ! subroutine [endpoint]: to search along the gradient trajectories
      subroutine endpoint(posx,posy,posz,endx,endy,endz,switch,
     . subdvx,subdvy,subdvz,endfind,jumpratio)
      implicit none
      double precision jumpratio
      double precision posx,posy,posz,nposx,nposy,nposz,endx,endy,endz
      double precision switch,angle(3),temp,temp1,temp2,temp3
      double precision subdvx(-1:0,-1:1,-1:1),subdvy(-1:1,-1:0,-1:1),
     . subdvz(-1:1,-1:1,-1:0)
      double precision normx,normy,normz,normt,jump
      integer ix,iy,iz,inner,endfind,force
      double precision xp,xm,yp,ym,zp,zm
      double precision vxp(3),vxm(3),vyp(3),vym(3),vzp(3),vzm(3)
      double precision diverg,totv,vsum(3)

      include "param.h" 
      
      inner=1
      force=0

      nposx=posx
      nposy=posy      
      nposz=posz
      
      normx=0.d0
      normy=0.d0
      normz=0.d0
      do ix=-1,0
      do iy=-1,1
      do iz=-1,1
      normx=max(dabs(subdvx(ix,iy,iz)),normx)
      enddo
      enddo
      enddo
      do ix=-1,1
      do iy=-1,0
      do iz=-1,1
      normy=max(dabs(subdvy(ix,iy,iz)),normy)
      enddo
      enddo
      enddo
      do ix=-1,1
      do iy=-1,1
      do iz=-1,0
      normz=max(dabs(subdvz(ix,iy,iz)),normz)
      enddo
      enddo
      enddo
      normt=dsqrt(normx*normx+normy*normy+normz*normz)+1.d-20
      
      
      do while(inner.eq.1)      
      ! if the searching distance along any direction is beyond 0.5 
      ! or the local extremal point has been found, then exit the [endpoint] subroutine
      if((abs(nposx).gt.0.5d0).or.(abs(nposy).gt.0.5d0).or.
     . (abs(nposz).gt.0.5d0)) then
      inner=0
      endx=nposx
      endy=nposy
      endz=nposz
      goto 54
      endif

      ! subroutine [findangle]: to get the direction of the gradient
      call findangle(nposx,nposy,nposz,switch,
     . subdvx,subdvy,subdvz,angle,totv)
     
      jump=min(pace,totv/normt/4.d0)

      ! if the norm of gradient is small enough, judge whether
      ! we find the extremal point or the saddle point
      if(jump.lt.jumpeps) then
      xp=nposx+pace
      xm=nposx-pace
      yp=nposy+pace
      ym=nposy-pace
      zp=nposz+pace
      zm=nposz-pace
      call findangle(xp,nposy,nposz,switch,subdvx,subdvy,subdvz,
     .  vxp,temp)
      call findangle(xm,nposy,nposz,switch,subdvx,subdvy,subdvz,
     .  vxm,temp)
      call findangle(nposx,yp,nposz,switch,subdvx,subdvy,subdvz,
     .  vyp,temp)
      call findangle(nposx,ym,nposz,switch,subdvx,subdvy,subdvz,
     .  vym,temp)
      call findangle(nposx,nposy,zp,switch,subdvx,subdvy,subdvz,
     .  vzp,temp)
      call findangle(nposx,nposy,zm,switch,subdvx,subdvy,subdvz,
     .  vzm,temp)
       diverg=vxp(1)-vxm(1)+vyp(2)-vym(2)+vzp(3)-vzm(3)       
       if(diverg.gt.-3.5) then
       call random_number(vsum(1))
       call random_number(vsum(2))        
       call random_number(vsum(3))       
       if(abs(angle(1)).ge.abs(angle(2)).
     .  and.abs(angle(1)).ge.abs(angle(3))) then
        temp2=vsum(2)**2*angle(2)
        temp3=vsum(3)**2*angle(3)
        temp1=10*(temp2*angle(2)+temp3*angle(3))/angle(1)
        temp=1.d-20+dsqrt(temp1**2+temp2**2+temp3**2)
        nposx=nposx-temp1/temp*pace
        nposy=nposy+temp2/temp*pace
        nposz=nposz+temp3/temp*pace
        goto 55
        elseif(abs(angle(2)).ge.abs(angle(1)).
     .  and.abs(angle(2)).ge.abs(angle(3))) then
        temp1=vsum(1)**2*angle(1)
        temp3=vsum(3)**2*angle(3)
        temp2=10*(temp1*angle(1)+temp3*angle(3))/angle(2)
        temp=1.d-20+dsqrt(temp1**2+temp2**2+temp3**2)
        nposx=nposx+temp1/temp*pace
        nposy=nposy-temp2/temp*pace
        nposz=nposz+temp3/temp*pace
        goto 55
        elseif(abs(angle(3)).ge.abs(angle(1)).
     .  and.abs(angle(3)).ge.abs(angle(2))) then
        temp1=vsum(1)**2*angle(1)
        temp2=vsum(2)**2*angle(2)
        temp3=10*(temp1*angle(1)+temp2*angle(2))/angle(3)
        temp=1.d-20+dsqrt(temp1**2+temp2**2+temp3**2)
        nposx=nposx+temp1/temp*pace
        nposy=nposy+temp2/temp*pace
        nposz=nposz-temp3/temp*pace
        goto 55
        endif
       endif
       endfind=1     
       inner=0
       endx = nposx
       endy = nposy
       endz = nposz
       goto 54    
      else
      !  if the norm of gradient is not small enough, trace along the gradient trajectories
       nposx=nposx+angle(1)*jump                 
       nposy=nposy+angle(2)*jump
       nposz=nposz+angle(3)*jump
      endif

  55  force=force+1                              
      if(force.gt.10000) then
      jumpratio=jumpratio+1
      endfind=1
      inner=0
      endx=nposx
      endy=nposy
      endz=nposz
      goto 54
      endif  
      enddo
      
  54  return
      
      end

      
      ! subroutine [findangle]: to get the direction of the gradient
      subroutine findangle(nposx,nposy,nposz,switch,
     . subdvx,subdvy,subdvz,angl,totv)
      implicit none
      double precision nposx,nposy,nposz,totv
      double precision switch,angl(3)
      double precision subdvx(-1:0,-1:1,-1:1),subdvy(-1:1,-1:0,-1:1),
     . subdvz(-1:1,-1:1,-1:0)
      double precision sx,sy,sz
      integer ix,iy,iz

      ix=-1
      iy=0
      iz=0
      sx=nposx+0.5d0
      sy=nposy
      sz=nposz
      if(nposy.lt.0.d0) then
       sy=1.d0+nposy
       iy=-1
      endif
      if(nposz.lt.0.d0) then
       sz=1.d0+nposz
       iz=-1
      endif
      angl(1)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvx(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvx(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvx(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvx(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvx(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvx(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvx(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvx(ix+1,iy+1,iz+1)    ! linear interpolation
      angl(1)=angl(1)*switch

      ix=0
      iy=-1
      iz=0
      sx=nposx
      sy=nposy+0.5d0
      sz=nposz
      if(nposx.lt.0.d0) then
       sx=1.d0+nposx      
       ix=-1
      endif
      if(nposz.lt.0.d0) then
       sz=1.d0+nposz
       iz=-1
      endif
      angl(2)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvy(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvy(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvy(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvy(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvy(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvy(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvy(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvy(ix+1,iy+1,iz+1)
      angl(2)=angl(2)*switch

      ix=0
      iy=0
      iz=-1
      sx=nposx
      sy=nposy
      sz=nposz+0.5d0
      if(nposx.lt.0.d0) then
       sx=1.d0+nposx      
       ix=-1
      endif
      if(nposy.lt.0.d0) then
       sy=1.d0+nposy
       iy=-1
      endif
      angl(3)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvz(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvz(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvz(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvz(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvz(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvz(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvz(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvz(ix+1,iy+1,iz+1)
      angl(3)=angl(3)*switch
      
      totv=1.d-20+dsqrt(angl(1)**2+angl(2)**2+angl(3)**2)
      angl(1)=1.d-20+angl(1)/totv
      angl(2)=1.d-20+angl(2)/totv
      angl(3)=1.d-20+angl(3)/totv
      
      return
      end

      
      ! subroutine [numberfind]: to match the extreme point and its corresponding grid point
      subroutine numberfind(idmm,i,j,k,numends,sign,tx,ty,tz)
      implicit none
      include "param.h" 

      integer idmm,i,j,k
      double precision sign(5000000)
      double precision distance,tx,ty,tz
      integer k0,numends
      
      double precision, dimension(:,:),pointer:: points
      integer, dimension(:,:,:,:),pointer:: numbering
      
      common/large6/ numbering
      common/large5/ points     

      numbering(i,j,k,idmm)=0

      if(numends.eq.0) then
      numends=numends+1
      numbering(i,j,k,idmm)=numends
      points(numends,1)=tx
      points(numends,2)=ty
      points(numends,3)=tz
      sign(numends)=idmm
      else
      k0=1
      do while(k0.le.numends)
      if(dabs(tx-points(k0,1)).le.1) then
      if(dabs(ty-points(k0,2)).le.1) then          
      if(dabs(tz-points(k0,3)).le.1) then
      distance=dsqrt( (tx-points(k0,1))**2+(ty-points(k0,2))**2 +
     . (tz-points(k0,3))**2 )
      if(distance.lt.err) then    ! to ignore the error of the same extremal point tracing from different grids  
      numbering(i,j,k,idmm)=k0
      return
      endif
      endif
      endif
      endif
      k0=k0+1
      enddo
      numends=numends+1
      numbering(i,j,k,idmm)=numends
      points(numends,1)=tx
      points(numends,2)=ty      
      points(numends,3)=tz
      sign(numends)=idmm
      endif
      
      return            
      end

      
      ! subroutine [process]: to match the maximal and minimal point in the same DE
      subroutine process(numends,numpairs,sc) 
      implicit none
      include "param.h"
      
      character(len=50):: sc
      integer status,numends
      integer numpairs,i,j,k,k0,px,py,pz
      integer sum,locmax,locmin
      integer ni,nj
      integer pairing(1000000,3)   ! pairing(*,1) - index of maximal point in *-index DE
                                                               ! pairing(*,2) - index of minimal point in *-index DE
                                                               ! pairing(*,3) - volume of *-index DE
      integer, dimension(:,:,:),pointer:: note    ! note(I,J,K) - index of DE which grid point (I,J,K) corresponds to
      double precision pairing_ave(1000000)   ! the average value of all spatial points within each DE
      double precision dzmean,dlmean,zmaxmean,zminmean,zmidmean
      double precision zmax,zmin,dx,dy,dz
      double precision jpdf(0:100,0:100),mpdf(0:100)
      double precision jpdfw(0:100,0:100),mpdfw(0:100)
      
      double precision jpdfmin(0:100,0:100),jpdfmax(0:100,0:100)
      double precision jpdfmid(0:100,0:100)
            
      double precision, dimension(:,:,:),pointer:: psi
      double precision, dimension(:,:),pointer:: points
      integer, dimension(:,:,:,:),pointer:: numbering

      common/large0/ psi
      common/large1/ pairing
      common/large6/ numbering
      common/large5/ points     
      
      allocate(note(ndimx,ndimy,ndimz),STAT=status) 
            
      write(*,*) 'PROCESS READY!'
      numpairs = 0
      pairing = 0
      pairing_ave = 0.d0

      ! matching process
      do i=1,ndimx
      do j=1,ndimy
      do k=1,ndimz
      note(i,j,k)=0
      if(numbering(i,j,k,1).gt.0.and.numbering(i,j,k,2).gt.0.and.
     . numbering(i,j,k,1).ne.numbering(i,j,k,2)) then
      locmax=numbering(i,j,k,1)
      locmin=numbering(i,j,k,2)

      if(numpairs.eq.0) then
      numpairs=1
      pairing(1,1)=locmax
      pairing(1,2)=locmin
      pairing(1,3)=pairing(1,3)+1                   ! DE volume
      pairing_ave(1) = pairing_ave(1) + psi(i,j,k)
      note(i,j,k)=1
      goto 35
      else
      do k0=1,numpairs
      if(locmax.eq.pairing(k0,1).and.locmin.eq.pairing(k0,2)) then
      note(i,j,k)=k0
      pairing(k0,3)=pairing(k0,3)+1
      pairing_ave(k0) = pairing_ave(k0) + psi(i,j,k)
      go to 35
      endif
      enddo
      endif

      numpairs=numpairs+1
      pairing(numpairs,1)=locmax
      pairing(numpairs,2)=locmin
      pairing(numpairs,3)=pairing(numpairs,3)+1
      pairing_ave(numpairs) = pairing_ave(numpairs) + psi(i,j,k)
      note(i,j,k)=numpairs
      endif

35    enddo
      enddo  
      enddo
      
      write(*,*) 'numpairs=',numpairs
      
      open(20,file = './note'//trim(adjustl(sc))//'.dat')
      do i = 1,ndimx
      do j = 1,ndimy
      do k = 1,ndimz
          write(20,*) note(i,j,k)
      enddo
      enddo
      enddo
      close(20)
      
      open(20,file = './ave_in_DE'//Trim(AdjustL(sc))//'.dat')
      do i = 1,numpairs
          write(20,'(f18.8)') pairing_ave(i)
      enddo
      close(20)
      
      zmaxmean=0.d0
      zminmean=0.d0
      zmidmean=0.d0
      dzmean=0.d0
      dlmean=0.d0
      sum=0

      ! to get the mean value of maximal point value, minimal point value of DE
      ! and linear distance between two extremal points of DE
      do k0=1,numpairs

      locmax=pairing(k0,1)
      locmin=pairing(k0,2)
     
      px=points(locmax,1)
      py=points(locmax,2)
      pz=points(locmax,3)
      dx=points(locmax,1)-px
      dy=points(locmax,2)-py
      dz=points(locmax,3)-pz      
      zmax=(1.d0-dx)*(1.d0-dy)*(1.d0-dz)*psi(px,py,pz)+
     . dx*(1.d0-dy)*(1.d0-dz)*psi(px+1,py,pz)+(1.d0-dx)*dy*(1.d0-dz)
     . *psi(px,py+1,pz)+(1.d0-dx)*(1.d0-dy)*dz*psi(px,py,pz+1)+
     . (1.d0-dx)*dy*dz*psi(px,py+1,pz+1)+dx*(1.d0-dy)*dz
     . *psi(px+1,py,pz+1)+dx*dy*(1.d0-dz)*psi(px+1,py+1,pz)+
     . dx*dy*dz*psi(px+1,py+1,pz+1)
      
      px=points(locmin,1)
      py=points(locmin,2)
      pz=points(locmin,3)
      dx=points(locmin,1)-px
      dy=points(locmin,2)-py
      dz=points(locmin,3)-pz
      zmin=(1.d0-dx)*(1.d0-dy)*(1.d0-dz)*psi(px,py,pz)+
     . dx*(1.d0-dy)*(1.d0-dz)*psi(px+1,py,pz)+(1.d0-dx)*dy*(1.d0-dz)
     . *psi(px,py+1,pz)+(1.d0-dx)*(1.d0-dy)*dz*psi(px,py,pz+1)+
     . (1.d0-dx)*dy*dz*psi(px,py+1,pz+1)+dx*(1.d0-dy)*dz
     . *psi(px+1,py,pz+1)+dx*dy*(1.d0-dz)*psi(px+1,py+1,pz)+
     . dx*dy*dz*psi(px+1,py+1,pz+1)
           
      dx=(points(locmax,1)-points(locmin,1))*dimx/(ndimx-1)
      dy=(points(locmax,2)-points(locmin,2))*dimy/(ndimy-1)
      dz=(points(locmax,3)-points(locmin,3))*dimz/(ndimz-1)

      zmaxmean=zmaxmean+zmax
      zminmean=zminmean+zmin
      zmidmean=zmidmean+(zmax+zmin)/2.d0
      dzmean=dzmean+dabs(zmax-zmin)
      dlmean=dlmean+dsqrt(dx**2+dy**2)
      sum=sum+1
      enddo

      zmaxmean=zmaxmean/sum
      zminmean=zminmean/sum
      zmidmean=zmidmean/sum
      dzmean=dzmean/sum
      dlmean=dlmean/sum


      write(*,*) 'sum=',sum
      write(*,*) 'dzmean,dlmean,numpairs=',dzmean,dlmean,numpairs
      write(*,*) 'zmaxmean,zminmean,zmidmean=',
     . zmaxmean,zminmean,zmidmean
      
      do i=0,100
      mpdf(i)=0.d0
      mpdfw(i)=0.d0
      do j=0,100
      jpdf(i,j)=0.d0
      jpdfw(i,j)=0.d0
      jpdfmin(i,j)=0.d0
      jpdfmax(i,j)=0.d0
      jpdfmid(i,j)=0.d0
      enddo
      enddo
      
      ! get the joint PDF between the linear distance and 
      ! a. scalar difference between two extremal points
      ! b. maxima
      ! c. minima
      ! d. scalar average between two extremal points of DEs
      do k0=1,numpairs
      
      locmax=pairing(k0,1)
      locmin=pairing(k0,2)

      px=points(locmax,1)
      py=points(locmax,2)
      pz=points(locmax,3)
      dx=points(locmax,1)-px
      dy=points(locmax,2)-py
      dz=points(locmax,3)-pz
      zmax=(1.d0-dx)*(1.d0-dy)*(1.d0-dz)*psi(px,py,pz)+
     . dx*(1.d0-dy)*(1.d0-dz)*psi(px+1,py,pz)+(1.d0-dx)*dy*(1.d0-dz)
     . *psi(px,py+1,pz)+(1.d0-dx)*(1.d0-dy)*dz*psi(px,py,pz+1)+
     . (1.d0-dx)*dy*dz*psi(px,py+1,pz+1)+dx*(1.d0-dy)*dz
     . *psi(px+1,py,pz+1)+dx*dy*(1.d0-dz)*psi(px+1,py+1,pz)+
     . dx*dy*dz*psi(px+1,py+1,pz+1)
      
      px=points(locmin,1)
      py=points(locmin,2)
      pz=points(locmin,3)
      dx=points(locmin,1)-px
      dy=points(locmin,2)-py
      dz=points(locmin,3)-pz
      zmin=(1.d0-dx)*(1.d0-dy)*(1.d0-dz)*psi(px,py,pz)+
     . dx*(1.d0-dy)*(1.d0-dz)*psi(px+1,py,pz)+(1.d0-dx)*dy*(1.d0-dz)
     . *psi(px,py+1,pz)+(1.d0-dx)*(1.d0-dy)*dz*psi(px,py,pz+1)+
     . (1.d0-dx)*dy*dz*psi(px,py+1,pz+1)+dx*(1.d0-dy)*dz
     . *psi(px+1,py,pz+1)+dx*dy*(1.d0-dz)*psi(px+1,py+1,pz)+
     . dx*dy*dz*psi(px+1,py+1,pz+1)
     
      dx=(points(locmax,1)-points(locmin,1))*dimx/(ndimx-1)
      dy=(points(locmax,2)-points(locmin,2))*dimy/(ndimy-1)
      dz=(points(locmax,3)-points(locmin,3))*dimz/(ndimz-1)
      
      ni=25*dsqrt(dx**2+dy**2+dz**2)/dlmean
      ni=min(ni,100)
      ni=max(ni,0)
      nj=30*(zmax-zmin)/dzmean
      nj=min(nj,100)
      nj=max(nj,0)
      jpdf(ni,nj)=jpdf(ni,nj)+1.d0
      jpdfw(ni,nj)=jpdfw(ni,nj)+pairing(k0,3)
      
      nj=30*zmax/zmaxmean
      nj=min(nj,100)
      nj=max(nj,0)
      jpdfmax(ni,nj)=jpdfmax(ni,nj)+1.d0
      
      nj=30*zmin/zminmean
      nj=min(nj,100)
      nj=max(nj,0)
      jpdfmin(ni,nj)=jpdfmin(ni,nj)+1.d0

      nj=30*(zmax+zmin)/2/zmidmean
      nj=min(nj,100)
      nj=max(nj,0)
      jpdfmid(ni,nj)=jpdfmid(ni,nj)+1.d0
      enddo
      
      do i=0,100
      do j=0,100
      mpdf(i)=mpdf(i)+jpdf(i,j)
      mpdfw(i)=mpdfw(i)+jpdfw(i,j)
      enddo
      enddo
      
55    format(2(es14.4,1x),es14.4)
      write(*,*) 'PROCESS END!'  
      end
