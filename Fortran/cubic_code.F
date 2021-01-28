c this code is used to test the results from 512^3
c here we want to study the average shape of elements
c from scalar k and eps.

      program main
      implicit none
      include "param.h"

      character(len=50) sc
            
      integer i,j,k,ni,nj,nk,endfind,i0,j0,k0,numpairs,c,m
      integer numends,findex,loop,force0,status
      integer pairing(1000000,3)
      integer, dimension(:,:,:),pointer:: numbering
      
      double precision jumpratio,tmain0,tmain1,tdiff0,tdiff1
      double precision reed,st,rtime,dzmean,dlmean,vis
      double precision sign(5000000)
      double precision switch,posx,posy,posz,endx,endy,endz
      double precision subdvx(-1:0,-1:1),subdvy(-1:1,-1:0)
      double precision, dimension(:,:),pointer:: points
      double precision, dimension(:,:,:), pointer:: u0
      double precision, dimension(:,:),pointer:: psi
      double precision, dimension(:,:),pointer:: coordx,coordy
      
      common/large0/ psi
      common/large1/ pairing
      common/large2/ coordx,coordy
      common/large6/ numbering
      common/large5/ points      
      
      allocate(psi(ndimx,ndimy),STAT=status)
      allocate(points(500000,2),STAT=status)
      allocate(numbering(ndimx,ndimy,2),STAT=status)
      allocate(coordx(ndimx,ndimy),STAT=status)
      allocate(coordy(ndimx,ndimy),STAT=status)
      
      do m = 1,2922!1,240!do m = 1,31
      write(sc,*) m
      write(*,*) 'we are begining snapshot '//Trim(AdjustL(sc))
      write(*,*) '======= begin reading ======='
!      open(20,file='../../../data2010-2017/
!     .ocean_summer'//Trim(AdjustL(sc))//'.dat')
	  !!!======= data input =======
      open(20,file='../../../data_allover/
     .ocean_summer'//Trim(AdjustL(sc))//'.dat')
      read(20,*)
      read(20,*)
      do j = 1,301
          do i = 1,229
              read(20,*) reed,reed,psi(i,j)
          enddo
      enddo
      close(20)
      write(*,*) '======= end reading ======='
      
      call CPU_TIME(tmain0)   
      numends = 0
      numbering = 0
      sign = 0.d0
      points = 0.d0
      jumpratio=0.d0
      dzmean=0.d0
      dlmean=0.d0
      loop=0
      write(*,*) '======= begin coord calculation ======='
      do j=1,ndimy
      do i=1,ndimx
      coordx(i,j)=(i-1)*2*pi/(ndimx-1)
      coordy(i,j)=(j-1)*2*pi/(ndimy-1)
      enddo
      enddo
      write(*,*) '======= end coord calculation ======='
          
      !!! algorithm to find local extremal points and 0+ level DE
      do i=2,ndimx-1
       !write(*,*) 'i=',i,numends
       do j=2,ndimy-1             !from 2 to ndimx-1, for 1 and ndimx are on the boundary         
      if (psi(i,j).ge.1.d0) then
       force0=0
       endfind=0
       posx=0      !distance between searching point and grid point,ranging 0~1
       posy=0
       ni=i
       nj=j           

c switch=1 means we do calculation in positive direction
       switch=1.d0
       do while(endfind.eq.0)       

       if((ni.eq.1).or.(ni.eq.ndimx).or.(nj.eq.1).or.(nj.eq.ndimy)) then
       goto 25
       endif
       
       do i0=-1,0
       do j0=-1,1
       subdvx(i0,j0)= (psi(ni+i0+1,nj+j0)-psi(ni+i0,nj+j0))/(2*pi/ndimx)
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       subdvy(i0,j0)= (psi(ni+i0,nj+j0+1)-psi(ni+i0,nj+j0))/(2*pi/ndimy)
       enddo
       enddo
       
       call endpoint(posx,posy,endx,endy,switch,subdvx,subdvy,endfind,
     . jumpratio)
	 
       force0=force0+1                  !force0 for local vortex situation
       if(force0.gt.500) then
       endfind=1
       goto 25       
       endif

       if(endfind.eq.1) then
       call numberfind(1,i,j,numends,sign,ni+endx,nj+endy)
       endif
       
       call tellnew(ni,nj,endx,endy,posx,posy)

       enddo
       
   25  force0=0
       endfind=0
       posx=0
       ni=i
       nj=j

c switch=-1 means we do calculation in negative direction
       switch=-1.d0
       do while(endfind.eq.0)

       if((ni.eq.1).or.(ni.eq.ndimx).or.(nj.eq.1).or.(nj.eq.ndimy)) then
       goto 30
       endif

       do i0=-1,0
       do j0=-1,1
       subdvx(i0,j0)= (psi(ni+i0+1,nj+j0)-psi(ni+i0,nj+j0))/(2*pi/ndimx)
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       subdvy(i0,j0)= (psi(ni+i0,nj+j0+1)-psi(ni+i0,nj+j0))/(2*pi/ndimy)
       enddo
       enddo
       
       call endpoint(posx,posy,endx,endy,switch,subdvx,subdvy,endfind,
     . jumpratio)

       force0=force0+1
       if(force0.gt.500) then
       endfind=1
       goto 30
       endif
       
       if(endfind.eq.1) then
       call numberfind(2,i,j,numends,sign,ni+endx,nj+endy)
       endif

       call tellnew(ni,nj,endx,endy,posx,posy)

       enddo
      
       endif
   30  enddo
       enddo
          
      write(*,*) 'jumpration=',jumpratio/ndimx/ndimy
      write(*,*) 'numends=',numends     
      call process(numends,numpairs,sc)
      write(*,*) 'numpairs=',numpairs

      call CPU_TIME(tmain1)
      write(*,*) 'main duration=',tmain1-tmain0
        
	  !!! output the info. of 0+ DE for multi-level algorithm code
      !open(12,file='./datas/numends'//Trim(AdjustL(sc))//'.dat')
      open(12,file='./datas_overall/numends'//Trim(AdjustL(sc))//'.dat')
      write(12,*) numends
      close(12)
      
      !open(12,file='./datas/numpairs'//Trim(AdjustL(sc))//'.dat')
      open(12,file='./datas_overall/numpairs'//Trim(AdjustL(sc))//
     .'.dat')
      write(12,*) numpairs
      close(12)   
      
      !open(12,file='./datas/pairing'//Trim(AdjustL(sc))//'.dat')
      open(12,file='./datas_overall/pairing'//Trim(AdjustL(sc))//'.dat')
      do i=1,numpairs
          write(12,*) pairing(i,1),pairing(i,2),pairing(i,3)
      enddo
      close(12)   
      
      !open(12,file='./datas/points'//Trim(AdjustL(sc))//'.dat')
      open(12,file='./datas_overall/points'//Trim(AdjustL(sc))//'.dat')
      write(12,*) ' variables = "x","y","sign" '
      do i=1,numends          
          write(12,'(3f12.6)') points(i,1),points(i,2),sign(i)
      enddo
      close(12)
      
      enddo
      
      end

      subroutine tellnew(ni,nj,endx,endy,posx,posy)   ! to judge which grid point the searching point correspond to 
      implicit none
      integer ni,nj
      double precision endx,endy,posx,posy           
      !   end(x y z) relative coordinate ranging -1~1;pos(x y z)simpliest coordinate ranging -0.5~0.5.
      !   e.g:for real(1.7,2.3,2.4): end(0.7,0.3,0.4)--(1,2,2)   pos(-0.3,0.3,0.4)--(2,2,2)

      if(abs(endx).le.0.5.and.abs(endy).le.0.5) then
      ni=ni
      nj=nj
      posx=endx
      posy=endy
      elseif(endx.gt.0.5.and.abs(endy).le.0.5) then
      ni=ni+1
      nj=nj
      posx=endx-1.d0            
      posy=endy
      elseif(endx.lt.-0.5.and.abs(endy).le.0.5) then
      ni=ni-1
      nj=nj
      posx=endx+1.d0
      posy=endy
      elseif(abs(endx).le.0.5.and.endy.gt.0.5) then
      ni=ni
      nj=nj+1
      posx=endx
      posy=endy-1.d0
      elseif(abs(endx).le.0.5.and.endy.lt.-0.5) then
      ni=ni
      nj=nj-1
      posx=endx
      posy=endy+1.d0
      elseif(endx.gt.0.5.and.endy.gt.0.5) then
      ni=ni+1
      nj=nj+1
      posx=endx-1.d0
      posy=endy-1.d0
      elseif(endx.lt.-0.5.and.endy.gt.0.5) then
      ni=ni-1
      nj=nj+1
      posx=endx+1.d0
      posy=endy-1.d0
      elseif(endx.lt.-0.5.and.endy.lt.-0.5) then
      ni=ni-1
      nj=nj-1
      posx=endx+1.d0
      posy=endy+1.d0
      elseif(endx.gt.0.5.and.endy.lt.-0.5) then
      ni=ni+1
      nj=nj-1
      posx=endx-1.d0
      posy=endy+1.d0
      endif
      
      return
      end

      subroutine endpoint(posx,posy,endx,endy,switch,
     .subdvx,subdvy,endfind,jumpratio)    ! search extrema along gradient trajectory
      implicit none
      double precision jumpratio
      double precision posx,posy,nposx,nposy,endx,endy
      double precision switch,angle(3),temp,temp1,temp2,temp3
      double precision subdvx(-1:0,-1:1),subdvy(-1:1,-1:0)
      double precision normx,normy,normt,jump
      integer ix,iy,iz,inner,endfind,force
      double precision xp,xm,yp,ym,zp,zm
      double precision vxp(3),vxm(3),vyp(3),vym(3)
      double precision diverg,totv,vsum(3)

      include "param.h" 
      
      inner=1
      force=0

      nposx=posx
      nposy=posy
      
      normx=0.d0
      normy=0.d0
      do ix=-1,0
      do iy=-1,1
      normx=max(dabs(subdvx(ix,iy)),normx)
      enddo
      enddo
      do ix=-1,1
      do iy=-1,0
      normy=max(dabs(subdvy(ix,iy)),normy)
      enddo
      enddo      
      normt=dsqrt(normx*normx+normy*normy)+1.d-20
      
      do while(inner.eq.1)
       
      if(abs(nposx).gt.0.5d0.or.abs(nposy).gt.0.5d0) then
      inner=0
      endx=nposx
      endy=nposy
      goto 54
      endif

      call findangle(nposx,nposy,switch,subdvx,subdvy,angle,totv)
     
      jump=min(pace,totv/normt/4.d0)

      if(jump.lt.jumpeps) then

       xp=nposx+pace
       xm=nposx-pace
       yp=nposy+pace
       ym=nposy-pace
       call findangle(xp,nposy,switch,subdvx,subdvy,vxp,temp)
       call findangle(xm,nposy,switch,subdvx,subdvy,vxm,temp)
       call findangle(nposx,yp,switch,subdvx,subdvy,vyp,temp)
       call findangle(nposx,ym,switch,subdvx,subdvy,vym,temp)
       diverg=vxp(1)-vxm(1)+vyp(2)-vym(2)
       if(diverg.gt.-2.5) then
       call random_number(vsum(1))
       call random_number(vsum(2)) 
       
       if(abs(angle(1)).ge.abs(angle(2))) then
       temp2=vsum(2)**2*angle(2)
        temp1=10*temp2*angle(2)/angle(1)
        temp=1.d-20+dsqrt(temp1**2+temp2**2)
        nposx=nposx-temp1/temp*pace
        nposy=nposy+temp2/temp*pace
        goto 55
        else
        temp1=vsum(1)**2*angle(1)
        temp2=10*temp1*angle(1)/angle(2)
        temp=1.d-20+dsqrt(temp1**2+temp2**2)
        nposx=nposx+temp1/temp*pace
        nposy=nposy-temp2/temp*pace
        goto 55
        endif
       endif

       endfind=1     
       inner=0
       endx=nposx
       endy=nposy
       goto 54
      
      else
       nposx=nposx+angle(1)*jump                 
       nposy=nposy+angle(2)*jump
      endif

  55  force=force+1                              
      if(force.gt.10000) then
      jumpratio=jumpratio+1
      endfind=1
      inner=0
      endx=nposx
      endy=nposy
      goto 54
      endif
  
      enddo
      
  54  return
      
      end

      subroutine findangle(nposx,nposy,switch,subdvx,subdvy,angl,totv) ! get the angle of gradient vector
       implicit none
      double precision nposx,nposy,totv
      double precision switch,angl(2)
      double precision subdvx(-1:0,-1:1),subdvy(-1:1,-1:0)
      double precision sx,sy
      integer ix,iy

      ix=-1
      iy=0
      sx=nposx+0.5d0
      sy=nposy
      if(nposy.lt.0.d0) then
       sy=1.d0+nposy
       iy=-1
      endif
      angl(1)=(1.d0-sx)*(1.d0-sy)*subdvx(ix,iy)+  ! linear interpolation
     . sx*(1.d0-sy)*subdvx(ix+1,iy)+
     . (1.d0-sx)*sy*subdvx(ix,iy+1)+
     . sx*sy*subdvx(ix+1,iy+1)
      angl(1)=angl(1)*switch

      ix=0
      iy=-1
      sx=nposx
      sy=nposy+0.5d0
      if(nposx.lt.0.d0) then
       sx=1.d0+nposx      
       ix=-1
      endif
      angl(2)=(1.d0-sx)*(1.d0-sy)*subdvy(ix,iy)+
     . sx*(1.d0-sy)*subdvy(ix+1,iy)+
     . (1.d0-sx)*sy*subdvy(ix,iy+1)+
     . sx*sy*subdvy(ix+1,iy+1)
      angl(2)=angl(2)*switch
      
      totv=1.d-20+dsqrt(angl(1)**2+angl(2)**2)
      angl(1)=1.d-20+angl(1)/totv
      angl(2)=1.d-20+angl(2)/totv
      
      return
      end

      subroutine numberfind(idmm,i,j,numends,sign,tx,ty)
      implicit none
      include "param.h" 

      integer idmm,i,j
      double precision sign(5000000)
      double precision distance,tx,ty
      integer k0,numends
      
      double precision, dimension(:,:),pointer:: points
      integer, dimension(:,:,:),pointer:: numbering
      
      common/large6/ numbering
      common/large5/ points     

! meaning of numbering(,,,2) is: 1.num of max point 2.num of min point

      numbering(i,j,idmm)=0

      if(numends.eq.0) then
      numends=numends+1
      numbering(i,j,idmm)=numends
      points(numends,1)=tx
      points(numends,2)=ty
      sign(numends)=idmm
      else
      k0=1
      do while(k0.le.numends)
      if(dabs(tx-points(k0,1)).le.1) then
      if(dabs(ty-points(k0,2)).le.1) then

      distance=dsqrt((tx-points(k0,1))**2+(ty-points(k0,2))**2)
      if(distance.lt.err) then        !to ignore the error of the same extrema from different grids
      numbering(i,j,idmm)=k0
      return  ! go back to program main
      endif

      endif
      endif
      k0=k0+1
      enddo
      numends=numends+1
      numbering(i,j,idmm)=numends
      points(numends,1)=tx
      points(numends,2)=ty
      sign(numends)=idmm
      endif
      
      return            
      end

      subroutine process(numends,numpairs,sc) 
      implicit none
      include "param.h"
      
      character(len=50):: sc
      integer status,numends
      integer numpairs,i,j,k0,px,py
      integer sum,locmax,locmin
      integer ni,nj
      integer pairing(1000000,3)
      integer, dimension(:,:),pointer:: note
      double precision pairing_tao(1000000)
      double precision dzmean,dlmean,zmaxmean,zminmean,zmidmean
      double precision zmax,zmin,dx,dy,dz
      double precision jpdf(0:100,0:100),mpdf(0:100)
      double precision jpdfw(0:100,0:100),mpdfw(0:100)
      
      double precision jpdfmin(0:100,0:100),jpdfmax(0:100,0:100)
      double precision jpdfmid(0:100,0:100)
      double precision rad(18,2,6),span(18,2,6),newi,newj,vol,tmp
      double precision cxmean,cymean
      integer pos0
            
      double precision, dimension(:,:),pointer:: psi
      double precision, dimension(:,:),pointer:: points
      integer, dimension(:,:,:),pointer:: numbering

      common/large0/ psi
      common/large1/ pairing
      common/large6/ numbering
      common/large5/ points     
      
      allocate(note(ndimx,ndimy),STAT=status) 
            
      write(*,*) 'PROCESS READY!'
      numpairs = 0
      pairing = 0
      pairing_tao = 0.d0

      do i=1,ndimx
      do j=1,ndimy
      note(i,j)=0
      if(numbering(i,j,1).gt.0.and.numbering(i,j,2).gt.0.and.
     . numbering(i,j,1).ne.numbering(i,j,2)) then
      locmax=numbering(i,j,1)
      locmin=numbering(i,j,2)

      if(numpairs.eq.0) then
      numpairs=1
      pairing(1,1)=locmax
      pairing(1,2)=locmin
      pairing(1,3)=pairing(1,3)+1                   ! DE volume
      pairing_tao(1) = pairing_tao(1) + psi(i,j)
      note(i,j)=1
      goto 35
      else
      do k0=1,numpairs
      if(locmax.eq.pairing(k0,1).and.locmin.eq.pairing(k0,2)) then
      note(i,j)=k0
      pairing(k0,3)=pairing(k0,3)+1
      pairing_tao(k0) = pairing_tao(k0) + psi(i,j)
      go to 35
      endif
      enddo
      endif

      numpairs=numpairs+1
      pairing(numpairs,1)=locmax
      pairing(numpairs,2)=locmin
      pairing(numpairs,3)=pairing(numpairs,3)+1
      pairing_tao(numpairs) = pairing_tao(numpairs) + psi(i,j)
      note(i,j)=numpairs

      endif

   35 enddo
      enddo  
      
      write(*,*) 'numpairs=',numpairs
      
      !open(20,file = 'note'//trim(adjustl(sc))//'.dat')
      !do i = 1,ndimx
      !do j = 1,ndimy
      !    write(20,*) note(i,j)
      !enddo
      !enddo
      !close(20)
      
      !open(20,file = './datas/tao_in_DE'//Trim(AdjustL(sc))//'.dat')
      open(20,file = './datas_overall/tao_in_DE'//Trim(AdjustL(sc))//
     .'.dat')
      do i = 1,numpairs
          write(20,'(f18.8)') pairing_tao(i)
      enddo
      close(20)
      
      zmaxmean=0.d0
      zminmean=0.d0
      zmidmean=0.d0
      dzmean=0.d0
      dlmean=0.d0
      sum=0

      do k0=1,numpairs

      locmax=pairing(k0,1)
      locmin=pairing(k0,2)
     
      px=points(locmax,1)
      py=points(locmax,2)
      dx=points(locmax,1)-px
      dy=points(locmax,2)-py
      !write(*,*) px,py,pz
      zmax=(1.d0-dx)*(1.d0-dy)*psi(px,py)+dx*(1.d0-dy)*psi(px+1,py)
     .+(1.d0-dx)*dy*psi(px,py+1)+dx*dy*psi(px+1,py+1)
      
      px=points(locmin,1)
      py=points(locmin,2)
      dx=points(locmin,1)-px
      dy=points(locmin,2)-py
      zmin=(1.d0-dx)*(1.d0-dy)*psi(px,py)+dx*(1.d0-dy)*psi(px+1,py)
     . +(1.d0-dx)*dy*psi(px,py+1)+dx*dy*psi(px+1,py+1)
           
      dx=(points(locmax,1)-points(locmin,1))*2*pi/(ndimx-1)
      dy=(points(locmax,2)-points(locmin,2))*2*pi/(ndimy-1)

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
      
      do k0=1,numpairs
      
      locmax=pairing(k0,1)
      locmin=pairing(k0,2)

      px=points(locmax,1)
      py=points(locmax,2)
      dx=points(locmax,1)-px
      dy=points(locmax,2)-py
      zmax=(1.d0-dx)*(1.d0-dy)*psi(px,py)+dx*(1.d0-dy)*psi(px+1,py)
     .+(1.d0-dx)*dy*psi(px,py+1)+dx*dy*psi(px+1,py+1)
      
      px=points(locmin,1)
      py=points(locmin,2)
      dx=points(locmin,1)-px
      dy=points(locmin,2)-py
      zmin=(1.d0-dx)*(1.d0-dy)*psi(px,py)+dx*(1.d0-dy)*psi(px+1,py)
     . +(1.d0-dx)*dy*psi(px,py+1)+dx*dy*psi(px+1,py+1)
     
      dx=(points(locmax,1)-points(locmin,1))*2*pi/(ndimx-1)
      dy=(points(locmax,2)-points(locmin,2))*2*pi/(ndimy-1)
      
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
      
   55 format(2(es14.4,1x),es14.4)
      write(*,*) 'PROCESS END!'  
      end
