       program main
       implicit none
       include "param.h"
       double precision,external::interp
      
       character(len=50):: R_str,sc
       integer i,j,windex,numends,zx,zy,zz,endfind,ni,nj,i0,j0,k0,kk,c
       integer windowfind,force0,status,meetbound,ii,jj,n1,n2,n3,n,m,Rball
       integer locmax,locmin,numpairs,pairsize,count,R,note(ndimx,ndimy)
       integer endid,backindex,oldmax,oldmin,numpairsm(winnum)
       integer pairing(1000000,3),k1,notem(winnum,1000000)

       double precision subdvx(-1:0,-1:1),subdvy(-1:1,-1:0),reed
       double precision maxxcor,maxycor,maxpsi,minxcor,minycor,minpsi,lenmean,dpsimean,spsimean,sum,len2mean
       double precision jpdf(0:100,0:200),tdiff0,tdiff1,mpdf(0:100),ypdf(0:100),my(0:100)
       double precision jpdf2(0:100,0:200),mpdf2(0:100),ypdf2(0:100),my2(0:100)
       double precision switch,endx,endy,jumpratio,window
       double precision xmax,xmin,ymax,ymin,px,py,x(1000000),y(1000000)
       double precision posx,posy,psimax,epsilon,inter(1000000)
       double precision minxco,minyco,mingrid,len1,deltapsi,len2
       double precision xcorball,ycorball,psiball
       double precision pairmsize(winnum,200000)
       double precision R2_mpdf(0:100),R2_my(0:100),R2_mpdf2(0:100),R2_my2(0:100)
       double precision R4_mpdf(0:100),R4_my(0:100),R4_mpdf2(0:100),R4_my2(0:100)
       double precision R6_mpdf(0:100),R6_my(0:100),R6_mpdf2(0:100),R6_my2(0:100)
       double precision R8_mpdf(0:100),R8_my(0:100),R8_mpdf2(0:100),R8_my2(0:100)
       double precision R10_mpdf(0:100),R10_my(0:100),R10_mpdf2(0:100),R10_my2(0:100)       
       double precision R12_mpdf(0:100),R12_my(0:100),R12_mpdf2(0:100),R12_my2(0:100)
       double precision R14_mpdf(0:100),R14_my(0:100),R14_mpdf2(0:100),R14_my2(0:100)
       double precision Rall_mpdf(0:100),Rall_my(0:100),Rall_mpdf2(0:100),Rall_my2(0:100)
       double precision Rall_jpdf(0:100,0:100),Rall_jpdf2(0:100,0:100)
       double precision R1_jpdf(0:100,0:100),R1_jpdf2(0:100,0:100)
       double precision R2_jpdf(0:100,0:100),R2_jpdf2(0:100,0:100)
       double precision R3_jpdf(0:100,0:100),R3_jpdf2(0:100,0:100)
       double precision, dimension(:,:),pointer:: points
       double precision, dimension(:,:),pointer:: psi
       double precision, dimension(:,:,:), pointer:: u0
       double precision, dimension(:),pointer:: dx,dy
       double precision, dimension(:,:),pointer:: coordx,coordy
       double precision, dimension(:,:,:),pointer:: multiscale
       double precision sign(5000000),pairingm(winnum,10000,2,2)
       double precision pairing_tao(1000000)
       character(len = 512)::cFile,ccFile,cccFile

      common/large0/ psi
      common/large1/ pairing
      common/large2/ coordx,coordy
      common/large5/ points

      allocate(psi(ndimx,ndimy),STAT = status)  
      allocate(points(5000000,3),STAT = status)
      allocate(coordx(ndimx,ndimy),STAT = status)
      allocate(coordy(ndimx,ndimy),STAT = status)
      allocate(dx(ndimx-1),STAT = status)
      allocate(dy(ndimy-1),STAT = status)

      
      
      R2_my = 0.d0
      R2_my2 = 0.d0
      R4_my = 0.d0
      R4_my2 = 0.d0
      R6_my = 0.d0
      R6_my2 = 0.d0
      R8_my = 0.d0
      R8_my2 = 0.d0
      R10_my = 0.d0
      R10_my2 = 0.d0
      R12_my = 0.d0
      R12_my2 = 0.d0
      R14_my = 0.d0
      R14_my2 = 0.d0
      Rall_my = 0.d0
      Rall_my2 = 0.d0
      Rall_jpdf = 0.d0
      Rall_jpdf2 = 0.d0
      R1_jpdf = 0.d0
      R1_jpdf2 = 0.d0
      R2_jpdf = 0.d0
      R2_jpdf2 = 0.d0
      R3_jpdf = 0.d0
      R3_jpdf2 = 0.d0
      
	  ! multi-level DE algorithm
      do m = 1,240!do m = 1,1
        !R2_my = 0.d0
        !R2_my2 = 0.d0
        !R4_my = 0.d0
        !R4_my2 = 0.d0
        !R6_my = 0.d0
        !R6_my2 = 0.d0
        !R8_my = 0.d0
        !R8_my2 = 0.d0
        !R10_my = 0.d0
        !R10_my2 = 0.d0
        !R12_my = 0.d0
        !R12_my2 = 0.d0
        !R14_my = 0.d0
        !R14_my2 = 0.d0
        !Rall_my = 0.d0
        !Rall_my2 = 0.d0
        !Rall_jpdf = 0.d0
        !Rall_jpdf2 = 0.d0
        !R1_jpdf = 0.d0
        !R1_jpdf2 = 0.d0
        !R2_jpdf = 0.d0
        !R2_jpdf2 = 0.d0
        !R3_jpdf = 0.d0
        !R3_jpdf2 = 0.d0
      write(sc,*) m
      write(*,*) 'we are begining snapshot '//Trim(adjustl(sc))
      write(*,*) '======= begin reading ======='
      open(20,file='../../../data2010-2017/ocean_summer'//TRim(adjustl(sc))//'.dat')
      !open(20,file='../../../data_allover/ocean_summer'//TRim(adjustl(sc))//'.dat')
      read(20,*)
      read(20,*)
      do j = 1,301
          do i = 1,229
              read(20,*) reed,reed,psi(i,j)
          enddo
      enddo
      close(20)
      write(*,*) '======= end reading ======='
      
      write(*,*) '======= begin coord calculation ======='
      do j = 1,ndimy
      do i = 1,ndimx
          coordx(i,j) = (i-1)*2*pi/(ndimx-1)
          coordy(i,j) = (j-1)*2*300/228*pi/(ndimy-1)
      enddo
      enddo
      write(*,*) '======= end coord calculation ======='
      
      numends = 0
      numpairs = 0
      sign = 0.d0
      pairing = 0
      points = 0.d0
      
      open(12,file = './datas/numends'//TRim(adjustl(sc))//'.dat')
      !open(12,file = './datas_overall/numends'//TRim(adjustl(sc))//'.dat')
      read(12,*) numends
      close(12)
      
      open(12,file = './datas/numpairs'//TRim(adjustl(sc))//'.dat')
      !open(12,file = './datas_overall/numpairs'//TRim(adjustl(sc))//'.dat')
      read(12,*) numpairs
      close(12)
      
      open(12,file = './datas/pairing'//TRim(adjustl(sc))//'.dat')
      !open(12,file = './datas_overall/pairing'//TRim(adjustl(sc))//'.dat')
      do i = 1,numpairs
          read(12,*) pairing(i,1),pairing(i,2),pairing(i,3)
      enddo
      close(12)   
      
      open(12,file = './datas/points'//TRim(adjustl(sc))//'.dat')
      !open(12,file = './datas_overall/points'//TRim(adjustl(sc))//'.dat')
      read(12,*) 
      do i = 1,numends
          read(12,*) points(i,1),points(i,2),sign(i)
      enddo
      close(12)
      
      open(12,file = './datas/tao_in_DE'//TRim(adjustl(sc))//'.dat')
      !open(12,file = './datas_overall/tao_in_DE'//TRim(adjustl(sc))//'.dat')
      do i = 1,numpairs
          read(12,*) pairing_tao(i)
      enddo
      close(12)
      
      !open(12,file = './datas/note'//TRim(adjustl(sc))//'.dat')
      !do i = 1,ndimx
      !do j = 1,ndimy
      !    read(12,*) note(i,j)
      !enddo
      !enddo
      !close(12)
      
      
      allocate(multiscale(winnum,numends,3),STAT = status)      
      
      multiscale = 0.d0
      numpairsm = 0.d0
      pairingm = 0.d0
      notem = 0.d0
      dx = 2*pi/(ndimx-1)
      dy = 2*pi/(ndimy-1)
      minxco = minval(dx)
      minyco = minval(dy)
      mingrid = min(minxco,minyco)
      epsilon = 1.d-4
      
      write(*,*)
      write(*,*) '======= begin multi-level calculation ======='
      do windex = 1,winnum    !!!!!!! high level searching process
          write(*,*) 'windex = ',windex
          do i = 1,numends
              if (windex.eq.1) then
                  multiscale(windex,i,1) = points(i,1)
                  multiscale(windex,i,2) = points(i,2)
                  multiscale(windex,i,3) = interp(points(i,1),points(i,2),psi)
              else
                  endid = 0        
                  windowfind = 0
                  !!!!!!!judge the extremal points is max or min
                  if (sign(i).eq.1) then       
                      switch = 1.d0
                  else
                      switch = -1.d0
                  endif
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ni = multiscale(windex-1,i,1)
                  nj = multiscale(windex-1,i,2)
                  endx = multiscale(windex-1,i,1)-ni
                  endy = multiscale(windex-1,i,2)-nj
                  zx = multiscale(windex-1,i,1)
                  zy = multiscale(windex-1,i,2)
                  if (zx.lt.0.or.zy.lt.0) then
                      multiscale(windex,i,:) = -100000.d0
                      goto 40
                  endif
                  if (switch.gt.0) then
                      psimax = maxval(psi(zx:zx+1,zy:zy+1))+epsilon
                  else
                      psimax = minval(psi(zx:zx+1,zy:zy+1))-epsilon
                  endif

                  Rball = 1 ! ======= searching R initialization =======
                  do while(.true.)
                      window = Rball*mingrid
                      call sphere(ni,nj,endx,endy,psimax,Rball,window,switch,psiball,xcorball,ycorball,meetbound)
60                  if(meetbound.eq.1) then
                          multiscale(windex,i,:) = -100000.d0
                          exit
                      endif
                      ! ======= If no jump, R increases =======
                      if ((meetbound.eq.0).and.(abs(psiball-psimax).le.1.d-5)) then
                          if (Rball.lt.windex) then     ! =======[If now R < windex, R++]=======
                              Rball = Rball + 1
                          else if (Rball.eq.windex) then    ! =======[If now R = windex, windowfind = 1]=======
                              windowfind = 1
                          endif
                      ! ======= If jump, then search along grad_traj =======
                      else if ((meetbound.eq.0).and.(abs(psiball-psimax).gt.1.d-5)) then
                          ! ======= search along gradient trajectory from grid/sphere-surface point
                          ni = xcorball
                          nj = ycorball
                          psimax = psiball
                          posx = xcorball - ni
                          posy = ycorball - nj
                          force0 = 0
                          endfind = 0
                          do while(endfind.eq.0)
                              if((ni.eq.1).or.(ni.eq.ndimx).or.(nj.eq.1).or.(nj.eq.ndimy)) then
                                  meetbound = 1
                                  go to 60
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
                              call endpoint(posx,posy,endx,endy,switch,subdvx,subdvy,endfind,jumpratio)                       
                              call tellnew(ni,nj,endx,endy,posx,posy)
                              force0=force0+1
                              if(force0.gt.500) then
                                  endfind=1
                                  exit
                              endif
                          enddo
                          Rball = 1
                      endif                      
                      if (windowfind.eq.1) then
                          multiscale(windex,i,1) = xcorball
                          multiscale(windex,i,2) = ycorball
                          multiscale(windex,i,3) = psimax
                          exit
                      endif
                  enddo      !!!!!!! loop when windowfind = 0
              endif      !!!!!!! loop for if windex = 1
40        enddo    !!!!!!! loop for i(numendsM)
      enddo  !!!!!!! loop for windex(winnum) high level searching process end
      write(*,*) '======= end multi-level calculation ======='
      write(*,*)
                
      write(*,*) '======= begin multi-level pairing ======='
      do windex = 1,winnum   !!!!!!! high level pairing process
       numpairsm(windex) = 0
       do k0 = 1,numpairs
           notem(windex,k0) = 0
           if(multiscale(windex,pairing(k0,1),1).gt.0.d0.and.multiscale(windex,pairing(k0,1),2).gt.0.d0.and. &
&           multiscale(windex,pairing(k0,2),1).gt.0.d0.and.multiscale(windex,pairing(k0,2),2).gt.0.d0 ) then          
           xmax = multiscale(windex,pairing(k0,1),1)
           ymax = multiscale(windex,pairing(k0,1),2)
           xmin = multiscale(windex,pairing(k0,2),1)
           ymin = multiscale(windex,pairing(k0,2),2)
           if (numpairsm(windex).eq.0) then
               numpairsm(windex) = 1
               pairingm(windex,1,1,1) = xmax
               pairingm(windex,1,1,2) = ymax
               pairingm(windex,1,2,1) = xmin
               pairingm(windex,1,2,2) = ymin
               notem(windex,k0) = 1
              goto 45
           else
               do k1 = 1,numpairsm(windex)
                   if(abs(xmax-pairingm(windex,k1,1,1)).lt.1.d0.and.abs(ymax-pairingm(windex,k1,1,2)).lt.1.d0.and. &
&                   abs(xmin-pairingm(windex,k1,2,1)).lt.1.d0.and.abs(ymin-pairingm(windex,k1,2,2)).lt.1.d0) then
                       notem(windex,k0) = k1
                       goto 45
                   endif
              enddo
           endif
           numpairsm(windex) = numpairsm(windex)+1
           pairingm(windex,numpairsm(windex),1,1) = xmax
           pairingm(windex,numpairsm(windex),1,2) = ymax
           pairingm(windex,numpairsm(windex),2,1) = xmin
           pairingm(windex,numpairsm(windex),2,2) = ymin
           notem(windex,k0) = numpairsm(windex)
           endif
45     enddo
      enddo    !!!!!!! high level pairing process end      
      write(*,*) '======= end multi-level pairing ======='
      write(*,*)
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      open(330,file = 'bound_5x.dat',action = 'write')
!      !open(330,file = 'bound_15x.dat',action = 'write')
!      write(330,*) ' variables = "x","y" '
!      do i = 2,ndimx-1
!      do j = 2,ndimy-1
!      if (psi(i,j).ge.1.d0) then
!      force0 = 0
!      endfind = 0
!      posx = 0
!      posy = 0
!      ni = i
!      nj = j
!      switch = 1.d0
!       
!       do while(endfind.eq.0)
!       if((ni.eq.1).or.(ni.eq.ndimx).or.(nj.eq.1).or.(nj.eq.ndimy)) then
!       goto 17
!       endif
!       
!       do i0 = -1,0
!       do j0 = -1,1
!       subdvx(i0,j0) = (psi(ni+i0+1,nj+j0)-psi(ni+i0,nj+j0))/(2*pi/(ndimx-1))
!       enddo
!       enddo
!
!       do i0 = -1,1
!       do j0 = -1,0
!       subdvy(i0,j0) = (psi(ni+i0,nj+j0+1)-psi(ni+i0,nj+j0))/(8*pi/(ndimx-1))
!       enddo
!       enddo
!
!       call endpoints(posx,posy,endx,endy,switch,subdvx,subdvy,endfind,jumpratio,ii,x,y)            
!
!	  do i0 = 1,ii
!           px = ni+x(i0)
!           py = nj+y(i0)
!      if ( px.gt.1.and.px.lt.ndimx.and.py.gt.1.and.py.lt.ndimy ) then
!      if ( notem(5,note(px,py)).eq.notem(5,note(px+1,py)).and.&
!          &notem(5,note(px,py)).eq.notem(5,note(px,py+1)).and.&
!          &notem(5,note(px,py)).eq.notem(5,note(px+1,py+1)) ) then
!          inter(i0) = 0.d0
!      else
!          inter(i0) = 1.d0
!      endif
!      endif !!!!!!! 1<R<ndim
!      enddo
!
!	  do i0 = 1,ii
!      if ( inter(i0).eq.1.d0 ) then
!      if ( i0.gt.1.and.i0.lt.ii) then
!      if ( MOD(i0,5).eq.1 ) then
!           inter(i0) = 1.d0
!           write(330,'(2(1x,f12.6))') ni+x(i0)-1,nj+y(i0)-1
!      endif
!      endif
!      endif
!      enddo
!                  
!      force0 = force0+1
!      if(force0.gt.500) then
!          endfind = 1
!          goto 17
!      endif
!       
!      call tellnew(ni,nj,endx,endy,posx,posy)
!      enddo
!      
!17   force0 = 0
!       endfind = 0
!       posx = 0
!       posy = 0
!       ni = i
!       nj = j
!       switch = -1.d0
!       do while(endfind.eq.0)
!       if((ni.eq.1).or.(ni.eq.ndimx).or.(nj.eq.1).or.(nj.eq.ndimy)) then
!       goto 18
!       endif
!
!       do i0 = -1,0
!       do j0 = -1,1
!       subdvx(i0,j0) = (psi(ni+i0+1,nj+j0)-psi(ni+i0,nj+j0))/(2*pi/(ndimx-1))
!       enddo
!       enddo
!
!       do i0 = -1,1
!       do j0 = -1,0
!       subdvy(i0,j0) = (psi(ni+i0,nj+j0+1)-psi(ni+i0,nj+j0))/(2*pi/(ndimx-1))
!       enddo
!       enddo
!
!       call endpoints(posx,posy,endx,endy,switch,subdvx,subdvy,endfind,jumpratio,ii,x,y)
!       
!	  do i0 = 1,ii
!           px = ni+x(i0)
!           py = nj+y(i0)
!      if ( px.gt.1.and.px.lt.ndimx.and.py.gt.1.and.py.lt.ndimy ) then
!      if ( notem(5,note(px,py)).eq.notem(5,note(px+1,py)).and.&
!          &notem(5,note(px,py)).eq.notem(5,note(px,py+1)).and.&
!          &notem(5,note(px,py)).eq.notem(5,note(px+1,py+1))  )then
!          inter(i0) = 0.d0
!      else
!          inter(i0) = 1.d0
!      endif
!      endif
!      enddo
!
!	  do i0 = 1,ii
!      if ( inter(i0).eq.1.d0 ) then
!      if ( i0.gt.1.and.i0.lt.ii) then
!      if ( MOD(i0,5).eq.1 ) then
!           inter(i0) = 1.d0
!           write(330,'(2(1x,f12.6))') ni+x(i0)-1,nj+y(i0)-1
!      endif
!      endif
!      endif
!      enddo
!
!       force0 = force0+1
!       if(force0.gt.500) then
!       endfind = 1
!       goto 18
!       endif
!
!       call tellnew(ni,nj,endx,endy,posx,posy)
!      enddo
!      
!       endif   
!18    enddo
!      enddo
!      close(330)
      

      write(*,*) '======= begin postprocess ======='
      write(*,*) '_______ begin R = 1~5 windowsize _______'
      !!!!!!!!!!!!!!!!!!!!!!!  let's postprocess  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! dpsimean = 0.d0
      ! lenmean = 0.d0
      ! spsimean = 0.d0
      ! sum = 0.d0  
      ! mpdf = 0.d0
      ! mpdf2 = 0.d0
      ! jpdf = 0.d0
      ! jpdf2 = 0.d0
      ! ypdf = 0.d0
      ! ypdf2 = 0.d0
      ! my = 0.d0
      ! my2 = 0.d0
      ! do jj = 1,5      !!!!!!! overall window size postprocess
      !     do ii = 1,numpairs
      !         locmax = pairing(ii,1)
      !         locmin = pairing(ii,2)
      !         pairsize = pairing(ii,3)
      !         if (locmax.gt.0.and.locmin.gt.0) then
      !             !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
      !             if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
      !                maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
      !                maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
      !                maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
      !                minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
      !                minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
      !                minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
      !                len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
      !                deltapsi = dabs(maxpsi-minpsi)
      !                dpsimean = dpsimean+dabs(deltapsi)
      !                spsimean = spsimean + pairing_tao(ii)/pairing(ii,3)
      !                lenmean = lenmean+len1
      !                sum = sum+1
      !            endif
      !        endif
      !    enddo
      !enddo
      !dpsimean = dpsimean/sum
      !lenmean = lenmean/sum
      !spsimean = spsimean/sum
      !write(*,*) 'dpsimean = ',dpsimean
      !write(*,*) 'lenmean = ',lenmean
      !write(*,*) 'spsimean = ',spsimean
      !write(*,*) 'sum = ',sum
      !
      !do jj = 1,5
      !     do ii = 1,numpairs
      !         locmax = pairing(ii,1)
      !         locmin = pairing(ii,2)
      !         pairsize = pairing(ii,3)
      !         if (locmax.gt.0.and.locmin.gt.0) then
      !             !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
      !             if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
      !                maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
      !                maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
      !                maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
      !                minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
      !                minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
      !                minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
      !                len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
      !                deltapsi = dabs(maxpsi-minpsi)
      !                n1 = 25*len1/lenmean
      !                n1 = min(n1,100)
      !                n1 = max(n1,0)
      !                n2 = 10*deltapsi!/dpsimean
      !                n2 = min(n2,100)
      !                n2 = max(n2,0)
      !                n3 = 5.d0*(pairing_tao(ii)/pairing(ii,3)-25.d0)
      !                n3 = min(n3,100)
      !                n3 = max(n3,0)
      !                jpdf(n1,n2) = jpdf(n1,n2) + pairing(ii,3)
      !                jpdf2(n1,n3) = jpdf2(n1,n3) + pairing(ii,3)
      !                R1_jpdf(n1,n2) = R1_jpdf(n1,n2) + pairing(ii,3)
      !                R1_jpdf2(n1,n3) = R1_jpdf2(n1,n3) + pairing(ii,3)
      !            endif
      !        endif
      !    enddo
      !enddo
      !
      !open(60,file = 'jpdf_deltatao_R1to5_summer.dat')
      !write(60,*) ' variables = "length","deltapsi","pdf" '
      !write(60,*) ' zone I = 101, J = 101, F = POINT '
      !open(61,file = 'jpdf_taomean_R1to5_summer.dat')
      !write(61,*) ' variables = "length","deltapsi","pdf" '
      !write(61,*) ' zone I = 101, J = 101, F = POINT '
      !do j = 0,100
      !    do i = 0,100
      !        write(60,'(3es14.4)') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0,j/10.d0,R1_jpdf(i,j)/m
      !        write(61,'(3es14.4)') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0,j/5.d0+25.d0,R1_jpdf2(i,j)/m
      !    enddo
      !enddo
      !close(60)
      !close(61)      
      !
      !
      write(*,*) '_______ end R = 1~5 windowsize _______'
      
      write(*,*) '_______ begin R = 10~15 windowsize _______'
      !!!!!!!!!!!!!!!!!!!!!!!  let's postprocess  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! dpsimean = 0.d0
      ! lenmean = 0.d0
      ! spsimean = 0.d0
      ! sum = 0.d0  
      ! mpdf = 0.d0
      ! mpdf2 = 0.d0
      ! jpdf = 0.d0
      ! jpdf2 = 0.d0
      ! ypdf = 0.d0
      ! ypdf2 = 0.d0
      ! my = 0.d0
      ! my2 = 0.d0
      ! do jj = 10,15     !!!!!!! overall window size postprocess
      !     do ii = 1,numpairs
      !         locmax = pairing(ii,1)
      !         locmin = pairing(ii,2)
      !         pairsize = pairing(ii,3)
      !         if (locmax.gt.0.and.locmin.gt.0) then
      !             !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
      !             if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
      !                maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
      !                maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
      !                maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
      !                minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
      !                minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
      !                minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
      !                len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
      !                deltapsi = dabs(maxpsi-minpsi)
      !                dpsimean = dpsimean+dabs(deltapsi)
      !                spsimean = spsimean + pairing_tao(ii)/pairing(ii,3)
      !                lenmean = lenmean+len1
      !                sum = sum+1
      !            endif
      !        endif
      !    enddo
      !enddo
      !dpsimean = dpsimean/sum
      !lenmean = lenmean/sum
      !spsimean = spsimean/sum
      !write(*,*) 'dpsimean = ',dpsimean
      !write(*,*) 'lenmean = ',lenmean
      !write(*,*) 'spsimean = ',spsimean
      !write(*,*) 'sum = ',sum
      !
      !do jj = 10,15
      !     do ii = 1,numpairs
      !         locmax = pairing(ii,1)
      !         locmin = pairing(ii,2)
      !         pairsize = pairing(ii,3)
      !         if (locmax.gt.0.and.locmin.gt.0) then
      !             !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
      !             if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
      !                maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
      !                maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
      !                maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
      !                minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
      !                minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
      !                minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
      !                len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
      !                deltapsi = dabs(maxpsi-minpsi)
      !                n1 = 25*len1/lenmean
      !                n1 = min(n1,100)
      !                n1 = max(n1,0)
      !                n2 = 10*deltapsi!/dpsimean
      !                n2 = min(n2,100)
      !                n2 = max(n2,0)
      !                n3 = 5.d0*(pairing_tao(ii)/pairing(ii,3)-25.d0)
      !                n3 = min(n3,100)
      !                n3 = max(n3,0)
      !                jpdf(n1,n2) = jpdf(n1,n2) + pairing(ii,3)
      !                jpdf2(n1,n3) = jpdf2(n1,n3) + pairing(ii,3)
      !                R2_jpdf(n1,n2) = R2_jpdf(n1,n2) + pairing(ii,3)
      !                R2_jpdf2(n1,n3) = R2_jpdf2(n1,n3) + pairing(ii,3)
      !            endif
      !        endif
      !    enddo
      !enddo
      !
      !open(60,file = 'jpdf_deltatao_R10to15_summer.dat')
      !write(60,*) ' variables = "length","deltapsi","pdf" '
      !write(60,*) ' zone I = 101, J = 101, F = POINT '
      !open(61,file = 'jpdf_taomean_R10to15_summer.dat')
      !write(61,*) ' variables = "length","deltapsi","pdf" '
      !write(61,*) ' zone I = 101, J = 101, F = POINT '
      !do j = 0,100
      !    do i = 0,100
      !        write(60,'(3es14.4)') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0,j/10.d0,R2_jpdf(i,j)/m
      !        write(61,'(3es14.4)') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0,j/5.d0+25.d0,R2_jpdf2(i,j)/m
      !    enddo
      !enddo
      !close(60)
      !close(61)           
      write(*,*) '_______ end R = 11~15 windowsize _______'   
      
      write(*,*) '_______ begin overall windowsize _______'
      !!!!!!!!!!!!!!!!!!!!!!!  let's postprocess  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       dpsimean = 0.d0
       lenmean = 0.d0
       spsimean = 0.d0
       sum = 0.d0  
       mpdf = 0.d0
       mpdf2 = 0.d0
       jpdf = 0.d0
       jpdf2 = 0.d0
       ypdf = 0.d0
       ypdf2 = 0.d0
       my = 0.d0
       my2 = 0.d0
       do jj = 1,winnum      !!!!!!! overall window size postprocess
           do ii = 1,numpairs
               locmax = pairing(ii,1)
               locmin = pairing(ii,2)
               pairsize = pairing(ii,3)
               if (locmax.gt.0.and.locmin.gt.0) then
                   !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
                   if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
                      maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
                      maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
                      maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
                      minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
                      minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
                      minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
                      len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
                      deltapsi = dabs(maxpsi-minpsi)
                      dpsimean = dpsimean+dabs(deltapsi)
                      spsimean = spsimean + pairing_tao(ii)/pairing(ii,3)
                      lenmean = lenmean+len1
                      sum = sum+1
                  endif
              endif
          enddo
      enddo
      dpsimean = dpsimean/sum
      lenmean = lenmean/sum
      spsimean = spsimean/sum
      write(*,*) 'dpsimean = ',dpsimean
      write(*,*) 'lenmean = ',lenmean
      write(*,*) 'spsimean = ',spsimean
      write(*,*) 'sum = ',sum
      
      do jj = 1,winnum
           do ii = 1,numpairs
               locmax = pairing(ii,1)
               locmin = pairing(ii,2)
               pairsize = pairing(ii,3)
               if (locmax.gt.0.and.locmin.gt.0) then
                   !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
                   if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
                      maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
                      maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
                      maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
                      minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
                      minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
                      minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
                      len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
                      deltapsi = dabs(maxpsi-minpsi)
                      n1 = 25*len1/lenmean
                      n1 = min(n1,100)
                      n1 = max(n1,0)
                      n2 = 10*deltapsi!/dpsimean
                      n2 = min(n2,100)
                      n2 = max(n2,0)
                      n3 = 5.d0*(pairing_tao(ii)/pairing(ii,3)-25.d0)
                      n3 = min(n3,100)
                      n3 = max(n3,0)
                      jpdf(n1,n2) = jpdf(n1,n2) + pairing(ii,3)
                      jpdf2(n1,n3) = jpdf2(n1,n3) + pairing(ii,3)
                      Rall_jpdf(n1,n2) = Rall_jpdf(n1,n2) + pairing(ii,3)
                      Rall_jpdf2(n1,n3) = Rall_jpdf2(n1,n3) + pairing(ii,3)
                  endif
              endif
          enddo
      enddo
      
      !open(60,file = 'jpdf_deltatao_summer.dat')
      !write(60,*) ' variables = "length","deltapsi","pdf" '
      !write(60,*) ' zone I = 101, J = 101, F = POINT '
      !open(61,file = 'jpdf_taomean_summer.dat')
      !write(61,*) ' variables = "length","deltapsi","pdf" '
      !write(61,*) ' zone I = 101, J = 101, F = POINT '
      !do j = 0,100
      !    do i = 0,100
      !        write(60,'(3es14.4)') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, j/10.d0, Rall_jpdf(i,j)/m
      !        write(61,'(3es14.4)') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0,j/5.d0+25.d0,Rall_jpdf2(i,j)/m
      !    enddo
      !enddo
      !close(60)
      !close(61)

      do i = 0,100
          do j = 0,99
              mpdf(i) = mpdf(i) + jpdf(i,j)
              mpdf2(i) = mpdf2(i) + jpdf2(i,j)
              !ypdf(i) = ypdf(i) + (j/10.d0)**2*jpdf(i,j)
              ypdf(i) = ypdf(i) + (j/10.d0)*jpdf(i,j)
              !ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)**2*jpdf2(i,j)
              ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)*jpdf2(i,j)
          enddo
          if ( mpdf(i).eq.0 ) then
              my(i) = 0.d0
          else
              my(i) = ypdf(i)/mpdf(i)
          endif
          if ( mpdf2(i).eq.0 ) then
              my2(i) = 0.d0
          else
              my2(i) = ypdf2(i)/mpdf2(i)
          endif
          Rall_mpdf(i) = Rall_mpdf(i) + mpdf(i)
          Rall_my(i) = Rall_my(i) + my(i)
          Rall_my2(i) = Rall_my2(i) + my2(i)
      enddo

      open(60,file= 'conmean_deltatao_Rall_summer.dat' ,action='write')
      write(60,*) ' variables = "length", "dtao|length" '
      write(60,*) ' ZONE I=101 , F=POINT '
      open(61,file= 'conmean_taomean_Rall_summer.dat' ,action='write')
      write(61,*) ' variables = "length", "taom|length" '
      write(61,*) ' ZONE I=101 , F=POINT '
      do i=0,100
          write(60,'(2(1x,f13.6,f13.6))') i/25.d0/(2*pi/(ndimx-1))*lenmean*8.d0, Rall_my(i)/m
          write(61,'(2(1x,f13.6,f13.6))') i/25.d0/(2*pi/(ndimx-1))*lenmean*8.d0, Rall_my2(i)/m
      enddo
      close(60)
      close(61)
          
      !open(60,file= './conmean_overall/conmean_deltatao_Rall_summer'//Trim(adjustl(sc))//'.dat' ,action='write')
      !write(60,*) ' variables = "length", "dtao|length" '
      !write(60,*) ' ZONE I=101 , F=POINT '
      !open(61,file= './conmean_overall/conmean_taomean_Rall_summer'//Trim(adjustl(sc))//'.dat' ,action='write')
      !write(61,*) ' variables = "length", "taom|length" '
      !write(61,*) ' ZONE I=101 , F=POINT '
      !do i=0,100
      !    write(60,'(2(1x,f13.6,f13.6))') i/25.d0/(2*pi/(ndimx-1))*lenmean*8.d0, Rall_my(i)
      !    write(61,'(2(1x,f13.6,f13.6))') i/25.d0/(2*pi/(ndimx-1))*lenmean*8.d0, Rall_my2(i)
      !enddo
      !close(60)
      !close(61)
           
      write(*,*) '_______ end overall windowsize _______'  
          
      write(*,*) '_______ begin R = 2dx level _______'
      !!!!!!!!!!!!!!!!!!!!!!!  let's postprocess  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       R = 2
       write(R_str,*) R
       dpsimean = 0.d0
       lenmean = 0.d0
       spsimean = 0.d0
       sum = 0.d0  
       mpdf = 0.d0
       mpdf2 = 0.d0
       jpdf = 0.d0
       jpdf2 = 0.d0
       ypdf = 0.d0
       ypdf2 = 0.d0
       my = 0.d0
       my2 = 0.d0
       do jj = R,R      !!!!!!! overall window size postprocess
           do ii = 1,numpairs
               locmax = pairing(ii,1)
               locmin = pairing(ii,2)
               pairsize = pairing(ii,3)
               if (locmax.gt.0.and.locmin.gt.0) then
                   !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
                   if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
                      maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
                      maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
                      maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
                      minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
                      minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
                      minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
                      len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
                      deltapsi = dabs(maxpsi-minpsi)
                      dpsimean = dpsimean+dabs(deltapsi)
                      spsimean = spsimean + pairing_tao(ii)/pairing(ii,3)
                      lenmean = lenmean+len1
                      sum = sum+1
                  endif
              endif
          enddo
      enddo
      dpsimean = dpsimean/sum
      lenmean = lenmean/sum
      spsimean = spsimean/sum
      write(*,*) 'dpsimean = ',dpsimean
      write(*,*) 'lenmean = ',lenmean
      write(*,*) 'spsimean = ',spsimean
      write(*,*) 'sum = ',sum
      
      do jj = 1,winnum
           do ii = 1,numpairs
               locmax = pairing(ii,1)
               locmin = pairing(ii,2)
               pairsize = pairing(ii,3)
               if (locmax.gt.0.and.locmin.gt.0) then
                   !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
                   if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
                      maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
                      maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
                      maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
                      minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
                      minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
                      minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
                      len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
                      deltapsi = dabs(maxpsi-minpsi)
                      n1 = 25*len1/lenmean
                      n1 = min(n1,100)
                      n1 = max(n1,0)
                      n2 = 10*deltapsi!/dpsimean
                      n2 = min(n2,100)
                      n2 = max(n2,0)
                      n3 = 5.d0*(pairing_tao(ii)/pairing(ii,3)-25.d0)
                      n3 = min(n3,100)
                      n3 = max(n3,0)
                      jpdf(n1,n2) = jpdf(n1,n2) + pairing(ii,3)
                      jpdf2(n1,n3) = jpdf2(n1,n3) + pairing(ii,3)
                  endif
              endif
          enddo
      enddo
      
      do i = 0,100
          do j = 0,99
              mpdf(i) = mpdf(i) + jpdf(i,j)
              mpdf2(i) = mpdf2(i) + jpdf2(i,j)
              !ypdf(i) = ypdf(i) + (j/10.d0)**2*jpdf(i,j)
              !ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)**2*jpdf2(i,j)
              ypdf(i) = ypdf(i) + (j/10.d0)*jpdf(i,j)
              ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)*jpdf2(i,j)
          enddo
          if ( mpdf(i).eq.0 ) then
              my(i) = 0.d0
          else
              my(i) = ypdf(i)/mpdf(i)
          endif
          if ( mpdf2(i).eq.0 ) then
              my2(i) = 0.d0
          else
              my2(i) = ypdf2(i)/mpdf2(i)
          endif
          R2_my(i) = R2_my(i) + my(i)
          R2_my2(i) = R2_my2(i) + my2(i)
      enddo      
      
      open(60,file= 'conmean_deltatao_R2_summer.dat' ,action='write')
      write(60,*) ' variables = "length", "dtao|length" '
      write(60,*) ' ZONE I=101 , F=POINT '
      open(61,file= 'conmean_taomean_R2_summer.dat' ,action='write')
      write(61,*) ' variables = "length", "taom|length" '
      write(61,*) ' ZONE I=101 , F=POINT '
      do i=0,100
          write(60,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R2_my(i)/m
          write(61,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R2_my2(i)/m
      enddo
      close(60)
      close(61)      
      
      !open(60,file= './conmean_overall/conmean_deltatao_R2_summer'//TRim(adjustl(sc))//'.dat' ,action='write')
      !write(60,*) ' variables = "length", "dtao|length" '
      !write(60,*) ' ZONE I=101 , F=POINT '
      !open(61,file= './conmean_overall/conmean_taomean_R2_summer'//TRim(adjustl(sc))//'.dat' ,action='write')
      !write(61,*) ' variables = "length", "taom|length" '
      !write(61,*) ' ZONE I=101 , F=POINT '
      !do i=0,100
      !    write(60,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R2_my(i)
      !    write(61,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R2_my2(i)
      !enddo
      !close(60)
      !close(61)
      
      write(*,*) '_______ end R = 2dx level _______'
      
      !write(*,*) '_______ begin R = 4dx level _______'
      !!!!!!!!!!!!!!!!!!!!!!!!  let's postprocess  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! R = 4
      ! write(R_str,*) R
      ! dpsimean = 0.d0
      ! lenmean = 0.d0
      ! spsimean = 0.d0
      ! sum = 0.d0  
      ! mpdf = 0.d0
      ! mpdf2 = 0.d0
      ! jpdf = 0.d0
      ! jpdf2 = 0.d0
      ! ypdf = 0.d0
      ! ypdf2 = 0.d0
      ! my = 0.d0
      ! my2 = 0.d0
      ! do jj = R,R      !!!!!!! overall window size postprocess
      !     do ii = 1,numpairs
      !         locmax = pairing(ii,1)
      !         locmin = pairing(ii,2)
      !         pairsize = pairing(ii,3)
      !         if (locmax.gt.0.and.locmin.gt.0) then
      !             !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
      !             if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
      !                maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
      !                maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
      !                maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
      !                minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
      !                minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
      !                minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
      !                len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
      !                deltapsi = dabs(maxpsi-minpsi)
      !                dpsimean = dpsimean+dabs(deltapsi)
      !                spsimean = spsimean + pairing_tao(ii)/pairing(ii,3)
      !                lenmean = lenmean+len1
      !                sum = sum+1
      !            endif
      !        endif
      !    enddo
      !enddo
      !dpsimean = dpsimean/sum
      !lenmean = lenmean/sum
      !spsimean = spsimean/sum
      !write(*,*) 'dpsimean = ',dpsimean
      !write(*,*) 'lenmean = ',lenmean
      !write(*,*) 'spsimean = ',spsimean
      !write(*,*) 'sum = ',sum
      !
      !do jj = 1,winnum
      !     do ii = 1,numpairs
      !         locmax = pairing(ii,1)
      !         locmin = pairing(ii,2)
      !         pairsize = pairing(ii,3)
      !         if (locmax.gt.0.and.locmin.gt.0) then
      !             !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
      !             if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
      !                maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
      !                maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
      !                maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
      !                minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
      !                minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
      !                minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
      !                len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
      !                deltapsi = dabs(maxpsi-minpsi)
      !                n1 = 25*len1/lenmean
      !                n1 = min(n1,100)
      !                n1 = max(n1,0)
      !                n2 = 10*deltapsi!/dpsimean
      !                n2 = min(n2,100)
      !                n2 = max(n2,0)
      !                n3 = 5.d0*(pairing_tao(ii)/pairing(ii,3)-25.d0)
      !                n3 = min(n3,100)
      !                n3 = max(n3,0)
      !                jpdf(n1,n2) = jpdf(n1,n2) + pairing(ii,3)
      !                jpdf2(n1,n3) = jpdf2(n1,n3) + pairing(ii,3)
      !            endif
      !        endif
      !    enddo
      !enddo
      !
      !do i = 0,100
      !    do j = 0,99
      !        mpdf(i) = mpdf(i) + jpdf(i,j)
      !        mpdf2(i) = mpdf2(i) + jpdf2(i,j)
      !        ypdf(i) = ypdf(i) + j/10.d0*jpdf(i,j)
      !        ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)*jpdf2(i,j)
      !    enddo
      !    if ( mpdf(i).eq.0 ) then
      !        my(i) = 0.d0
      !    else
      !        my(i) = ypdf(i)/mpdf(i)
      !    endif
      !    if ( mpdf2(i).eq.0 ) then
      !        my2(i) = 0.d0
      !    else
      !        my2(i) = ypdf2(i)/mpdf2(i)
      !    endif
      !    R4_my(i) = R4_my(i) + my(i)
      !    R4_my2(i) = R4_my2(i) + my2(i)
      !enddo      
      !
      !!open(60,file= 'conmean_deltatao_R4_summer.dat' ,action='write')
      !!write(60,*) ' variables = "length", "dtao|length" '
      !!write(60,*) ' ZONE I=101 , F=POINT '
      !!open(61,file= 'conmean_taomean_R4_summer.dat' ,action='write')
      !!write(61,*) ' variables = "length", "taom|length" '
      !!write(61,*) ' ZONE I=101 , F=POINT '
      !!do i=0,100
      !!    write(60,'(2(1x,f13.6,f13.6))') i/25.d0, R4_my(i)/m
      !!    write(61,'(2(1x,f13.6,f13.6))') i/25.d0, R4_my2(i)/m
      !!enddo
      !!close(60)
      !!close(61)     
      !write(*,*) '_______ end R = 4dx level _______'
      
      write(*,*) '_______ begin R = 6dx level _______'
      !!!!!!!!!!!!!!!!!!!!!!!  let's postprocess  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       R = 6
       write(R_str,*) R
       dpsimean = 0.d0
       lenmean = 0.d0
       spsimean = 0.d0
       sum = 0.d0  
       mpdf = 0.d0
       mpdf2 = 0.d0
       jpdf = 0.d0
       jpdf2 = 0.d0
       ypdf = 0.d0
       ypdf2 = 0.d0
       my = 0.d0
       my2 = 0.d0
       do jj = R,R      !!!!!!! overall window size postprocess
           do ii = 1,numpairs
               locmax = pairing(ii,1)
               locmin = pairing(ii,2)
               pairsize = pairing(ii,3)
               if (locmax.gt.0.and.locmin.gt.0) then
                   !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
                   if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
                      maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
                      maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
                      maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
                      minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
                      minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
                      minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
                      len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
                      deltapsi = dabs(maxpsi-minpsi)
                      dpsimean = dpsimean+dabs(deltapsi)
                      spsimean = spsimean + pairing_tao(ii)/pairing(ii,3)
                      lenmean = lenmean+len1
                      sum = sum+1
                  endif
              endif
          enddo
      enddo
      dpsimean = dpsimean/sum
      lenmean = lenmean/sum
      spsimean = spsimean/sum
      write(*,*) 'dpsimean = ',dpsimean
      write(*,*) 'lenmean = ',lenmean
      write(*,*) 'spsimean = ',spsimean
      write(*,*) 'sum = ',sum
      
      do jj = 1,winnum
           do ii = 1,numpairs
               locmax = pairing(ii,1)
               locmin = pairing(ii,2)
               pairsize = pairing(ii,3)
               if (locmax.gt.0.and.locmin.gt.0) then
                   !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
                   if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
                      maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
                      maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
                      maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
                      minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
                      minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
                      minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
                      len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
                      deltapsi = dabs(maxpsi-minpsi)
                      n1 = 25*len1/lenmean
                      n1 = min(n1,100)
                      n1 = max(n1,0)
                      n2 = 10*deltapsi!/dpsimean
                      n2 = min(n2,100)
                      n2 = max(n2,0)
                      n3 = 5.d0*(pairing_tao(ii)/pairing(ii,3)-25.d0)
                      n3 = min(n3,100)
                      n3 = max(n3,0)
                      jpdf(n1,n2) = jpdf(n1,n2) + pairing(ii,3)
                      jpdf2(n1,n3) = jpdf2(n1,n3) + pairing(ii,3)
                  endif
              endif
          enddo
      enddo
      
      do i = 0,100
          do j = 0,99
              mpdf(i) = mpdf(i) + jpdf(i,j)
              mpdf2(i) = mpdf2(i) + jpdf2(i,j)
              !ypdf(i) = ypdf(i) + (j/10.d0)**2*jpdf(i,j)
              !ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)**2*jpdf2(i,j)
              ypdf(i) = ypdf(i) + (j/10.d0)*jpdf(i,j)
              ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)*jpdf2(i,j)
          enddo
          if ( mpdf(i).eq.0 ) then
              my(i) = 0.d0
          else
              my(i) = ypdf(i)/mpdf(i)
          endif
          if ( mpdf2(i).eq.0 ) then
              my2(i) = 0.d0
          else
              my2(i) = ypdf2(i)/mpdf2(i)
          endif
          R6_my(i) = R6_my(i) + my(i)
          R6_my2(i) = R6_my2(i) + my2(i)
      enddo      
      
      open(60,file= 'conmean_deltatao_R6_summer.dat' ,action='write')
      write(60,*) ' variables = "length", "dtao|length" '
      write(60,*) ' ZONE I=101 , F=POINT '
      open(61,file= 'conmean_taomean_R6_summer.dat' ,action='write')
      write(61,*) ' variables = "length", "taom|length" '
      write(61,*) ' ZONE I=101 , F=POINT '
      do i=0,100
          write(60,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R6_my(i)/m
          write(61,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R6_my2(i)/m
      enddo
      close(60)
      close(61)  
      
      !open(60,file= './conmean_overall/conmean_deltatao_R6_summer'//TRim(adjustl(sc))//'.dat' ,action='write')
      !write(60,*) ' variables = "length", "dtao|length" '
      !write(60,*) ' ZONE I=101 , F=POINT '
      !open(61,file= './conmean_overall/conmean_taomean_R6_summer'//TRim(adjustl(sc))//'.dat' ,action='write')
      !write(61,*) ' variables = "length", "taom|length" '
      !write(61,*) ' ZONE I=101 , F=POINT '
      !do i=0,100
      !    write(60,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R6_my(i)
      !    write(61,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R6_my2(i)
      !enddo
      !close(60)
      !close(61)  
      
      write(*,*) '_______ end R = 6dx level _______'
      
      !write(*,*) '_______ begin R = 8dx level _______'
      !!!!!!!!!!!!!!!!!!!!!!!!  let's postprocess  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! R = 8
      ! write(R_str,*) R
      ! dpsimean = 0.d0
      ! lenmean = 0.d0
      ! spsimean = 0.d0
      ! sum = 0.d0  
      ! mpdf = 0.d0
      ! mpdf2 = 0.d0
      ! jpdf = 0.d0
      ! jpdf2 = 0.d0
      ! ypdf = 0.d0
      ! ypdf2 = 0.d0
      ! my = 0.d0
      ! my2 = 0.d0
      ! do jj = R,R      !!!!!!! overall window size postprocess
      !     do ii = 1,numpairs
      !         locmax = pairing(ii,1)
      !         locmin = pairing(ii,2)
      !         pairsize = pairing(ii,3)
      !         if (locmax.gt.0.and.locmin.gt.0) then
      !             !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
      !             if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
      !                maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
      !                maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
      !                maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
      !                minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
      !                minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
      !                minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
      !                len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
      !                deltapsi = dabs(maxpsi-minpsi)
      !                dpsimean = dpsimean+dabs(deltapsi)
      !                spsimean = spsimean + pairing_tao(ii)/pairing(ii,3)
      !                lenmean = lenmean+len1
      !                sum = sum+1
      !            endif
      !        endif
      !    enddo
      !enddo
      !dpsimean = dpsimean/sum
      !lenmean = lenmean/sum
      !spsimean = spsimean/sum
      !write(*,*) 'dpsimean = ',dpsimean
      !write(*,*) 'lenmean = ',lenmean
      !write(*,*) 'spsimean = ',spsimean
      !write(*,*) 'sum = ',sum
      !
      !do jj = 1,winnum
      !     do ii = 1,numpairs
      !         locmax = pairing(ii,1)
      !         locmin = pairing(ii,2)
      !         pairsize = pairing(ii,3)
      !         if (locmax.gt.0.and.locmin.gt.0) then
      !             !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
      !             if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
      !                maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
      !                maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
      !                maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
      !                minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
      !                minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
      !                minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
      !                len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
      !                deltapsi = dabs(maxpsi-minpsi)
      !                n1 = 25*len1/lenmean
      !                n1 = min(n1,100)
      !                n1 = max(n1,0)
      !                n2 = 10*deltapsi!/dpsimean
      !                n2 = min(n2,100)
      !                n2 = max(n2,0)
      !                n3 = 5.d0*(pairing_tao(ii)/pairing(ii,3)-25.d0)
      !                n3 = min(n3,100)
      !                n3 = max(n3,0)
      !                jpdf(n1,n2) = jpdf(n1,n2) + pairing(ii,3)
      !                jpdf2(n1,n3) = jpdf2(n1,n3) + pairing(ii,3)
      !            endif
      !        endif
      !    enddo
      !enddo
      !
      !do i = 0,100
      !    do j = 0,99
      !        mpdf(i) = mpdf(i) + jpdf(i,j)
      !        mpdf2(i) = mpdf2(i) + jpdf2(i,j)
      !        ypdf(i) = ypdf(i) + j/10.d0*jpdf(i,j)
      !        ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)*jpdf2(i,j)
      !    enddo
      !    if ( mpdf(i).eq.0 ) then
      !        my(i) = 0.d0
      !    else
      !        my(i) = ypdf(i)/mpdf(i)
      !    endif
      !    if ( mpdf2(i).eq.0 ) then
      !        my2(i) = 0.d0
      !    else
      !        my2(i) = ypdf2(i)/mpdf2(i)
      !    endif
      !    R8_my(i) = R8_my(i) + my(i)
      !    R8_my2(i) = R8_my2(i) + my2(i)
      !enddo      
      !
      !open(60,file= 'conmean_deltatao_R8_summer.dat' ,action='write')
      !write(60,*) ' variables = "length", "dtao|length" '
      !write(60,*) ' ZONE I=101 , F=POINT '
      !open(61,file= 'conmean_taomean_R8_summer.dat' ,action='write')
      !write(61,*) ' variables = "length", "taom|length" '
      !write(61,*) ' ZONE I=101 , F=POINT '
      !do i=0,100
      !    write(60,'(2(1x,f13.6,f13.6))') i/25.d0, R8_my(i)/m
      !    write(61,'(2(1x,f13.6,f13.6))') i/25.d0, R8_my2(i)/m
      !enddo
      !close(60)
      !close(61) 
      !write(*,*) '_______ end R = 8dx level _______'
      
      write(*,*) '_______ begin R = 10dx level _______'
      !!!!!!!!!!!!!!!!!!!!!!!  let's postprocess  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       R = 10
       write(R_str,*) R
       dpsimean = 0.d0
       lenmean = 0.d0
       spsimean = 0.d0
       sum = 0.d0  
       mpdf = 0.d0
       mpdf2 = 0.d0
       jpdf = 0.d0
       jpdf2 = 0.d0
       ypdf = 0.d0
       ypdf2 = 0.d0
       my = 0.d0
       my2 = 0.d0
       do jj = R,R      !!!!!!! overall window size postprocess
           do ii = 1,numpairs
               locmax = pairing(ii,1)
               locmin = pairing(ii,2)
               pairsize = pairing(ii,3)
               if (locmax.gt.0.and.locmin.gt.0) then
                   !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
                   if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
                      maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
                      maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
                      maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
                      minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
                      minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
                      minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
                      len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
                      deltapsi = dabs(maxpsi-minpsi)
                      dpsimean = dpsimean+dabs(deltapsi)
                      spsimean = spsimean + pairing_tao(ii)/pairing(ii,3)
                      lenmean = lenmean+len1
                      sum = sum+1
                  endif
              endif
          enddo
      enddo
      dpsimean = dpsimean/sum
      lenmean = lenmean/sum
      spsimean = spsimean/sum
      write(*,*) 'dpsimean = ',dpsimean
      write(*,*) 'lenmean = ',lenmean
      write(*,*) 'spsimean = ',spsimean
      write(*,*) 'sum = ',sum
      
      do jj = 1,winnum
           do ii = 1,numpairs
               locmax = pairing(ii,1)
               locmin = pairing(ii,2)
               pairsize = pairing(ii,3)
               if (locmax.gt.0.and.locmin.gt.0) then
                   !  if ((locmax.ne.oldmax).or.(locmin.ne.oldmin)) then
                   if (multiscale(jj,locmax,1).gt.0.and.multiscale(jj,locmin,1).gt.0) then
                      maxxcor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordx)
                      maxycor = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),coordy) 
                      maxpsi = interp(multiscale(jj,locmax,1),multiscale(jj,locmax,2),psi)
                      minxcor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordx)
                      minycor = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),coordy)
                      minpsi = interp(multiscale(jj,locmin,1),multiscale(jj,locmin,2),psi)
                      len1 = sqrt((minxcor-maxxcor)**2+(minycor-maxycor)**2)
                      deltapsi = dabs(maxpsi-minpsi)
                      n1 = 25*len1/lenmean
                      n1 = min(n1,100)
                      n1 = max(n1,0)
                      n2 = 10*deltapsi!/dpsimean
                      n2 = min(n2,100)
                      n2 = max(n2,0)
                      n3 = 5.d0*(pairing_tao(ii)/pairing(ii,3)-25.d0)
                      n3 = min(n3,100)
                      n3 = max(n3,0)
                      jpdf(n1,n2) = jpdf(n1,n2) + pairing(ii,3)
                      jpdf2(n1,n3) = jpdf2(n1,n3) + pairing(ii,3)
                  endif
              endif
          enddo
      enddo
      
      do i = 0,100
          do j = 0,99
              mpdf(i) = mpdf(i) + jpdf(i,j)
              mpdf2(i) = mpdf2(i) + jpdf2(i,j)
              !ypdf(i) = ypdf(i) + (j/10.d0)**2*jpdf(i,j)
              !ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)**2*jpdf2(i,j)
              ypdf(i) = ypdf(i) + (j/10.d0)*jpdf(i,j)
              ypdf2(i) = ypdf2(i) + (j/5.d0+25.d0)*jpdf2(i,j)
          enddo
          if ( mpdf(i).eq.0 ) then
              my(i) = 0.d0
          else
              my(i) = ypdf(i)/mpdf(i)
          endif
          if ( mpdf2(i).eq.0 ) then
              my2(i) = 0.d0
          else
              my2(i) = ypdf2(i)/mpdf2(i)
          endif
          R10_my(i) = R10_my(i) + my(i)
          R10_my2(i) = R10_my2(i) + my2(i)
      enddo  
      
      open(60,file= 'conmean_deltatao_R10_summer.dat' ,action='write')
      write(60,*) ' variables = "length", "dtao|length" '
      write(60,*) ' ZONE I=101 , F=POINT '
      open(61,file= 'conmean_taomean_R10_summer.dat' ,action='write')
      write(61,*) ' variables = "length", "taom|length" '
      write(61,*) ' ZONE I=101 , F=POINT '
      do i=0,100
          write(60,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R10_my(i)/m
          write(61,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R10_my2(i)/m
      enddo
      close(60)
      close(61) 
      
      !open(60,file= './conmean_overall/conmean_deltatao_R10_summer'//TRim(adjustl(sc))//'.dat' ,action='write')
      !write(60,*) ' variables = "length", "dtao|length" '
      !write(60,*) ' ZONE I=101 , F=POINT '
      !open(61,file= './conmean_overall/conmean_taomean_R10_summer'//TRim(adjustl(sc))//'.dat' ,action='write')
      !write(61,*) ' variables = "length", "taom|length" '
      !write(61,*) ' ZONE I=101 , F=POINT '
      !do i=0,100
      !    write(60,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R10_my(i)
      !    write(61,'(2(1x,f13.6,f13.6))') i/25.d0*lenmean/(2*pi/(ndimx-1))*8.d0, R10_my2(i)
      !enddo
      !close(60)
      !close(61)
      
      write(*,*) '_______ end R = 10dx level _______'

      enddo
      
      deallocate(multiscale)      
      deallocate(dx)
      deallocate(dy)      
      
56 format(i3,2es15.4)
      call CPU_TIME(tdiff1)
      write(*,*) 'diff_scale_duration = ',tdiff1-tdiff0
      end

       subroutine sphere(ni,nj,endx,endy,psimax,Rball,window,switch,psiball,xcorball,ycorball,meetbound) 
      implicit none
      include "param.h"
      double precision,external::lamda
      double precision psimax,psiball,xcorball,ycorball,R
      double precision epsilon,endx,endy,exmpsi,window,switch
      double precision dx,dy,exmx1,exmx2,lam,lampsi
      double precision distance(1:ndimx,1:ndimy)
      double precision, dimension(:,:),pointer:: coordx,coordy

      double precision, dimension(:,:),pointer:: psi
      double precision, dimension(:,:),pointer:: points
      integer pairing(1000000,3)
      integer status,numends,loop,px,py,zx,zy,i,j,ni,nj
      integer meetbound,Rball
      common/large0/ psi
      common/large1/ pairing
      common/large2/ coordx,coordy
      common/large5/ points
      epsilon = 1.d-4
      R = window
      meetbound = 0
      ! ======= psimax mean the local extremal(max/min) in the cubic
      psiball = psimax
      xcorball = ni+endx
      ycorball = nj+endy
      ! ======= psiball means the extremal point in the ball with radius R
      px = ni+endx
      py = nj+endy
      dx = ni+endx-px
      dy = nj+endy-py
      
      if (((xcorball-Rball).lt.1.d0).or.((ycorball-Rball).lt.1.d0).or.((ndimx-xcorball).lt.Rball).or.((ndimy-ycorball).lt.Rball)) then
          meetbound = 1
          go to 80
      else
          exmx1 = (ni+endx)*2*pi/(ndimx-1)
          exmx2 = (nj+endy)*2*pi/(ndimy-1)
          do j = max(1,nj-Rball),min(nj+Rball,ndimy-1)
          do i = max(1,ni-Rball),min(ni+Rball,ndimx-1)
              distance(i,j) = sqrt((coordx(i,j)-exmx1)**2+(coordy(i,j)-exmx2)**2)
          enddo
          enddo
          ! ======= calculate the max/min inside the ball
          do j = max(1,nj-Rball),min(nj+Rball+1,ndimy)
          do i = max(1,ni-Rball),min(ni+Rball+1,ndimx)
              if (distance(i,j).le.R) then
              if ((psi(i,j)*switch).gt.(psiball*switch)) then
                  psiball = psi(i,j)
                  xcorball = i
                  ycorball = j
              endif
              endif
          enddo
          enddo
          ! ======= calculate the min/max in the surface  
          do j = max(1,nj-Rball),min(nj+Rball,ndimy-1)
          do i = max(1,ni-Rball),min(ni+Rball,ndimx-1)
              ! ======= grid in x dirction
              if((distance(i,j)-R)*(distance(i+1,j)-R).le.0) then
                  lam = lamda(exmx1,exmx2,coordx(i,j),coordy(i,j),coordx(i+1,j),coordy(i+1,j),window)
                  if(lam.ge.0) then
                      lampsi = psi(i,j)+lam*(psi(i+1,j)-psi(i,j))
                      if ((switch*lampsi).gt.(switch*psiball)) then
                          psiball = lampsi
                          xcorball = i+lam
                          ycorball = j
                      endif
                  endif
              endif
              ! ======= grid in y dirction
              if((distance(i,j)-R)*(distance(i,j+1)-R).le.0) then
                  lam = lamda(exmx1,exmx2,coordx(i,j),coordy(i,j),coordx(i,j+1),coordy(i,j+1),window)
                  if(lam.ge.0) then
                      lampsi = psi(i,j)+lam*(psi(i,j+1)-psi(i,j))
                      if ((switch*lampsi).gt.(switch*psiball)) then
                          psiball = lampsi
                          xcorball = i
                          ycorball = j+lam
                      endif
                  endif
              endif
          enddo
          enddo
      endif

80      end subroutine

       function lamda(x0,y0,x1,y1,x2,y2,r)    !在不同尺度下搜寻时，球面（圆）上边界点距最近网格点距离 
       implicit none
       double precision x0,y0,x1,y1,x2,y2,r,d,l,cross
       double precision vect1x,vect1y,vect2x,vect2y,lamda1,lamda2,lamda
       d = dsqrt((x1-x0)**2+(y1-y0)**2)
       l = dsqrt((x2-x0)**2+(y2-y0)**2)
       if (((d.le.r).and.(l.gt.r)).or.((d.gt.r).and.(l.le.r))) then
           ! ======= vector OM
           vect1x = x1-x0
           vect1y = y1-y0
           ! ======= vector ON
           vect2x = x2-x0
           vect2y = y2-y0
           cross = vect1x*vect2x+vect1y*vect2y
           lamda1 = (d**2-cross+dsqrt((cross-d**2)**2-(d**2+l**2-2*cross)*(d**2-r**2)))/(d**2+l**2-2*cross)
           if ((lamda1.ge.0).and.(lamda1.le.1)) lamda = lamda1
           lamda2 = (d**2-cross-dsqrt((cross-d**2)**2-(d**2+l**2-2*cross)*(d**2-r**2)))/(d**2+l**2-2*cross)
           if ((lamda2.ge.0).and.(lamda2.le.1)) lamda = lamda2
       else
           lamda = -100.d0
       endif
       return
       end function

       function interp(tx,ty,phi)
       implicit none
       include "param.h"
       double precision tx,ty,dx,dy,interp
       double precision phi(ndimx,ndimy)
       integer px,py,status
       
       px = tx
       py = ty
       dx = tx-px
       dy = ty-py       
       interp=(1.d0-dx)*(1.d0-dy)*phi(px,py)+dx*(1.d0-dy)*phi(px+1,py)+(1.d0-dx)*dy*phi(px,py+1)+dx*dy*phi(px+1,py+1)
       
    end function

      subroutine endpoint(posx,posy,endx,endy,switch,subdvx,subdvy,endfind,jumpratio)
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
    
      subroutine endpoints(posx,posy,endx,endy,switch,subdvx,subdvy,endfind,jumpratio,ii,x,y)
       implicit none
      double precision jumpratio
      double precision posx,posy,nposx,nposy,endx,endy
      double precision switch,angle(2),temp,temp1,temp2
      double precision subdvx(-1:0,-1:1),subdvy(-1:1,-1:0),x(1000000),y(1000000)
      double precision normx,normy,normt,jump
      integer ix,iy,inner,endfind,force,ii
      double precision xp,xm,yp,ym
      double precision vxp(3),vxm(3),vyp(3),vym(3)
      double precision diverg,totv,vsum(3)

      include "param.h" 
      
      inner = 1
      force = 0

      nposx = posx
      nposy = posy
      
      normx = 0.d0
      normy = 0.d0
      do ix = -1,0
      do iy = -1,1
      normx = max(dabs(subdvx(ix,iy)),normx)
      enddo
      enddo
      do ix = -1,1
      do iy = -1,0
      normy = max(dabs(subdvy(ix,iy)),normy)
      enddo
      enddo
      normt = dsqrt(normx*normx+normy*normy)+1.d-20
      
      ii = 0
      do while(inner.eq.1)
      ii = ii+1
      x(ii) = nposx
      y(ii) = nposy
      if((abs(nposx).gt.0.5d0).or.(abs(nposy).gt.0.5d0)) then
      inner = 0
      endx = nposx
      endy = nposy
      goto 64
      endif

      call findangle(nposx,nposy,switch,subdvx,subdvy,angle,totv)

      jump = min(pace,totv/normt/4.d0)

      if(jump.lt.jumpeps) then
      xp = nposx+pace
      xm = nposx-pace
      yp = nposy+pace
      ym = nposy-pace
       call findangle(xp,nposy,switch,subdvx,subdvy,vxp,temp)
       call findangle(xm,nposy,switch,subdvx,subdvy,vxm,temp)
       call findangle(nposx,yp,switch,subdvx,subdvy,vyp,temp)
       call findangle(nposx,ym,switch,subdvx,subdvy,vym,temp)
      diverg = vxp(1)-vxm(1)+vyp(2)-vym(2)
      if(diverg.gt.-2.5) then
      call random_number(vsum(1))
      call random_number(vsum(2))
             
       if(abs(angle(1)).ge.abs(angle(2))) then
       temp2 = vsum(2)**2*angle(2)
       temp1 = 10*temp2*angle(2)/angle(1)
       temp = 1.d-20+dsqrt(temp1**2+temp2**2)
       nposx = nposx-temp1/temp*pace
       nposy = nposy+temp2/temp*pace       
       goto 65
       else if(abs(angle(2)).ge.abs(angle(1))) then
       temp1 = vsum(1)**2*angle(1)
       temp2 = 10*temp1*angle(1)/angle(2)
       temp = 1.d-20+dsqrt(temp1**2+temp2**2)
       nposx = nposx+temp1/temp*pace
       nposy = nposy-temp2/temp*pace
       goto 65
       endif
      endif

      endfind = 1     
      inner = 0
      endx = nposx
      endy = nposy
      goto 64
     
      else
       nposx = nposx+angle(1)*jump
       nposy = nposy+angle(2)*jump
      endif

  65  force = force+1
      if(force.gt.100) then
      jumpratio = jumpratio+1
      endfind = 1
      inner = 0
      endx = nposx
      endy = nposy
      goto 64
      endif
  
      enddo
      
  64  return
      
      end
    
      subroutine findangle(nposx,nposy,switch,subdvx,subdvy,angl,totv)
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
      angl(1)=(1.d0-sx)*(1.d0-sy)*subdvx(ix,iy)+sx*(1.d0-sy)*subdvx(ix+1,iy)+&
          & (1.d0-sx)*sy*subdvx(ix,iy+1)+sx*sy*subdvx(ix+1,iy+1)
      angl(1)=angl(1)*switch

      ix=0
      iy=-1
      sx=nposx
      sy=nposy+0.5d0
      if(nposx.lt.0.d0) then
       sx=1.d0+nposx      
       ix=-1
      endif
      angl(2)=(1.d0-sx)*(1.d0-sy)*subdvy(ix,iy)+sx*(1.d0-sy)*subdvy(ix+1,iy)+&
          & (1.d0-sx)*sy*subdvy(ix,iy+1)+sx*sy*subdvy(ix+1,iy+1)
      angl(2)=angl(2)*switch
      
      totv=1.d-20+dsqrt(angl(1)**2+angl(2)**2)
      angl(1)=1.d-20+angl(1)/totv
      angl(2)=1.d-20+angl(2)/totv
      
      return
      end
    
      subroutine tellnew(ni,nj,endx,endy,posx,posy)
      implicit none
      integer ni,nj
      double precision endx,endy,posx,posy

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
