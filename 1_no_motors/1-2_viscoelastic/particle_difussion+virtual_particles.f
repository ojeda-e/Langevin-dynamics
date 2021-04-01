	program eqsol

	implicit none
      
c-------parametros
	integer*8 nmod,nc,nream,nmaxl 
	parameter(nmod=9)		!numero de modos
	parameter(nc=60)  		!numero de puntos lag time
	parameter(nream=2000) 		!numero maximo de realizaciones
	parameter(nmaxl=20000000) 	!numero maximo de lineas a leer tray

c-------variables globales
	integer*4 idttr	
	integer*8 ir,nrea,numline
	integer*8 itmax,i,j
	real*4 dtsamp
	real*8 bi
	real*16 dt
	

c-------variables particulas virtuales
	real*8 alpha,b0,g0,nu0,c
        real*8 sqetdt(nmod),sqet0dt,sumeta,xsamp(nmaxl)
	real*8 x,u(nmod),KT,et(nmod),nu(nmod),k(nmod) !array con tamaño
	real*8 et0,rnd,gammln,s,gasdev
	integer*4 iseedg
	integer*8 it

c-------variables msd
	integer m,l3,lr,nr,p	
	real*8 corrtmp(3),avcorr(nc,3),corrtr(nream,nc,3),corr(nc,3)
	real*8 f(nmaxl)

c-------variables motores
	real*8 Fm

c-------Variables nombre archivos
	character*3 strnumF
      	character*4 namefile
      	character*12 filemsd,indat
      	character*16 filemsdtr
	character*21 ftrj



        open(2,file='traj.ini')
        read(2,*) namefile
        read(2,*) nrea
        read(2,*) dt
	read(2,*) iseedg
	read(2,*) KT
	read(2,*) alpha
	read(2,*) et0
	read(2,*) b0
	read(2,*) g0
	read(2,*) nu0
	read(2,*) c
	read(2,*) idttr
	read(2,*) x
	read(2,*) Fm
	read(2,*) itmax
	read(2,*) dtsamp
        close(2)
	
        indat=''//namefile//'ini.txt'
	
	open(2,file=indat)
 	  write(2,*) 'namefile=',namefile
          write(2,*) 'nrea= ',nrea
          write(2,*) 'iseedg=',iseedg
	  write(2,*) 'KT=',KT
    	  write(2,*) 'alpha=',alpha
	  write(2,*) 'et0=',et0
 	  write(2,*) 'bo=',b0
 	  write(2,*) 'g0=',g0
	  write(2,*) 'nu0=',nu0
	  write(2,*) 'c=',c
	  write(2,*) 'idttr=',idttr
	  write(2,*) 'x=',x
	  write(2,*) 'Fm=',Fm
	  write(2,*) 'itmax',itmax 
	  write(2,*) 'dtsampleo',dtsamp
	close(2)

	filemsd='msd-'//namefile//'.dat'
	filemsdtr='msd-tray'//namefile//'.dat'



	open(4,file=filemsd)
	open(5,file=filemsdtr)


	numline=0
	do ir=1,nrea

        write(strnumF,'(I3)') 100

        ftrj=''//namefile//'-'//'F'//strnumF//'.dat'
	 
	open(3,file=ftrj)
          sumeta=0
          sqet0dt=sqrt(2*KT/et0*dt)	   	
	  
	  do i=1,nmod
           nu(i)=nu0/b0**(i-1)   
           k(i)=g0/dexp(gammln(1.d0-alpha))*nu(i)**alpha*c
           et(i)=k(i)/nu(i)  
           sqetdt(i)=sqrt(2*KT/et(i)*dt)
           sumeta=sumeta+et(i)
	  end do

       	  
          x=0.d0        
          do it=1,itmax
            s=Fm
            do i=1,nmod
             s=s-k(i)*(x-u(i))
             rnd=gasdev(iseedg)
             u(i)=u(i)+sqetdt(i)*rnd+dt*nu(i)*(x-u(i))
            end do
            rnd=gasdev(iseedg)
            x=x+sqet0dt*rnd+s*dt/et0
  	     if(mod(it,idttr).eq.0) then     
		write(3,*) sngl(x)
		xsamp(numline)=x
		numline=numline+1
             end if	      
          end do   !it=1,itmax (ciclo temporal)
	close(3)
	write(*,*) numline
	
	if(numline.gt.nmaxl) then
           write(24,*) 'error, numline'
        end if

	call msd(ftrj,numline,nc,bi,nmaxl,dtsamp,xsamp,corr)
	
	do i=1,nc
          do j=1,3
           corrtr(ir,i,j)=corr(i,j)
	  end do
        end do

	end do !en do ir=1,nrea

	
	do m=1,nc
          p=int(1.2735**m)
          do l3=1,3
           avcorr(m,l3)=0.
	   do lr=1,nrea
             if (l3.eq.1) then
              write(5,*) sngl(p*dtsamp),lr,corrtr(lr,m,1)
             end if
             avcorr(m,l3)=avcorr(m,l3)+corrtr(lr,m,l3)
           end do !lr
           avcorr(m,l3)=avcorr(m,l3)/nrea
          end do !l3
          write(4,*) sngl(p*dtsamp),avcorr(m,1),avcorr(m,2),avcorr(m,3)
	end do !m

	close(4)
	close(5)
	
	end    !FIN DEL PROGRAM
	



      
c-------Subrutina MSD

	
	SUBROUTINE msd(fname,numline,nc,bi,nmaxl,dtsamp,x,corr)
	CHARACTER*21 fname
	INTEGER*8 nc,nmaxl,numline
	REAL*8 bi
	REAL*4 dtsamp
	INTEGER*4 k

        INTEGER j,i3,countx,i,h
	REAL*8 x(nmaxl)
        REAL*8 corrtmp(3),corr(nc,3),f(nmaxl)
	
	
        open(3,file=fname,status='old')
	  do i=1,numline
            read(3,*) x(i)
	  end do
	close(3)
		
	
	do j=1,nc
	  k=int(1.2735**j)
            do i3=1,3
             corrtmp(i3)=0.
             corr(j,i3)=0.
            end do !i3
            countx=0
            do i=1,numline-k
       	     countx=countx+1
             corrtmp(1)=corrtmp(1)+(x(i+k)-x(i))*(x(i+k)-x(i))  !MSD 
             corrtmp(2)=corrtmp(2)+x(i+k)*x(i)   !CORRELACION POSICION
             corrtmp(3)=corrtmp(3)+f(i+k)*f(i)	 !CORRELACION FUERZA
            end do 
            do i3=1,3	 
             corrtmp(i3)=corrtmp(i3)/countx
             corr(j,i3)=corr(j,i3)+corrtmp(i3)  
            end do         
	end do !j

      RETURN
      END



c-------------Funcion numero random diustribucion gaussiana gasdev
      
      FUNCTION gasdev(idum)
      INTEGER idum
      DOUBLE PRECISION gasdev
CU    USES ran1
      INTEGER iset
      DOUBLE PRECISION fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.d0*ran1(idum)-1.d0
        v2=2.d0*ran1(idum)-1.d0
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=sqrt(-2.d0*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END


c--------------------Funcion random

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END

c----------funcion logaritmo de funcion gamma

      FUNCTION gammln(xx)
      DOUBLE PRECISION gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
