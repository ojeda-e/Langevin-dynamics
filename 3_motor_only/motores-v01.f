
	program motor

	implicit none

c------variables motores
	integer*4 seed,xm,it
	integer*8 itm,Fs,v0,dx8,itmm
	real*4 dt,t,dtsmot,dts,tt
	real*8 Pstep,PstepA,Pdetach,PdetachA,rnd1,rnd2,ran1
	real*4 indatach,Fm
	character*8 namefm
	character*21 motortraj

	open(2,file='mot.ini')
	read (2,*) namefm
	read (2,*) dt
	read (2,*) dtsmot
	read (2,*) seed
	read (2,*) itmm
	read (2,*) Fm
	read (2,*) Fs
	read (2,*) v0
	read (2,*) dx8
	close(2)
	

	motortraj=''//namefm//'.dat'

	
	indatach=0
	open(6,file=motortraj)

	xm=0.d0
	t=0.d0

      it=0
      do while (t.lt.10.and.indatach.eq.0)
      it=it+1
	  t=it*dt
	  
		if(0.le.Fm.and.Fm.lt.Fs) then 	!asigna probabilidad
			Pstep=v0*(1-(Fm/Fs)**2)/dx8
            else
			if(Fm.ge.Fs) then 
				Pstep=0
                else 
				Pstep=v0/dx8
			end if
		end if
    
	  Pdetach=0.5*exp(Fm/4)				!prob de detach

	  PstepA=Pstep*dt					!evaluacion evento de paso // prob absoluta
	  PdetachA=Pdetach*dt				!evaluacion evento de detach // prob absoluta

	  rnd1=ran1(seed)					!numero aleatorio PdetachA

	    if(rnd1.lt.PdetachA) then		!si rnd1 es menor que pdetach abs se desatacha
	  	  indatach=1
		    else 	  
			rnd2=ran1(seed)				!numero aleatrio PstepA
			if (rnd2.lt.PstepA) then  	!si rnd2 es menor que pstep abs da un paso
	   	    xm=xm+dx8
			end if			
	    end if	
	  
  	   if(mod(it,1000).eq.0) then
		 write(6,*) t,xm
  	   end if

	  end do  !do-while

	 close(6)
	 end program


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
