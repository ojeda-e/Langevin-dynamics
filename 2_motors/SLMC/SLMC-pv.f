	
        program motorfext

        implicit none

c-------variables motores
        integer seed,ixmot,Nxmot
        integer Fs,v0,itmm,ixm,Nxm,contstep,instep
        real*8 t,dtstepm,itdesp,dtmot,dx8
        real*8 Pstep,PstepA,Pdetach,PdetachA,rnd1,rnd2,ran1,xm
        real*8 indatach,Fm,xm_old,r1,r2,argxm,tdesp
        character*6 namefm
        character*21 motpvtraj

c-------variables globales
        integer nrea
        integer i,it
        real*8 dt

c-------variable fuerza elastica
        real*4 felas
        real*8 kkelas,kelas,xnat

c-------parametros
        integer nmod,nc,nmaxl,nream,lextmax
        parameter(nmod=16)		!numero de modos
        parameter(nc=60)  		!numero de puntos lag time
        parameter(nmaxl=20000000) 	!numero maximo de lineas a leer tray		
        parameter(nream=1000)		!numero maximo de realizaciones
        parameter(lextmax=10)		!valor maximo de fuerza externa L
        
c-------estadistica
        integer nr
        real*8 trun,vprom,tau,xrun
        real*8 Rmedia,V,Vmedia,R,time
        character*22 fpotdata,finaldata
        
c-------variables fuerza externa
        integer NLext,iLext
        real*8 Lext,Lmax

c-------variables particulas virtuales
        real*8 alpha,b0,g0,nu0,c
        real*8 sqetdt(nmod),sqet0dt,sumeta
        real*8 x,x0,u(nmod),KT,et(nmod),nu(nmod),k(nmod),xdesp !array con tamaño
        real*8 et0,etw,Radorg,rnd,gammln,s,gasdev,x_old,rnd3
        integer iseedg

c-------variables trabajo y portencia
        real*8 p_run,f_prom,work,P,F,W,Pmedia,Fmedia,Wmedio


        open(2,file='cond-pv.ini')
        read(2,*) namefm
c-------estadistica
        read(2,*) nrea
c-------motores
        read(2,*) dt
        read(2,*) seed
        read(2,*) iseedg
        read(2,*) itmm
        read(2,*) Fs
        read(2,*) v0
        read(2,*) dx8
		read(2,*) tdesp
c-------organela
        read(2,*) KT
        read(2,*) etw
		read(2,*) Radorg
        read(2,*) alpha
        read(2,*) b0
        read(2,*) g0
        read(2,*) nu0
        read(2,*) c
c-------fuerza elastica
        read(2,*) kelas
        read(2,*) xnat
        read(2,*) Lmax
        read(2,*) NLext
        close(2)
        
        motpvtraj=''//namefm//'-xm-fm.dat'
        fpotdata=''//namefm//'pot-xm.dat'
        finaldata='tabla-'//namefm//'.dat'

c        open(3,file=motpvtraj)
        open(4,file=fpotdata)
        open(5,file=finaldata)


        itdesp=tdesp/dt

		et0=10*Radorg*etw/10**6*6*3.1416/10000

        sumeta=0
        sqet0dt=sqrt(2*KT/et0*dt)	   	
        do i=1,nmod
           nu(i)=nu0/b0**(i-1)   
           k(i)=g0/dexp(gammln(1.d0-alpha))*nu(i)**alpha*c
           et(i)=k(i)/nu(i)  
           sqetdt(i)=sqrt(2*KT/et(i)*dt)
           sumeta=sumeta+et(i)
        end do

        do iLext=0,NLext
        Lext=iLext*Lmax/NLext
c		write(*,*) Lext
       
        Rmedia=0.d0
        Vmedia=0.d0
        tau=0.d0	
        V=0.d0
        R=0.d0
        time=0.d0
        
            do nr=1,nrea
            it=0
            x0=0.d0
            x=0.d0
            t=0.d0
            indatach=0	
            vprom=0.d0
            trun=0.d0
            xrun=0.d0
            xm=0.d0
			argxm=ran1(seed)*220-110. !nint(2*r1*xnat/dx8)-nint(xnat/dx8)
            xm=xm+argxm
			xdesp=0.d0
			xm_old=xm
			x_old=x
			s=0.d0

                do while (t.lt.10.and.indatach.eq.0)
                 it=it+1
                 t=it*dt
                  
                 Fm=felas(x_old,xm_old,kelas,xnat)
c                  end if
c                 Borrador Fm=finterac(x,xm) calculo de fuerza entre organela y motor
c                 s=0.d0
                 

                 Pdetach=1.*exp(abs(Fm)/4.)				!prob de detach
                 PdetachA=Pdetach*dt				!evaluacion evento de detach // prob absoluta

                 rnd1=ran1(seed)					!numero aleatorio PdetachA

                    if(rnd1.lt.PdetachA.and.it.gt.itdesp) then		!si rnd1 es menor que pdetach abs se desatacha
                      indatach=1
                    else 	  
                        if(Fm.ge.0.and.Fm.lt.Fs) then 	!asigna probabilidad 
	
                        Pstep=v0*(1.-(Fm/Fs)**2)/dx8
                        else
                            if(Fm.ge.Fs) then 
                            Pstep=0
                            else 
                                Pstep=v0/dx8
                            end if
                        end if
                    PstepA=Pstep*dt					!evaluacion evento de paso // prob absoluta
                    rnd2=ran1(seed)				!numero aleatrio PstepA
c					write(*,*) xm,t
                        if (rnd2.lt.PstepA) then  	!si rnd2 es menor que pstep abs da un paso
c							do while(t.le.dt+dtmot)
							xm=xm+dx8
c							end do
                        end if			
                    end if	
c   					write(*,*) 0,xm
c                 if(mod(it,100).eq.0) then
c                   write(3,*) t,sngl(fm),sngl(xm)
c                 end if
			

c--------------Dinamica organela

c                 s=Fm
c              Borrador x=x+dt/gamma*(Fm)+ruido  (dinamica organela)
                 s=Fm
                 
                    do i=1,nmod
                    s=s-k(i)*(x-u(i))
                    rnd=gasdev(seed)
                    u(i)=u(i)+sqetdt(i)*rnd+dt*nu(i)*(x-u(i))
                    end do

                         
                 rnd3=gasdev(seed)
                 x=x+sqet0dt*rnd3+s*dt/et0-Lext*dt/et0
				 
				 	if (t.lt.tdesp) then
						xdesp=x
c						write(*,*) t,xdesp
					end if

c                    if(mod(it,1000).eq.0) then
c                    write(4,*) t,x,s,xm
c                    end if


                 x_old=x
                 xm_old=xm

                end do  !do-while
			
            xrun=x-xdesp
c			write(*,*) xrun,x,xdesp
            trun=t-tdesp
            vprom=xrun/trun
            p_run=Fm*xrun/trun
            work=Fm*xrun
           
            W=W+work
            P=P+p_run
            R=R+xrun
            time=time+trun
            V=V+vprom
            F=F+Fm

            end do !end do realizaciones 
                Wmedio=W/nrea
                Rmedia=R/nrea
                Vmedia=V/nrea
                tau=time/nrea
                Pmedia=P/nrea
                Fmedia=F/nrea

           write(4,*) sngl(Lext),sngl(Fmedia),sngl(Pmedia),sngl(Wmedio)
           write(5,*) sngl(Lext),sngl(Vmedia),sngl(Rmedia),sngl(tau)
        end do ! end do fuerza externa iLext


c        close(3)
        close(4)
        close(5)
        end program

c--------------Funcion fuerza elastica

      FUNCTION felas(x,y,kkelas,length)

      real*8 x,y,kkelas,dyx,sgn,length

       felas=0.d0
       dyx=abs(y-x)
       if (dyx.gt.length) then
          felas=kkelas*(dyx-length)*(y-x)/dyx
       end if

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
      
