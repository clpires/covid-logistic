      parameter (nmax=1000,nparmax=10)
      character*120 file1
      real z(nmax), gcost(nparmax),hcost(nparmax,nparmax),par(nparmax)
      real hcostinv(nparmax,nparmax)
      complex eval(nparmax),evec(nparmax,nparmax)
      real evalr(nparmax),evecr(nparmax,nparmax)
      real pari(nparmax),parf(nparmax)
      real z1(0:nmax),z2(0:nmax),z3(0:nmax)
      real zl1(0:nmax),zl2(0:nmax),zl3(0:nmax)

      real dz1(0:nmax),dz2(0:nmax),dz3(0:nmax)
      real dzl1(0:nmax),dzl2(0:nmax),dzl3(0:nmax)
c      real z1l(0:nmax),z2l(0:nmax),z3l(0:nmax)
      real z1l(-nmax:nmax),z2l(-nmax:nmax),z3l(-nmax:nmax)
      
      common n1,npar,scalet,t0,scalea,niter,cof,nini,z

      read*, n1
      read*, ns
      read*, scalet
      read*, scalea
      read*, t0
      read*, zscale
      read*, nini
      read*, file1
      read*, xl,ti,xk
      read*, xld,tid,xkd
      read*, lagid
      open(10,file=file1)

      
cccc  setting data at time t=-1  : n=0
      z1(0)=0
      z2(0)=0
      z3(0)=0
      zl1(0)=0
      zl2(0)=0
      zl3(0)=0
ccc read data cc      
      
      do i=1,n1
         read(10,*) k,it,ia,id

         z1(i)=float(it)
         z2(i)=float(ia)
         z3(i)=float(id)
         
         if(z3(i).gt.0) then
            zl3(i)=log(z3(i))
         else
            zl3(i)=0.
         endif
         
         if(z1(i).gt.0) then
            zl1(i)=log(z1(i))
         else
            zl1(i)=0.
         endif

         if(z2(i).gt.0) then
            zl2(i)=log(z2(i))
         else
            zl2(i)=0.
         endif

      enddo

      do i=1,n1

ccc relative daily variations cc         
         dzl2(i)=(z2(i)-z2(i-1))/z2(i-1)
         dzl3(i)=(z3(i)-z3(i-1))/z3(i-1)
         dzl1(i)=(z1(i)-z1(i-1))/z1(i-1)

ccc daily variations ccc
         
         dz2(i)=z2(i)-z2(i-1)
         dz3(i)=z3(i)-z3(i-1)
         dz1(i)=z1(i)-z1(i-1)
                  
      enddo
      
cccc italia  active 
c      xl=92850.7656
c      ti=36.0649452
c      xk=0.185417086
c       xa=exp(xk*ti)

cccc Spain active 
c      xl= 82373.3828
c      ti=39.1659431
c      xk=0.239941284 
c      xa=exp(xk*ti)
      
      nti=int(ti)+1
      tend=ti-log(5.)/xk
      tend=ti-log(1.)/xk
      xa=exp(xk*ti)
      
      tend=ti
      nend=int(tend)+1
      nti=int(ti)+1

cccc computation of the logistic for the active cccc      
      
      do i=-nmax,nmax
         t=float(i-1)
         z2l(i)= xlogi1d(t,xl,xa,xk)
      enddo
      
ccc   Italia dead

c      xld=16807.8828 
c      tid=39.9597511
c      xkd=0.195912778
c      xad=exp(xkd*tid)

cccc Spain dead 
c      xl= 10959.7158
c      ti=41.3946114
c      xk=0.297599077
c     xa=exp(xk*ti)

      
      xad=exp(xkd*tid)

cccc  computation of the logistic for the dead cccc
      
      do i=-nmax,nmax
         t=float(i-1)
         z3l(i)= xlogi1d(t,xld,xad,xkd)
      enddo

ccccc computation of the error dif 
      
      
      do lag=0,24
          call errcovid1(nti-nini+1,z2l(nini),z3l(nini),lag,dif,difm2)
          call errcovid1(nti-nini+1,dz2(nini),dz3(nini),lag,dif,difm1)
          print*, lag,sqrt(difm1),sqrt(difm2)
          write(22,*) lag,sqrt(difm1),sqrt(difm2)
      enddo

      stop

      end

cccc subroutine computing the SSE (dif) and the mean SSE (difm)      
      
      subroutine errcovid1(n,a,b,lag,dif,difm)

      real a(n),b(n)
      n1=n-lag
      err=0.
      do i=1,n-lag
         err=err+(a(i)-b(i+lag))**2
      enddo
      dif=err
      
      err=err/float(n-lag)
      difm=err


      return
      end

ccc logistic 
      
      real function xlogi1(t,xl,xa,xk)
      xlogi1=xl/(1+xa*exp(-xk*t))
      return
      end

ccc one step of the logistic 
      
      real function xlogi1d(t,xl,xa,xk)

      xlogi1d=xlogi1(t,xl,xa,xk)-xlogi1(t-1,xl,xa,xk)
      
      return
      end
      
