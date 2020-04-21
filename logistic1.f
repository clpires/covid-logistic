      parameter (nmax=1000,nparmax=10)
      character*120 file1
      real z(nmax), gcost(nparmax),hcost(nparmax,nparmax),par(nparmax)
      real hcostinv(nparmax,nparmax)
      complex eval(nparmax),evec(nparmax,nparmax)
      real evalr(nparmax),evecr(nparmax,nparmax)
      real pari(nparmax),parf(nparmax)
      real z1(nmax),z2(nmax),z3(nmax)
      
      common n1,npar,scalet,t0,scalea,niter,cof,nini,z

      
ccc  n1: number of available data in days      
      read*, n1
cccc ns: step to compute de first guess of the logistic function cc      
      read*, ns
ccccc scalet:  scale of time in the logistic function 
      read*, scalet
cccc scalea: scale for exp(r*ti) in the logistic         
      read*, scalea
ccccc time origin   or step at which time=t=0       
      read*, t0
ccccc zscale: scale of the carryimg capacity of the logistic      
      read*, zscale
ccccc initial instant included in the cost funtion = MSE      
      read*, nini
cccccc file with data      
      read*, file1
ccccc cof = coeffient that multilies the hessian matrix in the Newton method  
      read*, cof
ccc   niter = numebr of iterations
      read*, niter
cccc  n1min  = smallest value of te = last day included in the cost function
      read*, n1min
cccc  n1a  = largest value of te = last day included in the cost function
      read*, n1a
ccccccc

     
      open(10,file=file1)

cccc readinfg data ccc      
      do i=1,n1
         read(10,*) k,it,ia,id
cccc  it= total of infected, ia=total of actives ; id = total of deaths
         z(i)=float(ia)
c         z(i)=float(id)
         
ccc re scaling of the number cc
         z(i)=z(i)/zscale
      enddo
      
      do kkkk=n1min,n1a
         n1=kkkk
      
         npar=3
         dts=float(ns)/scalet
         
         c1=z(3*ns+1)/z(2*ns+1)
         c2=z(2*ns+1)/z(ns+1)
         
         beta=(c1*c2-1)/(c1*(c2-1))
         gama=(c1-1)/(c1*(c2-1))
         
         b=(beta-sqrt(beta**2-4*gama))/2.
         
         par(3)=-log(b)/dts
         par(2)=(c1-1)/(b**2-c1*b**3)
         par(1)=z(ns+1)*(1+b*par(2))
         
         xl=par(1)
         xa=par(2)/scalea
         xk=par(3)
         
         tin=log(par(2))/xk
         tind=tin*scalet
      
         f1=xlogi1(dts,xl,xa,xk)
         f2=xlogi1(2*dts,xl,xa,xk)
         f3=xlogi1(3*dts,xl,xa,xk)
         
         par(2)=par(2)/scalea

ccccc first guess of parameters cc

      do i=1,npar
         pari(i)=par(i)
      enddo

cccc  determination of the gradient and hessian matrix of the cost function to fit a logistic function ccc

      call scovid1(par,cost,gcost,hcost)

cccc  iteration 
     
      call scoviditer1(pari,parf,costi,costf,gcost,hcost)

cccc initail guess of parameters 
      
      xl1=pari(1)
      xa1=pari(2)
      xk1=pari(3)

cccccc final values of parameters cc
      
      xl2=parf(1)
      xa2=parf(2)
      xk2=parf(3)
            
      do i=1,n1
         t=(float(i-1)+t0)/scalet
         f1=xlogi1(t,xl1,xa1,xk1)*zscale
         f2=xlogi1(t,xl2,xa2,xk2)*zscale
         
ccc   write logistic for the fitting logistic using data until kkkk
         write(kkkk-1,*) t*scalet,f2
ccccc write observations          
         if(kkkk.eq.n1min)  write(15,*) t*scalet,z(i)*zscale
         
         z1(i)=xlogi1(t,xl2,xa2,xk2)
         z2(i)=log(xlogi1(t,xl2,xa2,xk2))
         z3(i)=log(z(i))
         
      enddo
      ntot=n1-nini+1
      call cor(ntot,z1(nini),z(nini),corfit)
      call cor(ntot,z2(nini),z3(nini),corfit1)

c      print*, corfit,corfit1
      
      zz=z(n1)
      th1=(parf(1)*zscale-100)/(parf(1)*zscale)
      th2=100./(parf(1)*zscale-100)
      tmax=-log(1./(scalea*parf(2))*th2)/parf(3)
      ff=xlogi1(tmax,xl2,xa2,xk2)*zscale
      n1max=int(tmax*scalet)+2

      do i=n1+1,n1max
         t=(float(i-1)+t0)/scalet
         f2=xlogi1(t,xl2,xa2,xk2)*zscale

cccc write logistic fit up to assimptotic value minus 100
         
         write(kkkk-1,*) t*scalet,f2
      enddo
      
      xl=parf(1)*zscale
      ti=log(parf(2)*scalea)/parf(3)*scalet
      xk=parf(3)/scalet
      
ccccc  writing parameters cccc
      
      
      write(91,*) kkkk-1, xl
      write(92,*) kkkk-1, ti
      write(93,*) kkkk-1, xk

      th=1.
      costf1=0.
      do i=nini,n1
         costf1=costf1+(z(i)*th)**2
      enddo

cccc write relative error of the fitting ccc
      
      write(94,*) kkkk-1,  sqrt(costf/costf1)
ccccc write correlations between (obs and logistic) and between logs of (obs and logistic) 
      
      write(95,*) kkkk-1,  corfit,corfit1
      
cccccccc end of the cicle of fittings ccc
      enddo
ccccccccccccccccccccccccc
      
      stop

      
      end

cccc subroutine to compute the gradient and hessian matrix of the costfunction to fit a the logistic function cccc
      
      subroutine scovid1(par,cost,gcost,hcost)
      parameter (nmax=1000,nparmax=10)
      real z(nmax), gcost(nparmax),hcost(nparmax,nparmax),par(nparmax)
      real pt(nmax)
      
      common n1,npar,scalet,t0,scalea,niter,cof,nini,z

      cost=0.
      xl=par(1)
      xa=par(2)
      xk=par(3)
      npar=3
      xtot=float(n1-nini+1)
      csig=0.1
      zvar=0.01
      
      do i=nini,n1
         t=(float(i-1)+t0)/scalet
         ft=xlogi1(t,xl,xa,xk)
         pt(i)=1.
         cost=cost+(ft-z(i))**2*pt(i)

      enddo
      
      do k1=1,npar
         gcost(k1)=0.
         do k2=1,npar
            hcost(k1,k2)=0.
         enddo
      enddo
      
      do i=nini,n1
         t=(float(i-1)+t0)/scalet
         ft=xlogi1(t,xl,xa,xk)
         ht=ft/xl
         zt=z(i)
         ekt=exp(-xk*t)*scalea
         
         gcost(1)=gcost(1)+2*(ft-zt)*ht*pt(i)
         gcost(2)=gcost(2)-2*xl*(ft-zt)*ht**2*ekt*pt(i)
         gcost(3)=gcost(3)+2*xl*xa*(ft-zt)*t*ht**2*ekt*pt(i)
         
         hcost(1,1)=hcost(1,1)+2*ht**2*pt(i)
         hcost(1,2)=hcost(1,2)-2*ht**2*ekt*(2*ft-zt)*pt(i)
         hcost(1,3)=hcost(1,3)+2*ht**2*xa*t*ekt*(2*ft-zt)*pt(i)
         hcost(2,2)=hcost(2,2)+2*xl*ht**3*ekt**2*(3*ft-2*zt)*pt(i)

         fzk=(ft-zt)*(2*t*ht**3*ekt**2*xa-t*ht**2*ekt)+ht**4*xl*t*xa*ekt**2
         
         hcost(2,3)=hcost(2,3)-2*xl*fzk*pt(i)
         hcost(3,3)=hcost(3,3)+2*xl*xa*t*fzk*pt(i)
                  
      enddo

      hcost(2,1)=hcost(1,2)
      hcost(3,1)=hcost(1,3)
      hcost(3,2)=hcost(2,3)
      
      
      return
      end
      
ccccc comupte the logistic for certain parameters ccc

      real function xlogi1(t,xl,xa,xk)
      common n1,npar,scalet,t0,scalea,niter,cof,nini,z


      xlogi1=xl/(1+scalea*xa*exp(-xk*t))

      
      return
      end

ccccc  iteration subroutine ccc

      
      subroutine  scoviditer1(pari,parf,costi,costf,gcost,hcost)
      parameter (nmax=1000,nparmax=10)

      real par(nparmax),gcost(nparmax),hcost(nparmax,nparmax)
      real hcostinv(nparmax,nparmax)
      real pari(nparmax),parf(nparmax),parn(nparmax)
      
      common n1,npar,scalet,t0,scalea,niter,cof,nini,z
      
      do i=1,npar
         par(i)=pari(i)
      enddo
      call scovid1(par,cost,gcost,hcost)
      costi=cost
      print*, 0,cost,(par(i),i=1,npar)
      
      do k1=1,niter
         call scovid1(par,cost,gcost,hcost)

ccc calls the imsl subroutine linrg to compute the inverse matrix 
         
         call linrg(npar,hcost,nparmax,hcostinv,nparmax)

         do k=1,npar
            parn(k)=par(k)
            do kk=1,npar
               parn(k)=parn(k)-cof*hcostinv(k,kk)*gcost(kk)
            enddo
         enddo

         do k=1,npar
            par(k)=parn(k)
         enddo
c         print*, k1,cost,(par(i),i=1,npar)
         
      enddo
      
      do k=1,npar
         parf(k)=par(k)
      enddo
      call scovid1(par,cost,gcost,hcost)
      costf=cost
c      print*, 'costf  ',cost,(par(i),i=1,npar)
     
      return
      end

      
