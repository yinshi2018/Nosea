program QM

  implicit none

  real(8) pi,hc
  parameter(pi=3.141592653589793d0)
  parameter(hc=197.33)
  real(8) T,mu !temperature and chemical potential
  real(8) l,lb !polyakov loop
  real(8) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,c,kappa
  real(8) sigma_UV,kappa_UV_i,kappa_UV_i_mu,kappa_UV
  real(8) Vtotal0
  real(8) Ti,dT,mu_down,mu_up
  integer i,iTmax,j,jmumax,m,mm,i_mm,j1
  parameter(iTmax=299,jmumax=51)
  !parameter(iTmax=0,jmumax=0)
  real(8) pre_res(0:jmumax,0:iTmax),T_res(0:iTmax),mu_res(0:jmumax,0:iTmax),pre_com(jmumax)
  real(8) mu_bound(0:iTmax),mu_bound_low,mu_bound_high
  real(8) T_MeV,muB_MeV
  real(8) mui,muBi,muB
  integer mmax
  parameter(mmax=10)
!order of chebyshev polynomial
  real(8) dcd0mu(mmax),dcd1mu(mmax),dcd2mu(mmax),dcd3mu(mmax),dcd4mu(mmax),chi(mmax,0:iTmax)
  real(8) factorial
  real(8) fpi_res(0:jmumax,0:iTmax),mPion_res(0:jmumax,0:iTmax),mSigma_res(0:jmumax,0:iTmax),mf_res(0:jmumax,0:iTmax)
  integer iT,iv
  real(8) l_i,lb_i,l_i_mu,lb_i_mu
  real(8) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(8) nfl,nfl0,nfl1,lset
  real(8) chebev
  external chebev


  common /Tmu/ T,mu
  common /prefit/ pre_com
  common /iTiv/ iT,iv



  Ti=1./hc                 !initial temperature
  dT=1./hc                 !stepsize of temperature

  do i=0, iTmax
	T_res(i)=Ti+dT*real(i)    !fm**(-1)
  end do

  mu_bound_low=200./hc                !fm**(-1)
  mu_bound_high=200./hc                !fm**(-1)

  do i=0, iTmax
    if(T_res(i)<mu_bound_low)then
      mu_bound(i)=mu_bound_low
    else
      mu_bound(i)=(mu_bound_high-mu_bound_low)/(T_res(iTmax)-mu_bound_low)*(T_res(i)-mu_bound_low)+mu_bound_low
    endif
  end do

!  mui=0.0/hc
!  muBi=3.*mui
  muBi=0./hc

  do i=0, iTmax
    mu_down=-mu_bound(i)/2.
    mu_up=mu_bound(i)/2.

    do j=0, jmumax
  	  if(j==0) then
        mu_res(j,i)=0.
	  else
        mu_res(j,i)=(cos(pi*(j-0.5d0)/jmumax)*(0.5d0*(mu_up-mu_down))+0.5d0*(mu_up+mu_down))   !fm**(-1)
	  end if
    end do

  end do



  sigma_UV=100./hc
  kappa_UV_i=sigma_UV**2/2.

  l_i=1.e-10
  lb_i=1.e-10

  do j1=0, jmumax
!  do j1=1, 1

    if(j1==0)then
      j=j1
    else
      j=jmumax+1-j1
    end if



    if(j/=0)then
      kappa_UV_i=kappa_UV_i_mu
      l_i=l_i_mu
      lb_i=lb_i_mu
    end if

!    do i=0, iTmax
    iT=iTmax
    do i=0, iT
      iv=i
      T=T_res(i)
      T_MeV=T*hc


      muB=muBi+mu_res(j,i)
      mu=1./3.*muB
      muB_MeV=muB*hc


      call selfEQ(kappa_UV_i,l_i,lb_i,kappa_UV,l,lb,rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,c,kappa)
      kappa_UV_i=kappa_UV
      l_i=l
      lb_i=lb

      if(i==0)then
        kappa_UV_i_mu=kappa_UV
        l_i_mu=l
        lb_i_mu=lb
      end if

      if(i==0.and.j==0)then
        Vtotal0=Vtotal
      end if


      pre_res(j,i)=-(Vtotal-Vtotal0)
      fpi_res(j,i)=fpi
      mPion_res(j,i)=mPion
      mSigma_res(j,i)=mSigma
      mf_res(j,i)=mf

      if(j==0)then


        nfl=Fnf0(mf,T,l,lb) + lb*Fnf1(mf,T,l,lb) + l*Fnf2(mf,T,l,lb)
        lset=0.
        nfl0=Fnf0(mf,T,lset,lset) + lset*Fnf1(mf,T,lset,lset) + lset*Fnf2(mf,T,lset,lset)
        lset=1.
        nfl1=Fnf0(mf,T,lset,lset) + lset*Fnf1(mf,T,lset,lset) + lset*Fnf2(mf,T,lset,lset)



        open(unit=51,file='./buffer/TMeV.dat',position='append')
        write(51, "(e21.14)")T_MeV
        close(51)

!  goto 210

        open(unit=51,file='./buffer/l.dat',position='append')
        write(51, "(e21.14)")l
        close(51)

        open(unit=51,file='./buffer/lb.dat',position='append')
        write(51, "(e21.14)")lb
        close(51)

        open(unit=51,file='./buffer/kappaUV.dat',position='append')
        write(51, "(e21.14)")kappa_UV
        close(51)

        open(unit=51,file='./buffer/rho0.dat',position='append')
        write(51, "(e21.14)")rho0
        close(51)

        open(unit=51,file='./buffer/mPion.dat',position='append')
        write(51, "(e21.14)")mPion*hc
        close(51)

        open(unit=51,file='./buffer/mSigma.dat',position='append')
        write(51, "(e21.14)")mSigma*hc
        close(51)

        open(unit=51,file='./buffer/mf.dat',position='append')
        write(51, "(e21.14)")mf*hc
        close(51)

        open(unit=51,file='./buffer/Vtotal.dat',position='append')
        write(51, "(e21.14)")Vtotal
        close(51)

        open(unit=51,file='./buffer/fpi.dat',position='append')
        write(51, "(e21.14)")fpi*hc
        close(51)

        open(unit=51,file='./buffer/h.dat',position='append')
        write(51, "(e21.14)")h
        close(51)

        open(unit=51,file='./buffer/Zphi.dat',position='append')
        write(51, "(e21.14)")Zphi
        close(51)

        open(unit=51,file='./buffer/Zpsi.dat',position='append')
        write(51, "(e21.14)")Zpsi
        close(51)

        open(unit=51,file='./buffer/c.dat',position='append')
        write(51, "(e21.14)")c
        close(51)

        open(unit=51,file='./buffer/kappa.dat',position='append')
        write(51, "(e21.14)")kappa
        close(51)

        open(unit=51,file='./buffer/failL.dat',position='append')
        write(51, "(e21.14)")l*exp(2.*mf/T)
        close(51)

        open(unit=51,file='./buffer/nfRatiol0.dat',position='append')
        write(51, "(e21.14)")nfl/nfl0
        close(51)

        open(unit=51,file='./buffer/nfRatiol1.dat',position='append')
        write(51, "(e21.14)")nfl/nfl1
        close(51)

        write(*,"('nfl/nfl1=', (e20.9))")nfl/nfl1

        open(unit=51,file='./buffer/nfl.dat',position='append')
        write(51, "(e21.14)")nfl
        close(51)

        open(unit=51,file='./buffer/nfl0.dat',position='append')
        write(51, "(e21.14)")nfl0
        close(51)

        open(unit=51,file='./buffer/nfl1.dat',position='append')
        write(51, "(e21.14)")nfl1
        close(51)


210 continue

      end if

      write(*,"('j=', I4,  t25, 'i=', I4)")j, i
      write(*,"('muB_MeV=', f15.7, t25, 'T_MeV=', f15.7)")muB_MeV,T_MeV
      write(*,"('fpi=', f15.7, t25, 'mPion=', f15.7,'mf=',f15.7)")fpi*hc,mPion*hc,mf*hc
      !write(*,"('kappa_UV=', f15.7)")kappa_UV

    end do
!    stop

  end do

  goto 220

  open(unit=51,file='./buffer/fpi_res.dat')
  do j=0, jmumax
    do i=0, iTmax
      write(51, "(e21.14)")fpi_res(j,i)
    end do
  end do
  close(51)

  open(unit=51,file='./buffer/mPion_res.dat')
  do j=0, jmumax
    do i=0, iTmax
      write(51, "(e21.14)")mPion_res(j,i)
    end do
  end do
  close(51)

  open(unit=51,file='./buffer/mSigma_res.dat')
  do j=0, jmumax
    do i=0, iTmax
      write(51, "(e21.14)")mSigma_res(j,i)
    end do
  end do
  close(51)

  open(unit=51,file='./buffer/mf_res.dat')
  do j=0, jmumax
    do i=0, iTmax
      write(51, "(e21.14)")mf_res(j,i)
    end do
  end do
  close(51)

220 continue

  open(unit=51,file='./buffer/pre_res.dat')
  do j=0, jmumax
    do i=0, iTmax
      write(51, "(e21.14)")pre_res(j,i)
    end do
  end do
  close(51)


  do i=0, iTmax

    mu_down=-mu_bound(i)/2.
    mu_up=mu_bound(i)/2.

    T=T_res(i)

    do j=1, jmumax
      pre_com(j)=pre_res(j,i)
	end do
    call chebft(dcd0mu,jmumax,mmax)

    call chder(mu_down,mu_up,dcd0mu,dcd1mu,mmax)
    call chder(mu_down,mu_up,dcd1mu,dcd2mu,mmax)
    call chder(mu_down,mu_up,dcd2mu,dcd3mu,mmax)
    call chder(mu_down,mu_up,dcd3mu,dcd4mu,mmax)

    chi(1,i)=chebev(mu_down,mu_up,dcd0mu,mmax,0.d0)/T**4
    chi(2,i)=chebev(mu_down,mu_up,dcd1mu,mmax,0.d0)/T**3
    chi(3,i)=chebev(mu_down,mu_up,dcd2mu,mmax,0.d0)/T**2
    chi(4,i)=chebev(mu_down,mu_up,dcd3mu,mmax,0.d0)/T
    chi(5,i)=chebev(mu_down,mu_up,dcd4mu,mmax,0.d0)

  end do

  open(unit=51,file='./buffer/chi0.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(1,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi1.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(2,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi2.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(3,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi3.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(4,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi4.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(5,i)
  end do
  close(51)

  open(unit=51,file='./buffer/p.dat')
  do i=0, iTmax
    write(51, "(e21.14)")pre_res(0,i)*hc**4
  end do
  close(51)

  open(unit=51,file='./buffer/muBi.dat')
  write(51, "(e21.14)")muBi*hc
  close(51)




end






