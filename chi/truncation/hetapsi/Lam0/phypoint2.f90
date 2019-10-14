subroutine phypoint2(Nflow,yflow,rho0,mPion,mSigma,mf,Vall)
!calculate the physical point rho0

  implicit none
  real(8) hc
  parameter(hc=197.33)
  real(8) T,mu !temperature and chemical potential
  integer Nflow
  real(8) yflow(Nflow)
  real(8) rho0,rho
  real(8) mPion2,mPion,mSigma2,mSigma,mf
  real(8) Vall !total effective potential
  integer N_str(4) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  real(8) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(8) h
  real(8) Zphi,Zpsi
  real(8) c,kappa
  integer nn
  parameter(nn=1)
  real(8) x(nn)
  logical check
  external gapEq

  common /strucFun/ N_str
  common /gapPara/ lam1,lam2,lam3,lam4,lam5,lam6,lam7,c,kappa
  common /Tmu/ T,mu


  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

  lam1=yflow(1)
  lam2=yflow(2)
  lam3=yflow(3)
  lam4=yflow(4)
  lam5=yflow(5)
  lam6=yflow(6)
  lam7=yflow(7)
  lam0=yflow(Nv+1)
  h=yflow((Nv+1)+1)
  Zphi=yflow((Nv+1)+(Nh+1)+1)
  Zpsi=yflow((Nv+1)+(Nh+1)+2)
  c=yflow((Nv+1)+(Nh+1)+Nz+1)
  kappa=yflow((Nv+1)+(Nh+1)+Nz+2)

  if(T<200./hc)then
    x(1)=kappa
  else
    x(1)=kappa*1.1
  end if
  call newt(x, nn, check, gapEq)
  rho0=x(1)

  rho=rho0

  mPion2=lam1 + lam2*(-kappa + rho) + (lam3*(-kappa + rho)**2)/2. + &
        (lam4*(-kappa + rho)**3)/6. + (lam5*(-kappa + rho)**4)/24. + &
        (lam6*(-kappa + rho)**5)/120. + (lam7*(-kappa + rho)**6)/720.

  mPion=sqrt(mPion2)

  mSigma2=lam1 + lam2*(-kappa + rho) + (lam3*(-kappa + rho)**2)/2. + &
        (lam4*(-kappa + rho)**3)/6. + (lam5*(-kappa + rho)**4)/24. + &
        (lam6*(-kappa + rho)**5)/120. + (lam7*(-kappa + rho)**6)/720. + &
        2*rho*(lam2 + lam3*(-kappa + rho) + (lam4*(-kappa + rho)**2)/2. + &
        (lam5*(-kappa + rho)**3)/6. + (lam6*(-kappa + rho)**4)/24. + &
        (lam7*(-kappa + rho)**5)/120.)

  mSigma=sqrt(mSigma2)

  mf=h*sqrt(rho/2.)

  Vall=lam0 - Sqrt(2.)*c*Sqrt(rho) + lam1*(-kappa + rho) + &
       (lam2*(-kappa + rho)**2)/2. + (lam3*(-kappa + rho)**3)/6. + &
       (lam4*(-kappa + rho)**4)/24. + (lam5*(-kappa + rho)**5)/120. + &
       (lam6*(-kappa + rho)**6)/720. + (lam7*(-kappa + rho)**7)/5040.

end

