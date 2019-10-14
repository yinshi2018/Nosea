subroutine initial(Nflow,yflow,kappa)
!make the initialization

  implicit none
  integer Nflow
  real(8) yflow(Nflow) !sigma is the fixed expansion point at UV
  integer N_str(4) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  real(8) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(8) h,dr1dh,dr2dh,dr3dh,dr4dh,dr5dh
  real(8) Zphi,Zpsi
  real(8) c,kappa
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) lambda,nu

  common /strucFun/ N_str

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

  lambda=1.
  nu=(96./hc)**2

  h=6.5
!expansion coefficients of Yukawa coupling

  c=1.7e-3*(1000./hc)**3  !explicit chiral symmetry breaking term, in unit of fm**(-3)


  lam0=(kappa**2*lambda)/2.+nu*kappa
  lam1=kappa*lambda+nu
  lam2=lambda
  lam3=0.
  lam4=0.
  lam5=0.
  lam6=0.
  lam7=0.
!expansion coefficients of effective potential V

  Zphi=1. !meson wave function renormalization
  Zpsi=1. !quark wave function renormalization
  dr1dh=0.
  dr2dh=0.
  dr3dh=0.
  dr4dh=0.
  dr5dh=0.

  yflow(1)=lam1
  yflow(2)=lam2
  yflow(3)=lam3
  yflow(4)=lam4
  yflow(5)=lam5
  yflow(6)=lam6
  yflow(7)=lam7
  yflow(Nv+1)=lam0
  yflow((Nv+1)+1)=h
  yflow((Nv+1)+2)=dr1dh
  yflow((Nv+1)+3)=dr2dh
  yflow((Nv+1)+4)=dr3dh
  yflow((Nv+1)+5)=dr4dh
  yflow((Nv+1)+6)=dr5dh
  yflow((Nv+1)+(Nh+1)+1)=Zphi
  yflow((Nv+1)+(Nh+1)+2)=Zpsi
  yflow((Nv+1)+(Nh+1)+Nz+1)=c
  yflow((Nv+1)+(Nh+1)+Nz+2)=kappa

end





