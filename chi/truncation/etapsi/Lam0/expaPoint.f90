subroutine expaPoint(Nflow,yflow,kappa_UV)

  implicit none

  integer Nflow
  real(8) yflow(Nflow) !sigma is the fixed expansion point at UV
  integer N_str(4) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  real(8) kappa_UV, kappa_old, kappa_new
  real(8) T,mu
  real(8) k_UV,k_IR,t_UV,t_IR
  external derivs,rkqs
  real(8) eps_ode,h1,hmin !variables in subroutine odeint
  integer nok,nbad !variables in subroutine odeint
  INTEGER kmax,kount !variables in common block of subroutine odeint
  INTEGER KMAXX,NMAX
  PARAMETER (NMAX=50,KMAXX=200)
  real(8) dxsav,xp(KMAXX),yp(NMAX,KMAXX) !variables in common block of subroutine odeint
  real(8) rho0,mPion,mSigma,mf,Vall
  real(8) epsi_rho0,epsi
  real(8) rho0_gaug
  real(8) kappa
  logical stopp
  integer i
  integer imax  !maximal number of loops
  parameter(imax=600)
  real(8) rescal,delta_IR

  common /strucFun/ N_str
  common /Tmu/ T,mu
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /odeContr/ eps_ode,h1,hmin
  COMMON /path/ kmax,kount,dxsav,xp,yp

  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)


!  epsi_rho0=0.00001
  epsi_rho0=0.
  epsi=1.e-5
  rescal=5.

  kappa_old=kappa_UV
  kappa_new=kappa_UV

  i=0                    !start of loops
  stopp=.false.
  do while((.not. stopp).and.(i < imax))
    i=i+1

    kappa_old=kappa_new
    call initial(Nflow,yflow,kappa_old)
    call odeint(yflow,Nflow,t_UV,t_IR,eps_ode,h1,hmin,nok,nbad,derivs,rkqs)
    call phypoint(Nflow,yflow,rho0,mPion,mSigma,mf,Vall)

    kappa=yflow((Nv+1)+(Nh+1)+Nz+2)
!    rho0_gaug=rho0*(1.+epsi_rho0)
    rho0_gaug=rho0+epsi_rho0

    delta_IR=kappa-rho0_gaug
    if(abs(delta_IR)/kappa<epsi)then
      stopp=.true.
    else
      kappa_new=kappa_old-delta_IR/rescal
    end if


!    write(*,"('delta_IR=', f15.7)")delta_IR
!    write(*,"('abs(delta_IR)/kappa=', f15.7)")abs(delta_IR)/kappa

  end do

  kappa_UV=kappa_old

end



