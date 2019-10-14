subroutine derivs(x,y,dydx)
!Calculating the right hand side of differential equations

  implicit none

  integer NMAX !maximal number of differential equations
  parameter(NMAX=50)
  real(8) x,y(NMAX),dydx(NMAX)
  integer N_str(4) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  real(8) k ! IR cutoff in flow equations
  real(8) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(8) h
  real(8) Zphi,Zpsi
  real(8) c,kappa
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) etaphi,etapsi
!meson and quark anomanous dimension
  real(8) Nc,Nf
  parameter(Nc=3.,Nf=2.)
  real(8) v3
  parameter(v3=1./(2.*pi**2))
  real(8) rho !phi_a**2/2
  real(8) Fnb,Fnf0,Fnf1,Fnf2
  external Fnb,Fnf0,Fnf1,Fnf2
  real(8) zb,zf !distinguish the transverse and longituidanl wave function renormalization
  real(8) p0,p0c !temporal compontent of external momentum
  real(8) T,mu,Lambda2
  real(8) mu0
  real(8) l,lb !polyakov loop
  real(8) k_UV,k_IR,t_UV,t_IR
  real(8) mp2,ms2,mf2,mp2d1rho,mp2d2rho,mp2d3rho,mp2d4rho,mp2d5rho,ms2d1rho,ms2d2rho,ms2d3rho,ms2d4rho,ms2d5rho,mf2d1rho
  real(8) nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x
  real(8) nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,nbSigma,nbd1xSigma,nbd2xSigma, &
          nbd3xSigma,nbd4xSigma,nbd5xSigma
  real(8) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(8) nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x
  real(8) xff,xfa
  real(8) f2a,f3a,b2b2PS,b2f1aPion,b2f1aPionI,b2f1aSigma,b2f1aSigmaI,b1f2Pion,b1f2PionI,b1f2Sigma,b1f2SigmaI
  complex(8) b2f1aPionC,b2f1aSigmaC,b1f2PionC,b1f2SigmaC
  real(8) dr0dtV,dr1dtV,dr2dtV,dr3dtV,dr4dtV,dr5dtV
  real(8) dlam0dt,dlam1dt,dlam2dt,dlam3dt,dlam4dt,dlam5dt,dlam6dt,dlam7dt
  real(8) L11Pion,L11Sigma
  real(8) dth,dhdt
  real(8) dZphidt,dZpsidt,dcdt,dkappadt
  real(8) l_com,lb_com


  common /strucFun/ N_str
  common /Tmu/ T,mu
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /polyakov_com/ l_com,lb_com
  common /Lam/ Lambda2


  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

  k=k_UV*exp(x)
  lam1=y(1)
  lam2=y(2)
  lam3=y(3)
  lam4=y(4)
  lam5=y(5)
  lam6=y(6)
  lam7=y(7)
  lam0=y(Nv+1)
  h=y((Nv+1)+1)
  Zphi=y((Nv+1)+(Nh+1)+1)
  Zpsi=y((Nv+1)+(Nh+1)+2)
  c=y((Nv+1)+(Nh+1)+Nz+1)
  kappa=y((Nv+1)+(Nh+1)+Nz+2)


  rho=kappa !calculations are performed at expansion point kappa

  zb=1.
  zf=1.


  mu0=0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mass and their derivatives
  mp2=lam1/k**2

  ms2=(lam1 + 2*lam2*rho)/k**2

  mf2=(h**2*rho)/(k**2*Nf)

  mp2d1rho=lam2/k**2

  mp2d2rho=lam3/k**2

  mp2d3rho=lam4/k**2

  mp2d4rho=lam5/k**2

  mp2d5rho=lam6/k**2

  ms2d1rho=(3*lam2 + 2*lam3*rho)/k**2

  ms2d2rho=(5*lam3 + 2*lam4*rho)/k**2

  ms2d3rho=(7*lam4 + 2*lam5*rho)/k**2

  ms2d4rho=(9*lam5 + 2*lam6*rho)/k**2

  ms2d5rho=(11*lam6 + 2*lam7*rho)/k**2

  mf2d1rho=h**2/(k**2*Nf)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nb=Fnb((k*Sqrt(1 + mp2))/Sqrt(zb),T)
  call nbdx(nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x)
  nbPion=nbd0x
  nbd1xPion=nbd1x
  nbd2xPion=nbd2x
  nbd3xPion=nbd3x
  nbd4xPion=nbd4x
  nbd5xPion=nbd5x

  nb=Fnb((k*Sqrt(1 + ms2))/Sqrt(zb),T)
  call nbdx(nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x)
  nbSigma=nbd0x
  nbd1xSigma=nbd1x
  nbd2xSigma=nbd2x
  nbd3xSigma=nbd3x
  nbd4xSigma=nbd4x
  nbd5xSigma=nbd5x


  xff=-mu + (k*Sqrt(1 + mf2))/zf
  xfa=mu + (k*Sqrt(1 + mf2))/zf

  l=l_com
  lb=lb_com

  nf0=Fnf0(xff,T,l,lb)
  nf1=Fnf1(xff,T,l,lb)
  nf2=Fnf2(xff,T,l,lb)
  call nfdx(l,lb,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  nff=nfd0x
  nfd1xf=nfd1x
  nfd2xf=nfd2x
  nfd3xf=nfd3x
  nfd4xf=nfd4x
  nfd5xf=nfd5x



  nf0=Fnf0(xfa,T,lb,l)
  nf1=Fnf1(xfa,T,lb,l)
  nf2=Fnf2(xfa,T,lb,l)
  call nfdx(lb,l,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  nfa=nfd0x
  nfd1xa=nfd1x
  nfd2xa=nfd2x
  nfd3xa=nfd3x
  nfd4xa=nfd4x
  nfd5xa=nfd5x


  f2a=(k*Sqrt(1 + mf2)*nfd1xa + k*Sqrt(1 + mf2)*nfd1xf + zf - nfa*zf - nff*zf)/ &
  (4.*(1 + mf2)**1.5*zf**2)

  f3a=(-((k**2*(1 + mf2)*(nfd2xa + nfd2xf))/zf) +                               &
    3*(k*Sqrt(1 + mf2)*nfd1xa + k*Sqrt(1 + mf2)*nfd1xf + zf - nfa*zf -          &
       nff*zf))/(16.*(1 + mf2)**2.5*zf**2)

  b2b2PS=-((-2*(1 + 2*nbPion))/(Sqrt(1 + mp2)*(mp2 - ms2)**3) -                 &
     (1 + 2*nbPion)/(2.*(1 + mp2)**1.5*(mp2 - ms2)**2) -                        &
     (2*(1 + 2*nbSigma))/(Sqrt(1 + ms2)*(-mp2 + ms2)**3) -                      &
     (1 + 2*nbSigma)/(2.*(1 + ms2)**1.5*(-mp2 + ms2)**2) +                      &
     (k*nbd1xPion)/((1 + mp2)*(mp2 - ms2)**2*Sqrt(zb)) +                        &
     (k*nbd1xSigma)/((1 + ms2)*(-mp2 + ms2)**2*Sqrt(zb)))/(2.*Sqrt(zb))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  p0=pi*T
  p0c=sqrt((pi*(T-1./hc)*exp(-k/T/5.)+pi*1./hc)**2+k**2)

  b2f1aPionC=-(k**2*(-((k*nbPion*(-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mp2))/Sqrt(zb)))/  &
          ((1 + mp2)*((-mu0 + Complex(0,1)*p0 -                                 &
                  (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)      &
**2)) - (k*nbd1xPion)/                                                          &
        (2.*(1 + mp2)*((-mu0 + Complex(0,1)*p0 -                                &
               (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)) +      &
       (nbPion*Sqrt(zb))/                                                       &
        (2.*(1 + mp2)**1.5*((-mu0 + Complex(0,1)*p0 -                           &
               (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)) -      &
       (k*(-mu0 + Complex(0,1)*p0c - (k*Sqrt(1 + mp2))/Sqrt(zb)))/              &
        ((1 + mp2)*((-mu0 + Complex(0,1)*p0c -                                  &
                (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**2     &
) + Sqrt(zb)/(2.*(1 + mp2)**1.5*                                                &
          ((-mu0 + Complex(0,1)*p0c - (k*Sqrt(1 + mp2))/Sqrt(zb))**2 -          &
            (k**2*(1 + mf2))/zf**2)) +                                          &
       (k*nbPion*(-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mp2))/Sqrt(zb)))/        &
        ((1 + mp2)*((-mu0 + Complex(0,1)*p0 +                                   &
                (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**2     &
) - (k*nbd1xPion)/                                                              &
        (2.*(1 + mp2)*((-mu0 + Complex(0,1)*p0 +                                &
               (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)) +      &
       (nbPion*Sqrt(zb))/                                                       &
        (2.*(1 + mp2)**1.5*((-mu0 + Complex(0,1)*p0 +                           &
               (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)) -      &
       (k**2*zf)/                                                               &
        (Sqrt(1 + mf2)*zb*(-((k**2*(1 + mp2))/zb) +                             &
             (-mu0 + Complex(0,1)*p0c + (k*Sqrt(1 + mf2))/zf)**2)**2) +         &
       ((k**2*nfa*zf)/                                                          &
           (Sqrt(1 + mf2)*zb*(-((k**2*(1 + mp2))/zb) +                          &
                (-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)**2) +       &
          (k**2*nff*zf)/                                                        &
           (Sqrt(1 + mf2)*zb*(-((k**2*(1 + mp2))/zb) +                          &
                (mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)**2))/2.      &
+ ((k**2*nff*zf)/                                                               &
           (Sqrt(1 + mf2)*zb*(-((k**2*(1 + mp2))/zb) +                          &
                (-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)**2) +       &
          (k**2*nfa*zf)/                                                        &
           (Sqrt(1 + mf2)*zb*(-((k**2*(1 + mp2))/zb) +                          &
                (mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)**2))/2.))/   &
  (2.*zb*zf**2)


  b2f1aPion=real(b2f1aPionC)
  b2f1aPionI=aimag(b2f1aPionC)

  b2f1aSigmaC=-(k**2*(-((k*nbSigma*(-mu0 + Complex(0,1)*p0 -                    &
              (k*Sqrt(1 + ms2))/Sqrt(zb)))/                                     &
          ((1 + ms2)*((-mu0 + Complex(0,1)*p0 -                                 &
                  (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)      &
**2)) - (k*nbd1xSigma)/                                                         &
        (2.*(1 + ms2)*((-mu0 + Complex(0,1)*p0 -                                &
               (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)) +      &
       (nbSigma*Sqrt(zb))/                                                      &
        (2.*(1 + ms2)**1.5*((-mu0 + Complex(0,1)*p0 -                           &
               (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)) -      &
       (k*(-mu0 + Complex(0,1)*p0c - (k*Sqrt(1 + ms2))/Sqrt(zb)))/              &
        ((1 + ms2)*((-mu0 + Complex(0,1)*p0c -                                  &
                (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**2     &
) + Sqrt(zb)/(2.*(1 + ms2)**1.5*                                                &
          ((-mu0 + Complex(0,1)*p0c - (k*Sqrt(1 + ms2))/Sqrt(zb))**2 -          &
            (k**2*(1 + mf2))/zf**2)) +                                          &
       (k*nbSigma*(-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + ms2))/Sqrt(zb)))/       &
        ((1 + ms2)*((-mu0 + Complex(0,1)*p0 +                                   &
                (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**2     &
) - (k*nbd1xSigma)/                                                             &
        (2.*(1 + ms2)*((-mu0 + Complex(0,1)*p0 +                                &
               (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)) +      &
       (nbSigma*Sqrt(zb))/                                                      &
        (2.*(1 + ms2)**1.5*((-mu0 + Complex(0,1)*p0 +                           &
               (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)) -      &
       (k**2*zf)/                                                               &
        (Sqrt(1 + mf2)*zb*(-((k**2*(1 + ms2))/zb) +                             &
             (-mu0 + Complex(0,1)*p0c + (k*Sqrt(1 + mf2))/zf)**2)**2) +         &
       ((k**2*nfa*zf)/                                                          &
           (Sqrt(1 + mf2)*zb*(-((k**2*(1 + ms2))/zb) +                          &
                (-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)**2) +       &
          (k**2*nff*zf)/                                                        &
           (Sqrt(1 + mf2)*zb*(-((k**2*(1 + ms2))/zb) +                          &
                (mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)**2))/2.      &
+ ((k**2*nff*zf)/                                                               &
           (Sqrt(1 + mf2)*zb*(-((k**2*(1 + ms2))/zb) +                          &
                (-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)**2) +       &
          (k**2*nfa*zf)/                                                        &
           (Sqrt(1 + mf2)*zb*(-((k**2*(1 + ms2))/zb) +                          &
                (mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)**2))/2.))/   &
  (2.*zb*zf**2)


  b2f1aSigma=real(b2f1aSigmaC)
  b2f1aSigmaI=aimag(b2f1aSigmaC)

  b1f2PionC=-(k**2*((k*(-mu0 + Complex(0,1)*p0c + (k*Sqrt(1 + mf2))/zf))/       &
        ((1 + mf2)*(-((k**2*(1 + mp2))/zb) +                                    &
             (-mu0 + Complex(0,1)*p0c + (k*Sqrt(1 + mf2))/zf)**2)**2) -         &
       (k**2*nbPion*Sqrt(zb))/                                                  &
        (Sqrt(1 + mp2)*((-mu0 + Complex(0,1)*p0 -                               &
                (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**      &
           2*zf**2) - (k**2*Sqrt(zb))/                                          &
        (Sqrt(1 + mp2)*((-mu0 + Complex(0,1)*p0c -                              &
                (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**      &
           2*zf**2) - (k**2*nbPion*Sqrt(zb))/                                   &
        (Sqrt(1 + mp2)*((-mu0 + Complex(0,1)*p0 +                               &
                (k*Sqrt(1 + mp2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**      &
           2*zf**2) + zf/                                                       &
        (2.*(1 + mf2)**1.5*(-((k**2*(1 + mp2))/zb) +                            &
            (-mu0 + Complex(0,1)*p0c + (k*Sqrt(1 + mf2))/zf)**2)) +             &
       ((k*nfd1xa)/                                                             &
           (2.*(1 + mf2)*(-((k**2*(1 + mp2))/zb) +                              &
               (-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)) +           &
          (k*nfd1xf)/                                                           &
           (2.*(1 + mf2)*(-((k**2*(1 + mp2))/zb) +                              &
               (mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)) +            &
          (k*nfa*(-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf))/              &
           ((1 + mf2)*(-((k**2*(1 + mp2))/zb) +                                 &
                (-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)**2) +       &
          (k*nff*(mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf))/               &
           ((1 + mf2)*(-((k**2*(1 + mp2))/zb) +                                 &
                (mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)**2) -        &
          (nfa*zf)/                                                             &
           (2.*(1 + mf2)**1.5*                                                  &
             (-((k**2*(1 + mp2))/zb) +                                          &
               (-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)) -           &
          (nff*zf)/                                                             &
           (2.*(1 + mf2)**1.5*                                                  &
             (-((k**2*(1 + mp2))/zb) +                                          &
               (mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)))/2. +        &
       ((k*nfd1xf)/                                                             &
           (2.*(1 + mf2)*(-((k**2*(1 + mp2))/zb) +                              &
               (-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)) +           &
          (k*nfd1xa)/                                                           &
           (2.*(1 + mf2)*(-((k**2*(1 + mp2))/zb) +                              &
               (mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)) -            &
          (k*nff*(-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf))/              &
           ((1 + mf2)*(-((k**2*(1 + mp2))/zb) +                                 &
                (-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)**2) -       &
          (k*nfa*(mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf))/               &
           ((1 + mf2)*(-((k**2*(1 + mp2))/zb) +                                 &
                (mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)**2) -        &
          (nff*zf)/                                                             &
           (2.*(1 + mf2)**1.5*                                                  &
             (-((k**2*(1 + mp2))/zb) +                                          &
               (-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)) -           &
          (nfa*zf)/                                                             &
           (2.*(1 + mf2)**1.5*(-((k**2*(1 + mp2))/zb) +                         &
               (mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)))/2.))/       &
  (2.*zb*zf**2)


  b1f2Pion=real(b1f2PionC)
  b1f2PionI=aimag(b1f2PionC)

  b1f2SigmaC=-(k**2*((k*(-mu0 + Complex(0,1)*p0c + (k*Sqrt(1 + mf2))/zf))/      &
        ((1 + mf2)*(-((k**2*(1 + ms2))/zb) +                                    &
             (-mu0 + Complex(0,1)*p0c + (k*Sqrt(1 + mf2))/zf)**2)**2) -         &
       (k**2*nbSigma*Sqrt(zb))/                                                 &
        (Sqrt(1 + ms2)*((-mu0 + Complex(0,1)*p0 -                               &
                (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**      &
           2*zf**2) - (k**2*Sqrt(zb))/                                          &
        (Sqrt(1 + ms2)*((-mu0 + Complex(0,1)*p0c -                              &
                (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**      &
           2*zf**2) - (k**2*nbSigma*Sqrt(zb))/                                  &
        (Sqrt(1 + ms2)*((-mu0 + Complex(0,1)*p0 +                               &
                (k*Sqrt(1 + ms2))/Sqrt(zb))**2 - (k**2*(1 + mf2))/zf**2)**      &
           2*zf**2) + zf/                                                       &
        (2.*(1 + mf2)**1.5*(-((k**2*(1 + ms2))/zb) +                            &
            (-mu0 + Complex(0,1)*p0c + (k*Sqrt(1 + mf2))/zf)**2)) +             &
       ((k*nfd1xa)/                                                             &
           (2.*(1 + mf2)*(-((k**2*(1 + ms2))/zb) +                              &
               (-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)) +           &
          (k*nfd1xf)/                                                           &
           (2.*(1 + mf2)*(-((k**2*(1 + ms2))/zb) +                              &
               (mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)) +            &
          (k*nfa*(-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf))/              &
           ((1 + mf2)*(-((k**2*(1 + ms2))/zb) +                                 &
                (-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)**2) +       &
          (k*nff*(mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf))/               &
           ((1 + mf2)*(-((k**2*(1 + ms2))/zb) +                                 &
                (mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)**2) -        &
          (nfa*zf)/                                                             &
           (2.*(1 + mf2)**1.5*                                                  &
             (-((k**2*(1 + ms2))/zb) +                                          &
               (-mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)) -           &
          (nff*zf)/                                                             &
           (2.*(1 + mf2)**1.5*                                                  &
             (-((k**2*(1 + ms2))/zb) +                                          &
               (mu0 + Complex(0,1)*p0 - (k*Sqrt(1 + mf2))/zf)**2)))/2. +        &
       ((k*nfd1xf)/                                                             &
           (2.*(1 + mf2)*(-((k**2*(1 + ms2))/zb) +                              &
               (-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)) +           &
          (k*nfd1xa)/                                                           &
           (2.*(1 + mf2)*(-((k**2*(1 + ms2))/zb) +                              &
               (mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)) -            &
          (k*nff*(-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf))/              &
           ((1 + mf2)*(-((k**2*(1 + ms2))/zb) +                                 &
                (-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)**2) -       &
          (k*nfa*(mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf))/               &
           ((1 + mf2)*(-((k**2*(1 + ms2))/zb) +                                 &
                (mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)**2) -        &
          (nff*zf)/                                                             &
           (2.*(1 + mf2)**1.5*                                                  &
             (-((k**2*(1 + ms2))/zb) +                                          &
               (-mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)) -           &
          (nfa*zf)/                                                             &
           (2.*(1 + mf2)**1.5*(-((k**2*(1 + ms2))/zb) +                         &
               (mu0 + Complex(0,1)*p0 + (k*Sqrt(1 + mf2))/zf)**2)))/2.))/       &
  (2.*zb*zf**2)


  b1f2Sigma=real(b1f2SigmaC)
  b1f2SigmaI=aimag(b1f2SigmaC)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  etaphi=-((f2a*h**2*Nc*v3 - (8*f3a*h**2*Nc*v3)/3. -                            &
!      (4*b2b2PS*lam2**2*rho*v3)/(3.*k**2) +                                     &
!      ((b2f1aSigma*h**2 + b2f1aPion*h**2*(-1 + Nf**2))*v3*                      &
!         ((-2*f2a*h**2*Nc*v3)/3. + (4*f3a*h**2*Nc*v3)/3.))/(3.*Nf))/            &
!    (1 - ((b2f1aSigma*h**2 + b2f1aPion*h**2*(-1 + Nf**2))*v3*                   &
!         ((-2*f2a*h**2*Nc*v3)/3. + (4*f3a*h**2*Nc*v3)/3.))/(12.*Nf)))


!  etapsi=((-(b2f1aPion*h**2) + b2f1aSigma*h**2 + b2f1aPion*h**2*Nf**2)*v3*      &
!    (12*k**2 + 3*f2a*h**2*k**2*Nc*v3 - 8*f3a*h**2*k**2*Nc*v3 -                  &
!      4*b2b2PS*lam2**2*rho*v3))/                                                &
!  (2.*k**2*(18*Nf - b2f1aPion*f2a*h**4*Nc*v3**2 +                               &
!      b2f1aSigma*f2a*h**4*Nc*v3**2 + 2*b2f1aPion*f3a*h**4*Nc*v3**2 -            &
!      2*b2f1aSigma*f3a*h**4*Nc*v3**2 + b2f1aPion*f2a*h**4*Nc*Nf**2*v3**2 -      &
!      2*b2f1aPion*f3a*h**4*Nc*Nf**2*v3**2))
  etapsi=((4 - etaphi)*v3*(h**2*(-1 + Nf**2)*b2f1aPion + h**2*b2f1aSigma))/(12.*Nf)
  etaphi=0.
!  etapsi=0.



  dZphidt=-etaphi*Zphi
  dZpsidt=-etapsi*Zpsi

  dcdt=(1./2.)*etaphi*c
  dkappadt=-etaphi*kappa



  if(k>Lambda2)then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!nosea part
  dr0dtV=(k**4*v3*((2*(1 - etaphi/5.)*(0.5 + nbSigma))/(3.*Sqrt(1 + ms2)*Sqrt(zb)) +   &
      (2*(1 - etaphi/5.)*(0.5 + nbPion)*(-1 + Nf**2))/                          &
       (3.*Sqrt(1 + mp2)*Sqrt(zb)) -                                            &
      (4*(1 - etapsi/4.)*Nc*Nf*( - nfa - nff))/(3.*Sqrt(1 + mf2)*zf)))/2.

  dr1dtv=(k**4*v3*(((1 - etaphi/5.)*k*ms2d1rho*nbd1xSigma)/                     &
         (3.*(1 + ms2)*zb) + ((1 - etaphi/5.)*k*mp2d1rho*nbd1xPion*             &
           (-1 + Nf**2))/(3.*(1 + mp2)*zb) -                                    &
        ((1 - etaphi/5.)*ms2d1rho*(0.5 + nbSigma))/                             &
         (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                         &
        ((1 - etaphi/5.)*mp2d1rho*(0.5 + nbPion)*(-1 + Nf**2))/                 &
         (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                         &
        (2*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*(-nfa - nff))/                        &
         (3.*(1 + mf2)**1.5*zf) -                                               &
        (4*(1 - etapsi/4.)*Nc*Nf*                                               &
           (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                        &
             (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                       &
         (3.*Sqrt(1 + mf2)*zf)))/2.


  dr2dtv=(k**4*v3*(((1 - etaphi/5.)*k**2*ms2d1rho**2*nbd2xSigma)/               &
         (6.*(1 + ms2)**1.5*zb**1.5) +                                          &
        ((1 - etaphi/5.)*k**2*mp2d1rho**2*nbd2xPion*(-1 + Nf**2))/              &
         (6.*(1 + mp2)**1.5*zb**1.5) -                                          &
        ((1 - etaphi/5.)*k*ms2d1rho**2*nbd1xSigma)/(2.*(1 + ms2)**2*zb) +       &
        ((1 - etaphi/5.)*k*ms2d2rho*nbd1xSigma)/(3.*(1 + ms2)*zb) -             &
        ((1 - etaphi/5.)*k*mp2d1rho**2*nbd1xPion*(-1 + Nf**2))/                 &
         (2.*(1 + mp2)**2*zb) +                                                 &
        ((1 - etaphi/5.)*k*mp2d2rho*nbd1xPion*(-1 + Nf**2))/                    &
         (3.*(1 + mp2)*zb) + ((1 - etaphi/5.)*ms2d1rho**2*                      &
           (0.5 + nbSigma))/(2.*(1 + ms2)**2.5*Sqrt(zb)) -                      &
        ((1 - etaphi/5.)*ms2d2rho*(0.5 + nbSigma))/                             &
         (3.*(1 + ms2)**1.5*Sqrt(zb)) +                                         &
        ((1 - etaphi/5.)*mp2d1rho**2*(0.5 + nbPion)*(-1 + Nf**2))/              &
         (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                         &
        ((1 - etaphi/5.)*mp2d2rho*(0.5 + nbPion)*(-1 + Nf**2))/                 &
         (3.*(1 + mp2)**1.5*Sqrt(zb)) -                                         &
        ((1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*(-nfa - nff))/                       &
         ((1 + mf2)**2.5*zf) +                                                  &
        (4*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                                      &
           (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                        &
             (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                       &
         (3.*(1 + mf2)**1.5*zf) -                                               &
        (4*(1 - etapsi/4.)*Nc*Nf*                                               &
           (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                   &
             (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                   &
             (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                    &
             (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                   &
         (3.*Sqrt(1 + mf2)*zf)))/2.



  dr3dtv=(k**4*v3*(((1 - etaphi/5.)*k**3*ms2d1rho**3*nbd3xSigma)/               &
         (12.*(1 + ms2)**2*zb**2) +                                             &
        ((1 - etaphi/5.)*k**3*mp2d1rho**3*nbd3xPion*(-1 + Nf**2))/              &
         (12.*(1 + mp2)**2*zb**2) -                                             &
        ((1 - etaphi/5.)*k**2*ms2d1rho**3*nbd2xSigma)/                          &
         (2.*(1 + ms2)**2.5*zb**1.5) +                                          &
        ((1 - etaphi/5.)*k**2*ms2d1rho*ms2d2rho*nbd2xSigma)/                    &
         (2.*(1 + ms2)**1.5*zb**1.5) -                                          &
        ((1 - etaphi/5.)*k**2*mp2d1rho**3*nbd2xPion*(-1 + Nf**2))/              &
         (2.*(1 + mp2)**2.5*zb**1.5) +                                          &
        ((1 - etaphi/5.)*k**2*mp2d1rho*mp2d2rho*nbd2xPion*(-1 + Nf**2))/        &
         (2.*(1 + mp2)**1.5*zb**1.5) +                                          &
        (5*(1 - etaphi/5.)*k*ms2d1rho**3*nbd1xSigma)/                           &
         (4.*(1 + ms2)**3*zb) -                                                 &
        (3*(1 - etaphi/5.)*k*ms2d1rho*ms2d2rho*nbd1xSigma)/                     &
         (2.*(1 + ms2)**2*zb) +                                                 &
        ((1 - etaphi/5.)*k*ms2d3rho*nbd1xSigma)/(3.*(1 + ms2)*zb) +             &
        (5*(1 - etaphi/5.)*k*mp2d1rho**3*nbd1xPion*(-1 + Nf**2))/               &
         (4.*(1 + mp2)**3*zb) -                                                 &
        (3*(1 - etaphi/5.)*k*mp2d1rho*mp2d2rho*nbd1xPion*(-1 + Nf**2))/         &
         (2.*(1 + mp2)**2*zb) +                                                 &
        ((1 - etaphi/5.)*k*mp2d3rho*nbd1xPion*(-1 + Nf**2))/                    &
         (3.*(1 + mp2)*zb) - (5*(1 - etaphi/5.)*ms2d1rho**3*                    &
           (0.5 + nbSigma))/(4.*(1 + ms2)**3.5*Sqrt(zb)) +                      &
        (3*(1 - etaphi/5.)*ms2d1rho*ms2d2rho*(0.5 + nbSigma))/                  &
         (2.*(1 + ms2)**2.5*Sqrt(zb)) -                                         &
        ((1 - etaphi/5.)*ms2d3rho*(0.5 + nbSigma))/                             &
         (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                         &
        (5*(1 - etaphi/5.)*mp2d1rho**3*(0.5 + nbPion)*(-1 + Nf**2))/            &
         (4.*(1 + mp2)**3.5*Sqrt(zb)) +                                         &
        (3*(1 - etaphi/5.)*mp2d1rho*mp2d2rho*(0.5 + nbPion)*(-1 + Nf**2))/      &
         (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                         &
        ((1 - etaphi/5.)*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/                 &
         (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                         &
        (5*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*(-nfa - nff))/                     &
         (2.*(1 + mf2)**3.5*zf) -                                               &
        (3*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*                                   &
           (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                        &
             (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                       &
         ((1 + mf2)**2.5*zf) +                                                  &
        (2*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                                      &
           (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                   &
             (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                   &
             (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                    &
             (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                   &
         ((1 + mf2)**1.5*zf) -                                                  &
        (4*(1 - etapsi/4.)*Nc*Nf*                                               &
           (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -              &
             (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +              &
             (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +              &
             (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -              &
             (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                  &
             (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                 &
         (3.*Sqrt(1 + mf2)*zf)))/2.



  dr4dtv=(k**4*v3*(((1 - etaphi/5.)*k**4*ms2d1rho**4*nbd4xSigma)/               &
         (24.*(1 + ms2)**2.5*zb**2.5) +                                         &
        ((1 - etaphi/5.)*k**4*mp2d1rho**4*nbd4xPion*(-1 + Nf**2))/              &
         (24.*(1 + mp2)**2.5*zb**2.5) -                                         &
        (5*(1 - etaphi/5.)*k**3*ms2d1rho**4*nbd3xSigma)/                        &
         (12.*(1 + ms2)**3*zb**2) +                                             &
        ((1 - etaphi/5.)*k**3*ms2d1rho**2*ms2d2rho*nbd3xSigma)/                 &
         (2.*(1 + ms2)**2*zb**2) -                                              &
        (5*(1 - etaphi/5.)*k**3*mp2d1rho**4*nbd3xPion*(-1 + Nf**2))/            &
         (12.*(1 + mp2)**3*zb**2) +                                             &
        ((1 - etaphi/5.)*k**3*mp2d1rho**2*mp2d2rho*nbd3xPion*                   &
           (-1 + Nf**2))/(2.*(1 + mp2)**2*zb**2) +                              &
        (15*(1 - etaphi/5.)*k**2*ms2d1rho**4*nbd2xSigma)/                       &
         (8.*(1 + ms2)**3.5*zb**1.5) -                                          &
        (3*(1 - etaphi/5.)*k**2*ms2d1rho**2*ms2d2rho*nbd2xSigma)/               &
         ((1 + ms2)**2.5*zb**1.5) +                                             &
        ((1 - etaphi/5.)*k**2*ms2d2rho**2*nbd2xSigma)/                          &
         (2.*(1 + ms2)**1.5*zb**1.5) +                                          &
        (2*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d3rho*nbd2xSigma)/                  &
         (3.*(1 + ms2)**1.5*zb**1.5) +                                          &
        (15*(1 - etaphi/5.)*k**2*mp2d1rho**4*nbd2xPion*(-1 + Nf**2))/           &
         (8.*(1 + mp2)**3.5*zb**1.5) -                                          &
        (3*(1 - etaphi/5.)*k**2*mp2d1rho**2*mp2d2rho*nbd2xPion*                 &
           (-1 + Nf**2))/((1 + mp2)**2.5*zb**1.5) +                             &
        ((1 - etaphi/5.)*k**2*mp2d2rho**2*nbd2xPion*(-1 + Nf**2))/              &
         (2.*(1 + mp2)**1.5*zb**1.5) +                                          &
        (2*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d3rho*nbd2xPion*(-1 + Nf**2))/      &
         (3.*(1 + mp2)**1.5*zb**1.5) -                                          &
        (35*(1 - etaphi/5.)*k*ms2d1rho**4*nbd1xSigma)/                          &
         (8.*(1 + ms2)**4*zb) +                                                 &
        (15*(1 - etaphi/5.)*k*ms2d1rho**2*ms2d2rho*nbd1xSigma)/                 &
         (2.*(1 + ms2)**3*zb) -                                                 &
        (3*(1 - etaphi/5.)*k*ms2d2rho**2*nbd1xSigma)/                           &
         (2.*(1 + ms2)**2*zb) -                                                 &
        (2*(1 - etaphi/5.)*k*ms2d1rho*ms2d3rho*nbd1xSigma)/                     &
         ((1 + ms2)**2*zb) + ((1 - etaphi/5.)*k*ms2d4rho*nbd1xSigma)/           &
         (3.*(1 + ms2)*zb) - (35*(1 - etaphi/5.)*k*mp2d1rho**4*nbd1xPion*       &
           (-1 + Nf**2))/(8.*(1 + mp2)**4*zb) +                                 &
        (15*(1 - etaphi/5.)*k*mp2d1rho**2*mp2d2rho*nbd1xPion*                   &
           (-1 + Nf**2))/(2.*(1 + mp2)**3*zb) -                                 &
        (3*(1 - etaphi/5.)*k*mp2d2rho**2*nbd1xPion*(-1 + Nf**2))/               &
         (2.*(1 + mp2)**2*zb) -                                                 &
        (2*(1 - etaphi/5.)*k*mp2d1rho*mp2d3rho*nbd1xPion*(-1 + Nf**2))/         &
         ((1 + mp2)**2*zb) + ((1 - etaphi/5.)*k*mp2d4rho*nbd1xPion*             &
           (-1 + Nf**2))/(3.*(1 + mp2)*zb) +                                    &
        (35*(1 - etaphi/5.)*ms2d1rho**4*(0.5 + nbSigma))/                       &
         (8.*(1 + ms2)**4.5*Sqrt(zb)) -                                         &
        (15*(1 - etaphi/5.)*ms2d1rho**2*ms2d2rho*(0.5 + nbSigma))/              &
         (2.*(1 + ms2)**3.5*Sqrt(zb)) +                                         &
        (3*(1 - etaphi/5.)*ms2d2rho**2*(0.5 + nbSigma))/                        &
         (2.*(1 + ms2)**2.5*Sqrt(zb)) +                                         &
        (2*(1 - etaphi/5.)*ms2d1rho*ms2d3rho*(0.5 + nbSigma))/                  &
         ((1 + ms2)**2.5*Sqrt(zb)) -                                            &
        ((1 - etaphi/5.)*ms2d4rho*(0.5 + nbSigma))/                             &
         (3.*(1 + ms2)**1.5*Sqrt(zb)) +                                         &
        (35*(1 - etaphi/5.)*mp2d1rho**4*(0.5 + nbPion)*(-1 + Nf**2))/           &
         (8.*(1 + mp2)**4.5*Sqrt(zb)) -                                         &
        (15*(1 - etaphi/5.)*mp2d1rho**2*mp2d2rho*(0.5 + nbPion)*                &
           (-1 + Nf**2))/(2.*(1 + mp2)**3.5*Sqrt(zb)) +                         &
        (3*(1 - etaphi/5.)*mp2d2rho**2*(0.5 + nbPion)*(-1 + Nf**2))/            &
         (2.*(1 + mp2)**2.5*Sqrt(zb)) +                                         &
        (2*(1 - etaphi/5.)*mp2d1rho*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/      &
         ((1 + mp2)**2.5*Sqrt(zb)) -                                            &
        ((1 - etaphi/5.)*mp2d4rho*(0.5 + nbPion)*(-1 + Nf**2))/                 &
         (3.*(1 + mp2)**1.5*Sqrt(zb)) -                                         &
        (35*(1 - etapsi/4.)*mf2d1rho**4*Nc*Nf*(-nfa - nff))/                    &
         (4.*(1 + mf2)**4.5*zf) +                                               &
        (10*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*                                  &
           (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                        &
             (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                       &
         ((1 + mf2)**3.5*zf) -                                                  &
        (6*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*                                   &
           (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                   &
             (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                   &
             (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                    &
             (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                   &
         ((1 + mf2)**2.5*zf) +                                                  &
        (8*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                                      &
           (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -              &
             (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +              &
             (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +              &
             (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -              &
             (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                  &
             (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                 &
         (3.*(1 + mf2)**1.5*zf) -                                               &
        (4*(1 - etapsi/4.)*Nc*Nf*                                               &
           (-(k**4*mf2d1rho**4*nfd4xa)/(16.*(1 + mf2)**2*zf**4) -               &
             (k**4*mf2d1rho**4*nfd4xf)/(16.*(1 + mf2)**2*zf**4) +               &
             (3*k**3*mf2d1rho**4*nfd3xa)/(8.*(1 + mf2)**2.5*zf**3) +            &
             (3*k**3*mf2d1rho**4*nfd3xf)/(8.*(1 + mf2)**2.5*zf**3) -            &
             (15*k**2*mf2d1rho**4*nfd2xa)/(16.*(1 + mf2)**3*zf**2) -            &
             (15*k**2*mf2d1rho**4*nfd2xf)/(16.*(1 + mf2)**3*zf**2) +            &
             (15*k*mf2d1rho**4*nfd1xa)/(16.*(1 + mf2)**3.5*zf) +                &
             (15*k*mf2d1rho**4*nfd1xf)/(16.*(1 + mf2)**3.5*zf)))/               &
         (3.*Sqrt(1 + mf2)*zf)))/2.



  dr5dtv=(k**4*v3*(((1 - etaphi/5.)*k**5*ms2d1rho**5*nbd5xSigma)/               &
         (48.*(1 + ms2)**3*zb**3) +                                             &
        ((1 - etaphi/5.)*k**5*mp2d1rho**5*nbd5xPion*(-1 + Nf**2))/              &
         (48.*(1 + mp2)**3*zb**3) -                                             &
        (5*(1 - etaphi/5.)*k**4*ms2d1rho**5*nbd4xSigma)/                        &
         (16.*(1 + ms2)**3.5*zb**2.5) +                                         &
        (5*(1 - etaphi/5.)*k**4*ms2d1rho**3*ms2d2rho*nbd4xSigma)/               &
         (12.*(1 + ms2)**2.5*zb**2.5) -                                         &
        (5*(1 - etaphi/5.)*k**4*mp2d1rho**5*nbd4xPion*(-1 + Nf**2))/            &
         (16.*(1 + mp2)**3.5*zb**2.5) +                                         &
        (5*(1 - etaphi/5.)*k**4*mp2d1rho**3*mp2d2rho*nbd4xPion*                 &
           (-1 + Nf**2))/(12.*(1 + mp2)**2.5*zb**2.5) +                         &
        (35*(1 - etaphi/5.)*k**3*ms2d1rho**5*nbd3xSigma)/                       &
         (16.*(1 + ms2)**4*zb**2) -                                             &
        (25*(1 - etaphi/5.)*k**3*ms2d1rho**3*ms2d2rho*nbd3xSigma)/              &
         (6.*(1 + ms2)**3*zb**2) +                                              &
        (5*(1 - etaphi/5.)*k**3*ms2d1rho*ms2d2rho**2*nbd3xSigma)/               &
         (4.*(1 + ms2)**2*zb**2) +                                              &
        (5*(1 - etaphi/5.)*k**3*ms2d1rho**2*ms2d3rho*nbd3xSigma)/               &
         (6.*(1 + ms2)**2*zb**2) +                                              &
        (35*(1 - etaphi/5.)*k**3*mp2d1rho**5*nbd3xPion*(-1 + Nf**2))/           &
         (16.*(1 + mp2)**4*zb**2) -                                             &
        (25*(1 - etaphi/5.)*k**3*mp2d1rho**3*mp2d2rho*nbd3xPion*                &
           (-1 + Nf**2))/(6.*(1 + mp2)**3*zb**2) +                              &
        (5*(1 - etaphi/5.)*k**3*mp2d1rho*mp2d2rho**2*nbd3xPion*                 &
           (-1 + Nf**2))/(4.*(1 + mp2)**2*zb**2) +                              &
        (5*(1 - etaphi/5.)*k**3*mp2d1rho**2*mp2d3rho*nbd3xPion*                 &
           (-1 + Nf**2))/(6.*(1 + mp2)**2*zb**2) -                              &
        (35*(1 - etaphi/5.)*k**2*ms2d1rho**5*nbd2xSigma)/                       &
         (4.*(1 + ms2)**4.5*zb**1.5) +                                          &
        (75*(1 - etaphi/5.)*k**2*ms2d1rho**3*ms2d2rho*nbd2xSigma)/              &
         (4.*(1 + ms2)**3.5*zb**1.5) -                                          &
        (15*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d2rho**2*nbd2xSigma)/              &
         (2.*(1 + ms2)**2.5*zb**1.5) -                                          &
        (5*(1 - etaphi/5.)*k**2*ms2d1rho**2*ms2d3rho*nbd2xSigma)/               &
         ((1 + ms2)**2.5*zb**1.5) +                                             &
        (5*(1 - etaphi/5.)*k**2*ms2d2rho*ms2d3rho*nbd2xSigma)/                  &
         (3.*(1 + ms2)**1.5*zb**1.5) +                                          &
        (5*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d4rho*nbd2xSigma)/                  &
         (6.*(1 + ms2)**1.5*zb**1.5) -                                          &
        (35*(1 - etaphi/5.)*k**2*mp2d1rho**5*nbd2xPion*(-1 + Nf**2))/           &
         (4.*(1 + mp2)**4.5*zb**1.5) +                                          &
        (75*(1 - etaphi/5.)*k**2*mp2d1rho**3*mp2d2rho*nbd2xPion*                &
           (-1 + Nf**2))/(4.*(1 + mp2)**3.5*zb**1.5) -                          &
        (15*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d2rho**2*nbd2xPion*                &
           (-1 + Nf**2))/(2.*(1 + mp2)**2.5*zb**1.5) -                          &
        (5*(1 - etaphi/5.)*k**2*mp2d1rho**2*mp2d3rho*nbd2xPion*                 &
           (-1 + Nf**2))/((1 + mp2)**2.5*zb**1.5) +                             &
        (5*(1 - etaphi/5.)*k**2*mp2d2rho*mp2d3rho*nbd2xPion*(-1 + Nf**2))/      &
         (3.*(1 + mp2)**1.5*zb**1.5) +                                          &
        (5*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d4rho*nbd2xPion*(-1 + Nf**2))/      &
         (6.*(1 + mp2)**1.5*zb**1.5) +                                          &
        (315*(1 - etaphi/5.)*k*ms2d1rho**5*nbd1xSigma)/                         &
         (16.*(1 + ms2)**5*zb) -                                                &
        (175*(1 - etaphi/5.)*k*ms2d1rho**3*ms2d2rho*nbd1xSigma)/                &
         (4.*(1 + ms2)**4*zb) +                                                 &
        (75*(1 - etaphi/5.)*k*ms2d1rho*ms2d2rho**2*nbd1xSigma)/                 &
         (4.*(1 + ms2)**3*zb) +                                                 &
        (25*(1 - etaphi/5.)*k*ms2d1rho**2*ms2d3rho*nbd1xSigma)/                 &
         (2.*(1 + ms2)**3*zb) -                                                 &
        (5*(1 - etaphi/5.)*k*ms2d2rho*ms2d3rho*nbd1xSigma)/                     &
         ((1 + ms2)**2*zb) - (5*(1 - etaphi/5.)*k*ms2d1rho*ms2d4rho*            &
           nbd1xSigma)/(2.*(1 + ms2)**2*zb) +                                   &
        ((1 - etaphi/5.)*k*ms2d5rho*nbd1xSigma)/(3.*(1 + ms2)*zb) +             &
        (315*(1 - etaphi/5.)*k*mp2d1rho**5*nbd1xPion*(-1 + Nf**2))/             &
         (16.*(1 + mp2)**5*zb) -                                                &
        (175*(1 - etaphi/5.)*k*mp2d1rho**3*mp2d2rho*nbd1xPion*                  &
           (-1 + Nf**2))/(4.*(1 + mp2)**4*zb) +                                 &
        (75*(1 - etaphi/5.)*k*mp2d1rho*mp2d2rho**2*nbd1xPion*(-1 + Nf**2))/     &
         (4.*(1 + mp2)**3*zb) +                                                 &
        (25*(1 - etaphi/5.)*k*mp2d1rho**2*mp2d3rho*nbd1xPion*(-1 + Nf**2))/     &
         (2.*(1 + mp2)**3*zb) -                                                 &
        (5*(1 - etaphi/5.)*k*mp2d2rho*mp2d3rho*nbd1xPion*(-1 + Nf**2))/         &
         ((1 + mp2)**2*zb) - (5*(1 - etaphi/5.)*k*mp2d1rho*mp2d4rho*            &
           nbd1xPion*(-1 + Nf**2))/(2.*(1 + mp2)**2*zb) +                       &
        ((1 - etaphi/5.)*k*mp2d5rho*nbd1xPion*(-1 + Nf**2))/                    &
         (3.*(1 + mp2)*zb) - (315*(1 - etaphi/5.)*ms2d1rho**5*                  &
           (0.5 + nbSigma))/(16.*(1 + ms2)**5.5*Sqrt(zb)) +                     &
        (175*(1 - etaphi/5.)*ms2d1rho**3*ms2d2rho*(0.5 + nbSigma))/             &
         (4.*(1 + ms2)**4.5*Sqrt(zb)) -                                         &
        (75*(1 - etaphi/5.)*ms2d1rho*ms2d2rho**2*(0.5 + nbSigma))/              &
         (4.*(1 + ms2)**3.5*Sqrt(zb)) -                                         &
        (25*(1 - etaphi/5.)*ms2d1rho**2*ms2d3rho*(0.5 + nbSigma))/              &
         (2.*(1 + ms2)**3.5*Sqrt(zb)) +                                         &
        (5*(1 - etaphi/5.)*ms2d2rho*ms2d3rho*(0.5 + nbSigma))/                  &
         ((1 + ms2)**2.5*Sqrt(zb)) +                                            &
        (5*(1 - etaphi/5.)*ms2d1rho*ms2d4rho*(0.5 + nbSigma))/                  &
         (2.*(1 + ms2)**2.5*Sqrt(zb)) -                                         &
        ((1 - etaphi/5.)*ms2d5rho*(0.5 + nbSigma))/                             &
         (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                         &
        (315*(1 - etaphi/5.)*mp2d1rho**5*(0.5 + nbPion)*(-1 + Nf**2))/          &
         (16.*(1 + mp2)**5.5*Sqrt(zb)) +                                        &
        (175*(1 - etaphi/5.)*mp2d1rho**3*mp2d2rho*(0.5 + nbPion)*               &
           (-1 + Nf**2))/(4.*(1 + mp2)**4.5*Sqrt(zb)) -                         &
        (75*(1 - etaphi/5.)*mp2d1rho*mp2d2rho**2*(0.5 + nbPion)*                &
           (-1 + Nf**2))/(4.*(1 + mp2)**3.5*Sqrt(zb)) -                         &
        (25*(1 - etaphi/5.)*mp2d1rho**2*mp2d3rho*(0.5 + nbPion)*                &
           (-1 + Nf**2))/(2.*(1 + mp2)**3.5*Sqrt(zb)) +                         &
        (5*(1 - etaphi/5.)*mp2d2rho*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/      &
         ((1 + mp2)**2.5*Sqrt(zb)) +                                            &
        (5*(1 - etaphi/5.)*mp2d1rho*mp2d4rho*(0.5 + nbPion)*(-1 + Nf**2))/      &
         (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                         &
        ((1 - etaphi/5.)*mp2d5rho*(0.5 + nbPion)*(-1 + Nf**2))/                 &
         (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                         &
        (315*(1 - etapsi/4.)*mf2d1rho**5*Nc*Nf*(-nfa - nff))/                   &
         (8.*(1 + mf2)**5.5*zf) -                                               &
        (175*(1 - etapsi/4.)*mf2d1rho**4*Nc*Nf*                                 &
           (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                        &
             (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                       &
         (4.*(1 + mf2)**4.5*zf) +                                               &
        (25*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*                                  &
           (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                   &
             (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                   &
             (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                    &
             (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                   &
         ((1 + mf2)**3.5*zf) -                                                  &
        (10*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*                                  &
           (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -              &
             (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +              &
             (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +              &
             (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -              &
             (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                  &
             (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                 &
         ((1 + mf2)**2.5*zf) +                                                  &
        (10*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                                     &
           (-(k**4*mf2d1rho**4*nfd4xa)/(16.*(1 + mf2)**2*zf**4) -               &
             (k**4*mf2d1rho**4*nfd4xf)/(16.*(1 + mf2)**2*zf**4) +               &
             (3*k**3*mf2d1rho**4*nfd3xa)/(8.*(1 + mf2)**2.5*zf**3) +            &
             (3*k**3*mf2d1rho**4*nfd3xf)/(8.*(1 + mf2)**2.5*zf**3) -            &
             (15*k**2*mf2d1rho**4*nfd2xa)/(16.*(1 + mf2)**3*zf**2) -            &
             (15*k**2*mf2d1rho**4*nfd2xf)/(16.*(1 + mf2)**3*zf**2) +            &
             (15*k*mf2d1rho**4*nfd1xa)/(16.*(1 + mf2)**3.5*zf) +                &
             (15*k*mf2d1rho**4*nfd1xf)/(16.*(1 + mf2)**3.5*zf)))/               &
         (3.*(1 + mf2)**1.5*zf) -                                               &
        (4*(1 - etapsi/4.)*Nc*Nf*                                               &
           (-(k**5*mf2d1rho**5*nfd5xa)/(32.*(1 + mf2)**2.5*zf**5) -             &
             (k**5*mf2d1rho**5*nfd5xf)/(32.*(1 + mf2)**2.5*zf**5) +             &
             (5*k**4*mf2d1rho**5*nfd4xa)/(16.*(1 + mf2)**3*zf**4) +             &
             (5*k**4*mf2d1rho**5*nfd4xf)/(16.*(1 + mf2)**3*zf**4) -             &
             (45*k**3*mf2d1rho**5*nfd3xa)/(32.*(1 + mf2)**3.5*zf**3) -          &
             (45*k**3*mf2d1rho**5*nfd3xf)/(32.*(1 + mf2)**3.5*zf**3) +          &
             (105*k**2*mf2d1rho**5*nfd2xa)/(32.*(1 + mf2)**4*zf**2) +           &
             (105*k**2*mf2d1rho**5*nfd2xf)/(32.*(1 + mf2)**4*zf**2) -           &
             (105*k*mf2d1rho**5*nfd1xa)/(32.*(1 + mf2)**4.5*zf) -               &
             (105*k*mf2d1rho**5*nfd1xf)/(32.*(1 + mf2)**4.5*zf)))/              &
         (3.*Sqrt(1 + mf2)*zf)))/2.


  else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!withsea part
  dr0dtV=(k**4*v3*((2*(1 - etaphi/5.)*(0.5 + nbSigma))/(3.*Sqrt(1 + ms2)*Sqrt(zb)) +   &
      (2*(1 - etaphi/5.)*(0.5 + nbPion)*(-1 + Nf**2))/                          &
       (3.*Sqrt(1 + mp2)*Sqrt(zb)) -                                            &
      (4*(1 - etapsi/4.)*Nc*Nf*(1 - nfa - nff))/(3.*Sqrt(1 + mf2)*zf)))/2.

  dr1dtV=(k**4*v3*(((1 - etaphi/5.)*k*ms2d1rho*nbd1xSigma)/(3.*(1 + ms2)*zb) +    &
      ((1 - etaphi/5.)*k*mp2d1rho*nbd1xPion*(-1 + Nf**2))/                      &
       (3.*(1 + mp2)*zb) - ((1 - etaphi/5.)*ms2d1rho*(0.5 + nbSigma))/          &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*mp2d1rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                           &
      (2*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*(1 - nfa - nff))/                       &
       (3.*(1 + mf2)**1.5*zf) -                                                 &
      (4*(1 - etapsi/4.)*Nc*Nf*                                                 &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/(3.*Sqrt(1 + mf2)*zf)    &
))/2.

  dr2dtV=(k**4*v3*(((1 - etaphi/5.)*k**2*ms2d1rho**2*nbd2xSigma)/               &
       (6.*(1 + ms2)**1.5*zb**1.5) +                                            &
      ((1 - etaphi/5.)*k**2*mp2d1rho**2*nbd2xPion*(-1 + Nf**2))/                &
       (6.*(1 + mp2)**1.5*zb**1.5) -                                            &
      ((1 - etaphi/5.)*k*ms2d1rho**2*nbd1xSigma)/(2.*(1 + ms2)**2*zb) +         &
      ((1 - etaphi/5.)*k*ms2d2rho*nbd1xSigma)/(3.*(1 + ms2)*zb) -               &
      ((1 - etaphi/5.)*k*mp2d1rho**2*nbd1xPion*(-1 + Nf**2))/                   &
       (2.*(1 + mp2)**2*zb) + ((1 - etaphi/5.)*k*mp2d2rho*nbd1xPion*            &
         (-1 + Nf**2))/(3.*(1 + mp2)*zb) +                                      &
      ((1 - etaphi/5.)*ms2d1rho**2*(0.5 + nbSigma))/                            &
       (2.*(1 + ms2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*ms2d2rho*(0.5 + nbSigma))/                               &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) +                                           &
      ((1 - etaphi/5.)*mp2d1rho**2*(0.5 + nbPion)*(-1 + Nf**2))/                &
       (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*mp2d2rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) -                                           &
      ((1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*(1 - nfa - nff))/                      &
       ((1 + mf2)**2.5*zf) + (4*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                 &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                         &
       (3.*(1 + mf2)**1.5*zf) -                                                 &
      (4*(1 - etapsi/4.)*Nc*Nf*                                                 &
         (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                     &
           (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                     &
           (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                      &
           (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                     &
       (3.*Sqrt(1 + mf2)*zf)))/2.

  dr3dtV=(k**4*v3*(((1 - etaphi/5.)*k**3*ms2d1rho**3*nbd3xSigma)/               &
       (12.*(1 + ms2)**2*zb**2) +                                               &
      ((1 - etaphi/5.)*k**3*mp2d1rho**3*nbd3xPion*(-1 + Nf**2))/                &
       (12.*(1 + mp2)**2*zb**2) -                                               &
      ((1 - etaphi/5.)*k**2*ms2d1rho**3*nbd2xSigma)/                            &
       (2.*(1 + ms2)**2.5*zb**1.5) +                                            &
      ((1 - etaphi/5.)*k**2*ms2d1rho*ms2d2rho*nbd2xSigma)/                      &
       (2.*(1 + ms2)**1.5*zb**1.5) -                                            &
      ((1 - etaphi/5.)*k**2*mp2d1rho**3*nbd2xPion*(-1 + Nf**2))/                &
       (2.*(1 + mp2)**2.5*zb**1.5) +                                            &
      ((1 - etaphi/5.)*k**2*mp2d1rho*mp2d2rho*nbd2xPion*(-1 + Nf**2))/          &
       (2.*(1 + mp2)**1.5*zb**1.5) +                                            &
      (5*(1 - etaphi/5.)*k*ms2d1rho**3*nbd1xSigma)/(4.*(1 + ms2)**3*zb) -       &
      (3*(1 - etaphi/5.)*k*ms2d1rho*ms2d2rho*nbd1xSigma)/                       &
       (2.*(1 + ms2)**2*zb) + ((1 - etaphi/5.)*k*ms2d3rho*nbd1xSigma)/          &
       (3.*(1 + ms2)*zb) + (5*(1 - etaphi/5.)*k*mp2d1rho**3*nbd1xPion*          &
         (-1 + Nf**2))/(4.*(1 + mp2)**3*zb) -                                   &
      (3*(1 - etaphi/5.)*k*mp2d1rho*mp2d2rho*nbd1xPion*(-1 + Nf**2))/           &
       (2.*(1 + mp2)**2*zb) + ((1 - etaphi/5.)*k*mp2d3rho*nbd1xPion*            &
         (-1 + Nf**2))/(3.*(1 + mp2)*zb) -                                      &
      (5*(1 - etaphi/5.)*ms2d1rho**3*(0.5 + nbSigma))/                          &
       (4.*(1 + ms2)**3.5*Sqrt(zb)) +                                           &
      (3*(1 - etaphi/5.)*ms2d1rho*ms2d2rho*(0.5 + nbSigma))/                    &
       (2.*(1 + ms2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*ms2d3rho*(0.5 + nbSigma))/                               &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                           &
      (5*(1 - etaphi/5.)*mp2d1rho**3*(0.5 + nbPion)*(-1 + Nf**2))/              &
       (4.*(1 + mp2)**3.5*Sqrt(zb)) +                                           &
      (3*(1 - etaphi/5.)*mp2d1rho*mp2d2rho*(0.5 + nbPion)*(-1 + Nf**2))/        &
       (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                           &
      (5*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*(1 - nfa - nff))/                    &
       (2.*(1 + mf2)**3.5*zf) -                                                 &
      (3*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*                                     &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/((1 + mf2)**2.5*zf)      &
+ (2*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                                            &
         (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                     &
           (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                     &
           (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                      &
           (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                     &
       ((1 + mf2)**1.5*zf) - (4*(1 - etapsi/4.)*Nc*Nf*                          &
         (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -                &
           (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +                &
           (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +                &
           (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -                &
           (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                    &
           (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                   &
       (3.*Sqrt(1 + mf2)*zf)))/2.

  dr4dtV=(k**4*v3*(((1 - etaphi/5.)*k**4*ms2d1rho**4*nbd4xSigma)/               &
       (24.*(1 + ms2)**2.5*zb**2.5) +                                           &
      ((1 - etaphi/5.)*k**4*mp2d1rho**4*nbd4xPion*(-1 + Nf**2))/                &
       (24.*(1 + mp2)**2.5*zb**2.5) -                                           &
      (5*(1 - etaphi/5.)*k**3*ms2d1rho**4*nbd3xSigma)/                          &
       (12.*(1 + ms2)**3*zb**2) +                                               &
      ((1 - etaphi/5.)*k**3*ms2d1rho**2*ms2d2rho*nbd3xSigma)/                   &
       (2.*(1 + ms2)**2*zb**2) -                                                &
      (5*(1 - etaphi/5.)*k**3*mp2d1rho**4*nbd3xPion*(-1 + Nf**2))/              &
       (12.*(1 + mp2)**3*zb**2) +                                               &
      ((1 - etaphi/5.)*k**3*mp2d1rho**2*mp2d2rho*nbd3xPion*(-1 + Nf**2))/       &
       (2.*(1 + mp2)**2*zb**2) +                                                &
      (15*(1 - etaphi/5.)*k**2*ms2d1rho**4*nbd2xSigma)/                         &
       (8.*(1 + ms2)**3.5*zb**1.5) -                                            &
      (3*(1 - etaphi/5.)*k**2*ms2d1rho**2*ms2d2rho*nbd2xSigma)/                 &
       ((1 + ms2)**2.5*zb**1.5) +                                               &
      ((1 - etaphi/5.)*k**2*ms2d2rho**2*nbd2xSigma)/                            &
       (2.*(1 + ms2)**1.5*zb**1.5) +                                            &
      (2*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d3rho*nbd2xSigma)/                    &
       (3.*(1 + ms2)**1.5*zb**1.5) +                                            &
      (15*(1 - etaphi/5.)*k**2*mp2d1rho**4*nbd2xPion*(-1 + Nf**2))/             &
       (8.*(1 + mp2)**3.5*zb**1.5) -                                            &
      (3*(1 - etaphi/5.)*k**2*mp2d1rho**2*mp2d2rho*nbd2xPion*(-1 + Nf**2))/     &
       ((1 + mp2)**2.5*zb**1.5) +                                               &
      ((1 - etaphi/5.)*k**2*mp2d2rho**2*nbd2xPion*(-1 + Nf**2))/                &
       (2.*(1 + mp2)**1.5*zb**1.5) +                                            &
      (2*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d3rho*nbd2xPion*(-1 + Nf**2))/        &
       (3.*(1 + mp2)**1.5*zb**1.5) -                                            &
      (35*(1 - etaphi/5.)*k*ms2d1rho**4*nbd1xSigma)/(8.*(1 + ms2)**4*zb) +      &
      (15*(1 - etaphi/5.)*k*ms2d1rho**2*ms2d2rho*nbd1xSigma)/                   &
       (2.*(1 + ms2)**3*zb) - (3*(1 - etaphi/5.)*k*ms2d2rho**2*nbd1xSigma)/     &
       (2.*(1 + ms2)**2*zb) - (2*(1 - etaphi/5.)*k*ms2d1rho*ms2d3rho*           &
         nbd1xSigma)/((1 + ms2)**2*zb) +                                        &
      ((1 - etaphi/5.)*k*ms2d4rho*nbd1xSigma)/(3.*(1 + ms2)*zb) -               &
      (35*(1 - etaphi/5.)*k*mp2d1rho**4*nbd1xPion*(-1 + Nf**2))/                &
       (8.*(1 + mp2)**4*zb) + (15*(1 - etaphi/5.)*k*mp2d1rho**2*mp2d2rho*       &
         nbd1xPion*(-1 + Nf**2))/(2.*(1 + mp2)**3*zb) -                         &
      (3*(1 - etaphi/5.)*k*mp2d2rho**2*nbd1xPion*(-1 + Nf**2))/                 &
       (2.*(1 + mp2)**2*zb) - (2*(1 - etaphi/5.)*k*mp2d1rho*mp2d3rho*           &
         nbd1xPion*(-1 + Nf**2))/((1 + mp2)**2*zb) +                            &
      ((1 - etaphi/5.)*k*mp2d4rho*nbd1xPion*(-1 + Nf**2))/                      &
       (3.*(1 + mp2)*zb) + (35*(1 - etaphi/5.)*ms2d1rho**4*                     &
         (0.5 + nbSigma))/(8.*(1 + ms2)**4.5*Sqrt(zb)) -                        &
      (15*(1 - etaphi/5.)*ms2d1rho**2*ms2d2rho*(0.5 + nbSigma))/                &
       (2.*(1 + ms2)**3.5*Sqrt(zb)) +                                           &
      (3*(1 - etaphi/5.)*ms2d2rho**2*(0.5 + nbSigma))/                          &
       (2.*(1 + ms2)**2.5*Sqrt(zb)) +                                           &
      (2*(1 - etaphi/5.)*ms2d1rho*ms2d3rho*(0.5 + nbSigma))/                    &
       ((1 + ms2)**2.5*Sqrt(zb)) -                                              &
      ((1 - etaphi/5.)*ms2d4rho*(0.5 + nbSigma))/                               &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) +                                           &
      (35*(1 - etaphi/5.)*mp2d1rho**4*(0.5 + nbPion)*(-1 + Nf**2))/             &
       (8.*(1 + mp2)**4.5*Sqrt(zb)) -                                           &
      (15*(1 - etaphi/5.)*mp2d1rho**2*mp2d2rho*(0.5 + nbPion)*                  &
         (-1 + Nf**2))/(2.*(1 + mp2)**3.5*Sqrt(zb)) +                           &
      (3*(1 - etaphi/5.)*mp2d2rho**2*(0.5 + nbPion)*(-1 + Nf**2))/              &
       (2.*(1 + mp2)**2.5*Sqrt(zb)) +                                           &
      (2*(1 - etaphi/5.)*mp2d1rho*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/        &
       ((1 + mp2)**2.5*Sqrt(zb)) -                                              &
      ((1 - etaphi/5.)*mp2d4rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) -                                           &
      (35*(1 - etapsi/4.)*mf2d1rho**4*Nc*Nf*(1 - nfa - nff))/                   &
       (4.*(1 + mf2)**4.5*zf) +                                                 &
      (10*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*                                    &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/((1 + mf2)**3.5*zf)      &
- (6*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*                                         &
         (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                     &
           (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                     &
           (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                      &
           (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                     &
       ((1 + mf2)**2.5*zf) + (8*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                 &
         (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -                &
           (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +                &
           (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +                &
           (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -                &
           (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                    &
           (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                   &
       (3.*(1 + mf2)**1.5*zf) -                                                 &
      (4*(1 - etapsi/4.)*Nc*Nf*                                                 &
         (-(k**4*mf2d1rho**4*nfd4xa)/(16.*(1 + mf2)**2*zf**4) -                 &
           (k**4*mf2d1rho**4*nfd4xf)/(16.*(1 + mf2)**2*zf**4) +                 &
           (3*k**3*mf2d1rho**4*nfd3xa)/(8.*(1 + mf2)**2.5*zf**3) +              &
           (3*k**3*mf2d1rho**4*nfd3xf)/(8.*(1 + mf2)**2.5*zf**3) -              &
           (15*k**2*mf2d1rho**4*nfd2xa)/(16.*(1 + mf2)**3*zf**2) -              &
           (15*k**2*mf2d1rho**4*nfd2xf)/(16.*(1 + mf2)**3*zf**2) +              &
           (15*k*mf2d1rho**4*nfd1xa)/(16.*(1 + mf2)**3.5*zf) +                  &
           (15*k*mf2d1rho**4*nfd1xf)/(16.*(1 + mf2)**3.5*zf)))/                 &
       (3.*Sqrt(1 + mf2)*zf)))/2.

  dr5dtV=(k**4*v3*(((1 - etaphi/5.)*k**5*ms2d1rho**5*nbd5xSigma)/               &
       (48.*(1 + ms2)**3*zb**3) +                                               &
      ((1 - etaphi/5.)*k**5*mp2d1rho**5*nbd5xPion*(-1 + Nf**2))/                &
       (48.*(1 + mp2)**3*zb**3) -                                               &
      (5*(1 - etaphi/5.)*k**4*ms2d1rho**5*nbd4xSigma)/                          &
       (16.*(1 + ms2)**3.5*zb**2.5) +                                           &
      (5*(1 - etaphi/5.)*k**4*ms2d1rho**3*ms2d2rho*nbd4xSigma)/                 &
       (12.*(1 + ms2)**2.5*zb**2.5) -                                           &
      (5*(1 - etaphi/5.)*k**4*mp2d1rho**5*nbd4xPion*(-1 + Nf**2))/              &
       (16.*(1 + mp2)**3.5*zb**2.5) +                                           &
      (5*(1 - etaphi/5.)*k**4*mp2d1rho**3*mp2d2rho*nbd4xPion*(-1 + Nf**2))/     &
       (12.*(1 + mp2)**2.5*zb**2.5) +                                           &
      (35*(1 - etaphi/5.)*k**3*ms2d1rho**5*nbd3xSigma)/                         &
       (16.*(1 + ms2)**4*zb**2) -                                               &
      (25*(1 - etaphi/5.)*k**3*ms2d1rho**3*ms2d2rho*nbd3xSigma)/                &
       (6.*(1 + ms2)**3*zb**2) +                                                &
      (5*(1 - etaphi/5.)*k**3*ms2d1rho*ms2d2rho**2*nbd3xSigma)/                 &
       (4.*(1 + ms2)**2*zb**2) +                                                &
      (5*(1 - etaphi/5.)*k**3*ms2d1rho**2*ms2d3rho*nbd3xSigma)/                 &
       (6.*(1 + ms2)**2*zb**2) +                                                &
      (35*(1 - etaphi/5.)*k**3*mp2d1rho**5*nbd3xPion*(-1 + Nf**2))/             &
       (16.*(1 + mp2)**4*zb**2) -                                               &
      (25*(1 - etaphi/5.)*k**3*mp2d1rho**3*mp2d2rho*nbd3xPion*                  &
         (-1 + Nf**2))/(6.*(1 + mp2)**3*zb**2) +                                &
      (5*(1 - etaphi/5.)*k**3*mp2d1rho*mp2d2rho**2*nbd3xPion*(-1 + Nf**2))/     &
       (4.*(1 + mp2)**2*zb**2) +                                                &
      (5*(1 - etaphi/5.)*k**3*mp2d1rho**2*mp2d3rho*nbd3xPion*(-1 + Nf**2))/     &
       (6.*(1 + mp2)**2*zb**2) -                                                &
      (35*(1 - etaphi/5.)*k**2*ms2d1rho**5*nbd2xSigma)/                         &
       (4.*(1 + ms2)**4.5*zb**1.5) +                                            &
      (75*(1 - etaphi/5.)*k**2*ms2d1rho**3*ms2d2rho*nbd2xSigma)/                &
       (4.*(1 + ms2)**3.5*zb**1.5) -                                            &
      (15*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d2rho**2*nbd2xSigma)/                &
       (2.*(1 + ms2)**2.5*zb**1.5) -                                            &
      (5*(1 - etaphi/5.)*k**2*ms2d1rho**2*ms2d3rho*nbd2xSigma)/                 &
       ((1 + ms2)**2.5*zb**1.5) +                                               &
      (5*(1 - etaphi/5.)*k**2*ms2d2rho*ms2d3rho*nbd2xSigma)/                    &
       (3.*(1 + ms2)**1.5*zb**1.5) +                                            &
      (5*(1 - etaphi/5.)*k**2*ms2d1rho*ms2d4rho*nbd2xSigma)/                    &
       (6.*(1 + ms2)**1.5*zb**1.5) -                                            &
      (35*(1 - etaphi/5.)*k**2*mp2d1rho**5*nbd2xPion*(-1 + Nf**2))/             &
       (4.*(1 + mp2)**4.5*zb**1.5) +                                            &
      (75*(1 - etaphi/5.)*k**2*mp2d1rho**3*mp2d2rho*nbd2xPion*                  &
         (-1 + Nf**2))/(4.*(1 + mp2)**3.5*zb**1.5) -                            &
      (15*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d2rho**2*nbd2xPion*                  &
         (-1 + Nf**2))/(2.*(1 + mp2)**2.5*zb**1.5) -                            &
      (5*(1 - etaphi/5.)*k**2*mp2d1rho**2*mp2d3rho*nbd2xPion*(-1 + Nf**2))/     &
       ((1 + mp2)**2.5*zb**1.5) +                                               &
      (5*(1 - etaphi/5.)*k**2*mp2d2rho*mp2d3rho*nbd2xPion*(-1 + Nf**2))/        &
       (3.*(1 + mp2)**1.5*zb**1.5) +                                            &
      (5*(1 - etaphi/5.)*k**2*mp2d1rho*mp2d4rho*nbd2xPion*(-1 + Nf**2))/        &
       (6.*(1 + mp2)**1.5*zb**1.5) +                                            &
      (315*(1 - etaphi/5.)*k*ms2d1rho**5*nbd1xSigma)/                           &
       (16.*(1 + ms2)**5*zb) -                                                  &
      (175*(1 - etaphi/5.)*k*ms2d1rho**3*ms2d2rho*nbd1xSigma)/                  &
       (4.*(1 + ms2)**4*zb) + (75*(1 - etaphi/5.)*k*ms2d1rho*ms2d2rho**2*       &
         nbd1xSigma)/(4.*(1 + ms2)**3*zb) +                                     &
      (25*(1 - etaphi/5.)*k*ms2d1rho**2*ms2d3rho*nbd1xSigma)/                   &
       (2.*(1 + ms2)**3*zb) - (5*(1 - etaphi/5.)*k*ms2d2rho*ms2d3rho*           &
         nbd1xSigma)/((1 + ms2)**2*zb) -                                        &
      (5*(1 - etaphi/5.)*k*ms2d1rho*ms2d4rho*nbd1xSigma)/                       &
       (2.*(1 + ms2)**2*zb) + ((1 - etaphi/5.)*k*ms2d5rho*nbd1xSigma)/          &
       (3.*(1 + ms2)*zb) + (315*(1 - etaphi/5.)*k*mp2d1rho**5*nbd1xPion*        &
         (-1 + Nf**2))/(16.*(1 + mp2)**5*zb) -                                  &
      (175*(1 - etaphi/5.)*k*mp2d1rho**3*mp2d2rho*nbd1xPion*(-1 + Nf**2))/      &
       (4.*(1 + mp2)**4*zb) + (75*(1 - etaphi/5.)*k*mp2d1rho*mp2d2rho**2*       &
         nbd1xPion*(-1 + Nf**2))/(4.*(1 + mp2)**3*zb) +                         &
      (25*(1 - etaphi/5.)*k*mp2d1rho**2*mp2d3rho*nbd1xPion*(-1 + Nf**2))/       &
       (2.*(1 + mp2)**3*zb) - (5*(1 - etaphi/5.)*k*mp2d2rho*mp2d3rho*           &
         nbd1xPion*(-1 + Nf**2))/((1 + mp2)**2*zb) -                            &
      (5*(1 - etaphi/5.)*k*mp2d1rho*mp2d4rho*nbd1xPion*(-1 + Nf**2))/           &
       (2.*(1 + mp2)**2*zb) + ((1 - etaphi/5.)*k*mp2d5rho*nbd1xPion*            &
         (-1 + Nf**2))/(3.*(1 + mp2)*zb) -                                      &
      (315*(1 - etaphi/5.)*ms2d1rho**5*(0.5 + nbSigma))/                        &
       (16.*(1 + ms2)**5.5*Sqrt(zb)) +                                          &
      (175*(1 - etaphi/5.)*ms2d1rho**3*ms2d2rho*(0.5 + nbSigma))/               &
       (4.*(1 + ms2)**4.5*Sqrt(zb)) -                                           &
      (75*(1 - etaphi/5.)*ms2d1rho*ms2d2rho**2*(0.5 + nbSigma))/                &
       (4.*(1 + ms2)**3.5*Sqrt(zb)) -                                           &
      (25*(1 - etaphi/5.)*ms2d1rho**2*ms2d3rho*(0.5 + nbSigma))/                &
       (2.*(1 + ms2)**3.5*Sqrt(zb)) +                                           &
      (5*(1 - etaphi/5.)*ms2d2rho*ms2d3rho*(0.5 + nbSigma))/                    &
       ((1 + ms2)**2.5*Sqrt(zb)) +                                              &
      (5*(1 - etaphi/5.)*ms2d1rho*ms2d4rho*(0.5 + nbSigma))/                    &
       (2.*(1 + ms2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*ms2d5rho*(0.5 + nbSigma))/                               &
       (3.*(1 + ms2)**1.5*Sqrt(zb)) -                                           &
      (315*(1 - etaphi/5.)*mp2d1rho**5*(0.5 + nbPion)*(-1 + Nf**2))/            &
       (16.*(1 + mp2)**5.5*Sqrt(zb)) +                                          &
      (175*(1 - etaphi/5.)*mp2d1rho**3*mp2d2rho*(0.5 + nbPion)*                 &
         (-1 + Nf**2))/(4.*(1 + mp2)**4.5*Sqrt(zb)) -                           &
      (75*(1 - etaphi/5.)*mp2d1rho*mp2d2rho**2*(0.5 + nbPion)*                  &
         (-1 + Nf**2))/(4.*(1 + mp2)**3.5*Sqrt(zb)) -                           &
      (25*(1 - etaphi/5.)*mp2d1rho**2*mp2d3rho*(0.5 + nbPion)*                  &
         (-1 + Nf**2))/(2.*(1 + mp2)**3.5*Sqrt(zb)) +                           &
      (5*(1 - etaphi/5.)*mp2d2rho*mp2d3rho*(0.5 + nbPion)*(-1 + Nf**2))/        &
       ((1 + mp2)**2.5*Sqrt(zb)) +                                              &
      (5*(1 - etaphi/5.)*mp2d1rho*mp2d4rho*(0.5 + nbPion)*(-1 + Nf**2))/        &
       (2.*(1 + mp2)**2.5*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.)*mp2d5rho*(0.5 + nbPion)*(-1 + Nf**2))/                   &
       (3.*(1 + mp2)**1.5*Sqrt(zb)) +                                           &
      (315*(1 - etapsi/4.)*mf2d1rho**5*Nc*Nf*(1 - nfa - nff))/                  &
       (8.*(1 + mf2)**5.5*zf) -                                                 &
      (175*(1 - etapsi/4.)*mf2d1rho**4*Nc*Nf*                                   &
         (-(k*mf2d1rho*nfd1xa)/(2.*Sqrt(1 + mf2)*zf) -                          &
           (k*mf2d1rho*nfd1xf)/(2.*Sqrt(1 + mf2)*zf)))/                         &
       (4.*(1 + mf2)**4.5*zf) +                                                 &
      (25*(1 - etapsi/4.)*mf2d1rho**3*Nc*Nf*                                    &
         (-(k**2*mf2d1rho**2*nfd2xa)/(4.*(1 + mf2)*zf**2) -                     &
           (k**2*mf2d1rho**2*nfd2xf)/(4.*(1 + mf2)*zf**2) +                     &
           (k*mf2d1rho**2*nfd1xa)/(4.*(1 + mf2)**1.5*zf) +                      &
           (k*mf2d1rho**2*nfd1xf)/(4.*(1 + mf2)**1.5*zf)))/                     &
       ((1 + mf2)**3.5*zf) - (10*(1 - etapsi/4.)*mf2d1rho**2*Nc*Nf*             &
         (-(k**3*mf2d1rho**3*nfd3xa)/(8.*(1 + mf2)**1.5*zf**3) -                &
           (k**3*mf2d1rho**3*nfd3xf)/(8.*(1 + mf2)**1.5*zf**3) +                &
           (3*k**2*mf2d1rho**3*nfd2xa)/(8.*(1 + mf2)**2*zf**2) +                &
           (3*k**2*mf2d1rho**3*nfd2xf)/(8.*(1 + mf2)**2*zf**2) -                &
           (3*k*mf2d1rho**3*nfd1xa)/(8.*(1 + mf2)**2.5*zf) -                    &
           (3*k*mf2d1rho**3*nfd1xf)/(8.*(1 + mf2)**2.5*zf)))/                   &
       ((1 + mf2)**2.5*zf) + (10*(1 - etapsi/4.)*mf2d1rho*Nc*Nf*                &
         (-(k**4*mf2d1rho**4*nfd4xa)/(16.*(1 + mf2)**2*zf**4) -                 &
           (k**4*mf2d1rho**4*nfd4xf)/(16.*(1 + mf2)**2*zf**4) +                 &
           (3*k**3*mf2d1rho**4*nfd3xa)/(8.*(1 + mf2)**2.5*zf**3) +              &
           (3*k**3*mf2d1rho**4*nfd3xf)/(8.*(1 + mf2)**2.5*zf**3) -              &
           (15*k**2*mf2d1rho**4*nfd2xa)/(16.*(1 + mf2)**3*zf**2) -              &
           (15*k**2*mf2d1rho**4*nfd2xf)/(16.*(1 + mf2)**3*zf**2) +              &
           (15*k*mf2d1rho**4*nfd1xa)/(16.*(1 + mf2)**3.5*zf) +                  &
           (15*k*mf2d1rho**4*nfd1xf)/(16.*(1 + mf2)**3.5*zf)))/                 &
       (3.*(1 + mf2)**1.5*zf) -                                                 &
      (4*(1 - etapsi/4.)*Nc*Nf*                                                 &
         (-(k**5*mf2d1rho**5*nfd5xa)/(32.*(1 + mf2)**2.5*zf**5) -               &
           (k**5*mf2d1rho**5*nfd5xf)/(32.*(1 + mf2)**2.5*zf**5) +               &
           (5*k**4*mf2d1rho**5*nfd4xa)/(16.*(1 + mf2)**3*zf**4) +               &
           (5*k**4*mf2d1rho**5*nfd4xf)/(16.*(1 + mf2)**3*zf**4) -               &
           (45*k**3*mf2d1rho**5*nfd3xa)/(32.*(1 + mf2)**3.5*zf**3) -            &
           (45*k**3*mf2d1rho**5*nfd3xf)/(32.*(1 + mf2)**3.5*zf**3) +            &
           (105*k**2*mf2d1rho**5*nfd2xa)/(32.*(1 + mf2)**4*zf**2) +             &
           (105*k**2*mf2d1rho**5*nfd2xf)/(32.*(1 + mf2)**4*zf**2) -             &
           (105*k*mf2d1rho**5*nfd1xa)/(32.*(1 + mf2)**4.5*zf) -                 &
           (105*k*mf2d1rho**5*nfd1xf)/(32.*(1 + mf2)**4.5*zf)))/                &
       (3.*Sqrt(1 + mf2)*zf)))/2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end if



  dlam0dt=(dr0dtV)
  dlam1dt=(etaphi*lam1+dr1dtV)
  dlam2dt=(2.*etaphi*lam2+dr2dtV)
  dlam3dt=(3.*etaphi*lam3+dr3dtV)
  dlam4dt=(4.*etaphi*lam4+dr4dtV)
  dlam5dt=(5.*etaphi*lam5+dr5dtV)
!  dlam6dt=(6.*etaphi*lam6+dr6dtV)
!  dlam7dt=(7.*etaphi*lam7+dr7dtV)
  dlam6dt=0.
  dlam7dt=0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  L11Pion=(2*(b2f1aPion*(1 - etaphi/5.) + b1f2Pion*(1 - etapsi/4.)))/3.
  L11Sigma=(2*(b2f1aSigma*(1 - etaphi/5.) + b1f2Sigma*(1 - etapsi/4.)))/3.

  dth=(etaphi/2. + etapsi)*h + (((h**3*L11Sigma)/Nf -                           &
       (h**3*L11Pion*(-1 + Nf**2))/Nf)*v3)/2.

  dhdt=dth

  dydx(1)=dlam1dt
  dydx(2)=dlam2dt
  dydx(3)=dlam3dt
  dydx(4)=dlam4dt
  dydx(5)=dlam5dt
  dydx(6)=dlam6dt
  dydx(7)=dlam7dt
  dydx(Nv+1)=dlam0dt
  dydx((Nv+1)+1)=dhdt
!  dydx((Nv+1)+1)=0.
  dydx((Nv+1)+(Nh+1)+1)=dZphidt
  dydx((Nv+1)+(Nh+1)+2)=dZpsidt
  dydx((Nv+1)+(Nh+1)+Nz+1)=dcdt
  dydx((Nv+1)+(Nh+1)+Nz+2)=dkappadt

  goto 100

  open(unit=101,file='./buffer/k1.dat')
  write(101, "(e20.9)")k*hc

  open(unit=102,file='./buffer/etaphi.dat')
  write(102, "(e20.9)")etaphi

  open(unit=103,file='./buffer/etapsi.dat')
  write(103, "(e20.9)")etapsi

  open(unit=104,file='./buffer/y1.dat')
  write(104, "(e20.9)")y(1)

  open(unit=105,file='./buffer/y2.dat')
  write(105, "(e20.9)")y(2)

  open(unit=106,file='./buffer/y3.dat')
  write(106, "(e20.9)")y(3)

  open(unit=107,file='./buffer/y4.dat')
  write(107, "(e20.9)")y(4)

  open(unit=108,file='./buffer/y5.dat')
  write(108, "(e20.9)")y(5)

  open(unit=109,file='./buffer/y6.dat')
  write(109, "(e20.9)")y(6)

  open(unit=110,file='./buffer/y7.dat')
  write(110, "(e20.9)")y(7)

  open(unit=111,file='./buffer/y8.dat')
  write(111, "(e20.9)")y(8)

  open(unit=112,file='./buffer/y9.dat')
  write(112, "(e20.9)")y(9)

  open(unit=113,file='./buffer/y10.dat')
  write(113, "(e20.9)")y(10)

  open(unit=114,file='./buffer/y11.dat')
  write(114, "(e20.9)")y(11)

  open(unit=115,file='./buffer/y12.dat')
  write(115, "(e20.9)")y(12)

  open(unit=116,file='./buffer/y13.dat')
  write(116, "(e20.9)")y(13)

  open(unit=117,file='./buffer/t1.dat')
  write(117, "(e20.9)")x


100 continue

end



