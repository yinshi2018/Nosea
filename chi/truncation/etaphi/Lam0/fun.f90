function Fnb(x,T)          !boson distribution function
  implicit none

  real(8) x, T
  real(8) Fnb
  real(8) over

  over=x/T
  if(over > 60.)then
    Fnb=0.
  else
    Fnb=1./(exp(over)-1.)
  end if

  return
end function Fnb

function Fnf0(x,T,l,lb)         !fermion distribution function

  implicit none

  real(8) x,T,l,lb
  real(8) Fnf0
  real(8) over

  over=x/T
  if(over > 60.)then
    Fnf0=0.
  else
    Fnf0=1./(1.+3.*lb*exp(x/T)+3.*l*exp(2.*x/T)+exp(3.*x/T))
  end if

  return
end function Fnf0


function Fnf1(x,T,l,lb)         !fermion distribution function
  implicit none

  real(8) x,T,l,lb
  real(8) Fnf1
  real(8) over

  over=x/T
  if(over > 60.)then
    Fnf1=0.
  else
    Fnf1=(2*exp(x/T))/(1.+3.*lb*exp(x/T)+3.*l*exp(2.*x/T)+exp(3.*x/T))
  end if

  return
end function Fnf1

function Fnf2(x,T,l,lb)         !fermion distribution function
  implicit none

  real(8) x,T,l,lb
  real(8) Fnf2
  real(8) over

  over=x/T
  if(over > 60.)then
    Fnf2=0.
  else
    Fnf2=(exp(2.*x/T))/(1.+3.*lb*exp(x/T)+3.*l*exp(2.*x/T)+exp(3.*x/T))
  end if

  return
end function Fnf2


