!the function form of Morse and its paritial deviation of D,rm and alpha
real(kind=8) function dmorse_d(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dmorse_d=exp(a*(1.0-L/r))-2*exp(0.5*a*(1.0-L/r))
	

return
end function


real(kind=8) function dmorse_r(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dmorse_r=(d*a*L/(r*r))*(exp(a*(1.0-L/r))-exp(0.5*a*(1.0-L/r)))
	

return
end function

real(kind=8) function dmorse_a(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dmorse_a=d*(1.0-L/r)*(exp(a*(1.0-L/r))-exp(0.5*a*(1.0-L/r)))
	

return
end function

real(kind=8) function morse(i,j,L)
use globe
implicit none
real(kind=8)::		L
integer		::		i,j
	
	morse=d(i,j)*(exp(a(i,j)*(1.0-L/r(i,j)))-2.0*exp(0.5*a(i,j)*(1.0-L/r(i,j))))
	
	
return	
end function

!the function form of Morse-4-2 and its paritial deviation of D,rm and alpha
real(kind=8) function dmorse42_d(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dmorse42_d=exp(a*(1.0-L/r))-((L/r)**4-2.0*(L/r)**2+3.0)*exp(0.5*a*(1.0-L/r))

return
end function


real(kind=8) function dmorse42_r(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dmorse42_r=(d*a*L/(r*r))*(exp(a*(1.0-L/r))-((L/r)**4-2.0*(L/r)**2+3)*0.5*exp(0.5*a*(1.0-L/r)))-&
			& d/L*(-4.0*(L/r)**5+4.0*(L/r)**3)*exp(0.5*a*(1.0-L/r))

return
end function

real(kind=8) function dmorse42_a(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dmorse42_a=d*(1.0-L/r)*(exp(a*(1.0-L/r))-((L/r)**4-2*(L/r)**2+3)*0.5*exp(0.5*a*(1.0-L/r)))

return
end function

real(kind=8) function morse42(i,j,L)
use globe
implicit none
real(kind=8)::		L
integer		::		i,j
	
	morse42=d(i,j)*(exp(a(i,j)*(1.0-L/r(i,j)))-((L/r(i,j))**4-2*(L/r(i,j))**2+3)*exp(0.5*a(i,j)*(1.0-L/r(i,j))))
	
return	
end function


!the function form of Lennard Jones and its paritial deviation of D,rm and alpha
real(kind=8) function dlj_d(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dlj_d=(r/L)**12-2.0*(r/L)**6
	

return
end function


real(kind=8) function dlj_r(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dlj_r=d/r*12.0*((r/L)**12-(r/L)**6)
	

return
end function

real(kind=8) function dlj_a(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dlj_a=0
	

return
end function

real(kind=8) function lj(i,j,L)
use globe
implicit none
real(kind=8)::		L
integer		::		i,j
	
	lj=d(i,j)*((r(i,j)/L)**12-2.0*(r(i,j)/L)**6)
	
	
return	
end function

!the function form of Exp-6 and its paritial deviation of D,rm and alpha
real(kind=8) function dexp6_d(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dexp6_d=6.0/(a-6.0)*exp(a*(1.0-L/r))-a/(a-6.0)*(r/L)**6

return
end function


real(kind=8) function dexp6_r(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dexp6_r=(d*a/r)*6.0/(a-6.0)*(L/r*exp(a*(1.0-L/r))-(r/L)**6)
			

return
end function

real(kind=8) function dexp6_a(d,r,a,L)
reaL(kind=8) :: d,r,a
reaL(kind=8) :: L   !L represent distance

	dexp6_a=d*(6.0/(a-6.0)**2)*(-exp(a*(1.0-L/r))+(r/L)**6)+&
	        & 6.0/(a-6.0)*(1.0-L/r)*exp(a*(1.0-L/r))

return
end function

real(kind=8) function exp6(i,j,L)
use globe
implicit none
real(kind=8)::		L
integer		::		i,j
	
	exp6=d(i,j)*(6.0/(a(i,j)-6.0)*exp(a(i,j)*(1.0-L/r(i,j)))-a(i,j)/(a(i,j)-6.0)*(r(i,j)/L)**6)
	
return	
end function


!the combination rules function for D, Rm and alpha
real(kind=8) function cbr_d(d1,d2)
reaL(kind=8) :: d1,d2

	cbr_d=2*(d1*d2)/(d1+d2)

return
end function

real(kind=8) function cbr_r(r1,r2)
reaL(kind=8) :: r1,r2

	cbr_r=r1*r2*(r1+r2)/(r1**2+r2**2)
	!cbr_r=0.5*(r1+r2)

return
end function

real(kind=8) function cbr_a(a1,a2,r1,r2,r12)
reaL(kind=8) :: a1,a2,r1,r2,r12
reaL(kind=8) :: sig1,sig2,sig12

	sig1=r1*(1-2.0*log(2.0)/a1)
	sig2=r2*(1-2.0*log(2.0)/a2)
	sig12=sig1*sig2*(sig1+sig2)/(sig1**2+sig2**2)
	cbr_a=2.0*log(2.0)/(1.0-sig12/r12)
	!cbr_a=0.5*(a1+a2)
	
return
end function


