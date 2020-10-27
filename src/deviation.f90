subroutine deviation(Func,x0,i,j,L,T3)
use globe
implicit none
real(kind=8),external	::	Func
real(kind=8)			::	x0,L
real(kind=8)			::	yi(6),inix
real(kind=8)			::	T1(3),T2(2),T3
integer					::	i,j
	
	inix=x0
	
	x0=inix-h
	yi(1)=Func(i,j,L)
	
	x0=inix-h/2.0
	yi(2)=Func(i,j,L)
	
	x0=inix-h/4.0
	yi(3)=Func(i,j,L)
	
	x0=inix+h/4.0
	yi(4)=Func(i,j,L)
	
	x0=inix+h/2.0
	yi(5)=Func(i,j,L)
	
	x0=inix+h
	yi(6)=Func(i,j,L)
	
	T1(1)=(yi(4)-yi(3))*2.0/h
	T1(2)=(yi(5)-yi(2))/h
	T1(3)=(yi(6)-yi(1))/h/2.0
	
	T2(1)=(T1(1)-0.25*T1(2))/0.75
	T2(2)=(T1(2)-0.25*T1(3))/0.75
	T3	 =(T2(1)-0.0625*T2(2))/0.9375
	
	x0=inix
	
return
	
end subroutine