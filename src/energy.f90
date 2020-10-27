subroutine energy(Func,ifEQ)
use globe
implicit none
real (kind=8), external		::	func
integer 					:: 	i,j,k,l,m    ! N the total number of iterations
logical						::  ifEQ
                                                      												

delta(1)=0.0

if (ifEQ) then
	write(6,"(4(3xA))") "Etotal-Ecoulomb","Evdw","Etotal","Evdw+EQ(k)"
else 
	write(6,"(2(3xA))")"Etotal","Evdw"
end if 

Do k=1,Nmole

	f =0.0

	do i=1,Npara
		
		do j=1,Npara
		
			do l=1,p1(i,j,k)
			
				f=f+    (   Func(i,j,dist(i,j,l,k)))
				
			end do
			
		end do
	
	end do
	
	
	if (ifEQ) then
		write(6,*) E(k),f,E(k)+EQ(k),f+EQ(k)
	else 
		write(6,*)E(k),f
	end if 
	
	
	delta(1)=delta(1)+((f-E(k))**2)

	
 	
end do



delta(1)=delta(1)/real(Nmole)

	

	
	

!===========================================
! the total number of iterations


!==================================================

return

end subroutine energy