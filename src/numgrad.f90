subroutine numgrad(Func,sig,N)
use globe
implicit none
real (kind=8), external		::	func
real (kind=8), external		::	cbr_d,cbr_r,cbr_a
integer 					:: 	i,j,k,l,m,N    ! N the total number of iterations
real (kind=8)				::	T3
real (kind=8),intent(in)	::  sig                                                       
													   

N=0
do while(.true.)

	gd=0.0
	gr=0.0
	ga=0.0
	delta=0.0
	
	Do k=1,Nmole

	
		f =0.0
		fd=0.0
		fr=0.0
		fa=0.0
		
		
		do i=1,Nopt
			
			do l=1,p1(i,i,k)
			
				call deviation(Func,d(i,i),i,i,dist(i,i,l,k),T3)	
				fd(i)=fd(i)+T3
				call deviation(Func,r(i,i),i,i,dist(i,i,l,k),T3)
				fr(i)=fr(i)+T3
				call deviation(Func,a(i,i),i,i,dist(i,i,l,k),T3)
				fa(i)=fa(i)+T3
				
			end do
		
		end do
		

		do i=1,Npara
			
			do j=1,Npara
			
				do l=1,p1(i,j,k)
				
					f=f+    (   Func(i,j,dist(i,j,l,k)))
					
				end do
				
			end do
		
		end do
		
		
		
		gd=gd+(f-E(k))*fd*W(k)
		gr=gr+(f-E(k))*fr*W(k)
		ga=ga+(f-E(k))*fa*W(k)
		

		
	 	
	end do
	

	
	gd=gd/real(Nmole)
    gr=gr/real(Nmole)
    ga=ga/real(Nmole)
	

	do i=1,Nopt
		delta(2)=delta(2)+gd(i)*gd(i)+gr(i)*gr(i)+ga(i)*ga(i)
	end do	
	delta(2)=sqrt(delta(2))	
	if (delta(2)<=conver) then
		write(6,*)
		write(6,*)"The optimization was normally finished"
		write(6,*)
		exit
	end if
	
	do i=1,Nopt
	
		d(i,i)=d(i,i)-step*gd(i)
		r(i,i)=r(i,i)-step*gr(i)
		a(i,i)=a(i,i)-step*ga(i)
		
	end do
	
	do i=1,Npara-1
		do j=i+1,Npara
			
			d(i,j)=cbr_d(d(i,i),d(j,j))
			r(i,j)=cbr_r(r(i,i),r(j,j))
			a(i,j)=cbr_a(a(i,i),a(j,j),r(i,i),r(j,j),r(i,j))
			
			d(j,i)=d(i,j)
			r(j,i)=r(i,j)
	        a(j,i)=a(i,j)
			          
		end do        
	end do

!===========================================
! the total number of iterations
	N=N+1
!===========================================
	if(mod(N,100000)==0) then
		
	
		step=1.50*step
		
	else if (N>=1000000)then
		write(6,*)
		write(6,*)"The maximum number of iteration steps has been reached, &
		         the cycle will be ended, you should change input information"
		write(6,*)
		exit
		
	end if
		

end do

!==================================================
 

return

end subroutine numgrad