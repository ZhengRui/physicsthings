program main
implicit none
integer*8 L,i,j,isweep,totalsteps,balancesteps,ici,ien,count,n,seed
parameter(L=30,totalsteps=60000,balancesteps=20000)
integer*8 S(1:L,1:L),ip(1:L),im(1:L)
integer dt(8)
real*8 E,E1,E2,M,M1,M2,W(-4:4),Jc,K,T,x,c,kai,Mav,M4,U
real,external::random
parameter(Jc=1,K=1)
open(unit=10,file='L=30.txt')
	
	call date_and_time(values=dt)
	seed=dt(1)+70*(dt(2)+12*(dt(3)+31*(dt(5)+23*(dt(6)+59*dt(7)))))

	do i=1,L
		ip(i)=i+1
		im(i)=i-1
	end do
	ip(L)=1
	im(1)=L

	do i=1,L
		do j=1,L
			S(i,j)=-1
		end do
	end do

	E=0
	M=0
	do i=1,L
		do j=1,L
			ici=S(i,j)
			ien=S(ip(i),j)+S(im(i),j)+S(i,ip(j))+S(i,im(j))
			ien=ici*ien
			E=E-Jc*ien
			M=M+S(i,j)
		end do
	end do
	E=E/2

	do T=2.0,2.41,0.02
		do i=-4,4,2
			W(i)=1
			if(i>0) then
				W(i)=exp(-2*Jc*i/(K*T))
			end if
		end do
		count=0
		n=0
		E1=0
		E2=0
		M1=0
		M2=0
		M4=0
		do isweep=1,totalsteps
			do i=1,L
				do j=1,L
					ici=S(i,j)
					ien=S(ip(i),j)+S(im(i),j)+S(i,ip(j))+S(i,im(j))
					ien=ici*ien					
					x=random(seed)
					if(x<W(ien)) then
						S(i,j)=-ici
						E=E+2*Jc*ien
						M=M-2*ici
					end if					
				end do
			end do
			if(isweep>balancesteps) then
						count=count+1
						if(mod(count,5)==0) then
							n=n+1
							E1=E1+E
							E2=E2+E**2
							M1=M1+abs(M)
							M2=M2+M**2
							M4=M4+M**4
						end if
			end if
		end do

		c=((E2/n)-(E1/n)**2)/(K*(T**2)*(L**2))
		kai=((M2/n)-(M1/n)**2)/(K*T*(L**2))
		Mav=M1/(n*L**2)
		U=1-M4*n/(3*M2**2)
		write(10,"(e15.6,e15.6,e15.6,e15.6,e15.6)") T,c,kai,Mav,U
		write(*,"(e15.6,e15.6,e15.6,e15.6,e15.6)") T,c,kai,Mav,U
	end do

close(10)
stop
end

real function random(ibm)
	integer*8 mult,modulo,ibm
	real rmodulo
	mult=16807
	modulo=2147483647
	rmodulo=modulo
	ibm=ibm*mult
	ibm=mod(ibm,modulo)
	random=ibm/rmodulo
	return
end