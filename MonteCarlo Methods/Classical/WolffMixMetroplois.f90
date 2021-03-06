module pre
implicit none
	type::spin
		integer::belong,positioni,positionj
		real*8::sx,sy,sz
	end type spin
	type::cluster
		type(spin),pointer::p
	end type cluster
	
	integer,parameter::l=10
	real*8,parameter::Jc=1.0,Kb=1.0,Dhundred=5.0,pie=3.1415926535898
	type(spin),target::s(l,l),pstart
	type(cluster)::member(l**2)
	integer*8::seed
	integer::dt(8),im(l),ip(l)
	real*8::D,E,M,mx,my,mz,rx,ry,rz,T
	real*8,external::random
contains
	subroutine seedchoice()
		call date_and_time(values=dt)
		seed=dt(1)+70*(dt(2)+12*(dt(3)+31*(dt(5)+23*(dt(6)+59*dt(7)))))
		return
	end subroutine
	
	subroutine chooseR()
	implicit none
!		real*8::temp1,temp2,temp
!		temp1=1-2*random(seed)
!		temp2=1-2*random(seed)
!		temp=temp1**2+temp2**2
!		do while(temp>=1)
!			temp1=1-2*random(seed)
!			temp2=1-2*random(seed)
!			temp=temp1**2+temp2**2
!		end do
!		rx=2*temp1*sqrt(1-temp)
!		ry=2*temp2*sqrt(1-temp)
!		rz=1-2*temp
		integer::temp
		temp=2*random(seed)
		if(temp==0) then
			rx=0
			ry=0
			rz=1
		else 
			rx=dcos(2*pie*random(seed))
			ry=sqrt(1-rx**2)
			rz=0
		end if
		return
	end subroutine

	function judge(i1,j1,i2,j2)
	implicit none
		integer::i1,j1,i2,j2,judge
		real*8::cosangle1,cosangle2
		cosangle1=s(i1,j1)%sx*rx+s(i1,j1)%sy*ry+s(i1,j1)%sz*rz
		cosangle2=s(i2,j2)%sx*rx+s(i2,j2)%sy*ry+s(i2,j2)%sz*rz
		if(s(i1,j1)%belong*s(i2,j2)%belong==0.and.random(seed)<1-exp(min(0.0,2*(Jc/(Kb*T))*cosangle1*cosangle2))) then
			judge=1
		else 
			judge=0
		end if
		return
	end function
			
	subroutine rotate(i,j)
	implicit none
		integer::i,j
		real*8::sumx,sumy,sumz,cosangle,flipsx,flipsy,flipsz
		sumx=s(im(i),j)%sx+s(ip(i),j)%sx+s(i,im(j))%sx+s(i,ip(j))%sx
		sumy=s(im(i),j)%sy+s(ip(i),j)%sy+s(i,im(j))%sy+s(i,ip(j))%sy
		sumz=s(im(i),j)%sz+s(ip(i),j)%sz+s(i,im(j))%sz+s(i,ip(j))%sz
		cosangle=s(i,j)%sx*rx+s(i,j)%sy*ry+s(i,j)%sz*rz
		flipsx=s(i,j)%sx-2*cosangle*rx
		flipsy=s(i,j)%sy-2*cosangle*ry
		flipsz=s(i,j)%sz-2*cosangle*rz
		E=E-Jc*((flipsx-s(i,j)%sx)*sumx+(flipsy-s(i,j)%sy)*sumy+(flipsz-s(i,j)%sz)*sumz)
		mx=mx+flipsx-s(i,j)%sx
		my=my+flipsy-s(i,j)%sy
		mz=mz+flipsz-s(i,j)%sz
		s(i,j)%sx=flipsx
		s(i,j)%sy=flipsy
		s(i,j)%sz=flipsz
		s(i,j)%belong=1
		return
	end subroutine

	subroutine metripolisMCS()
	implicit none
		integer::i,j
		real*8::temp1,temp2,temp,deltaE,sxtrial,sytrial,sztrial,sumsx,sumsy,sumsz
		do i=1,l
			do j=1,l
				temp1=1-2*random(seed)
				temp2=1-2*random(seed)
				temp=temp1**2+temp2**2
				do while(temp>=1)
					temp1=1-2*random(seed)
					temp2=1-2*random(seed)
					temp=temp1**2+temp2**2
				end do
				sxtrial=2*temp1*sqrt(1-temp)
				sytrial=2*temp2*sqrt(1-temp)
				sztrial=1-2*temp
				sumsx=s(im(i),j)%sx+s(ip(i),j)%sx+s(i,im(j))%sx+s(i,ip(j))%sx
				sumsy=s(im(i),j)%sy+s(ip(i),j)%sy+s(i,im(j))%sy+s(i,ip(j))%sy
				sumsz=s(im(i),j)%sz+s(ip(i),j)%sz+s(i,im(j))%sz+s(i,ip(j))%sz
				deltaE=-Jc*((sxtrial-s(i,j)%sx)*sumsx+(sytrial-s(i,j)%sy)*sumsy+(sztrial-s(i,j)%sz)*sumsz)-D*(sztrial**2-s(i,j)%sz**2)
				if(random(seed)<min(1d0,exp(-deltaE/(Kb*T)))) then
					E=E+deltaE
					mx=mx+sxtrial-s(i,j)%sx
					my=my+sytrial-s(i,j)%sy
					mz=mz+sztrial-s(i,j)%sz
					s(i,j)%sx=sxtrial
					s(i,j)%sy=sytrial
					s(i,j)%sz=sztrial
				end if
			end do
		end do
		return
	end subroutine

end module


program main
use pre
implicit none
	integer*8,parameter::totalsteps=1200000,balancesteps=200000
	integer*8::isweep,n,first,last,last2,count,num
	integer::i,j,Tcycle
	real*8::E1,E2,M1,M2,M4,c,kai,Mav,U
	open(unit=10,file='l10WolffMixMetripolis.txt')
	D=Dhundred/100
	
	call seedchoice()
	
	do i=1,l
		im(i)=i-1
		ip(i)=i+1
	end do
	im(1)=l
	ip(l)=1
		
	do i=1,l
		do j=1,l
			s(i,j)%sx=0
			s(i,j)%sy=0
			s(i,j)%sz=-1
			s(i,j)%belong=0
			s(i,j)%positioni=i
			s(i,j)%positionj=j
		end do
	end do

	do n=1,l**2
		member(n)%p=>pstart
	end do
	
	E=0
	mx=0
	my=0
	mz=0
	do i=1,l
		do j=1,l
			mx=mx+s(i,j)%sx
			my=my+s(i,j)%sy
			mz=mz+s(i,j)%sz
			E=E-Jc*(s(i,j)%sx*(s(im(i),j)%sx+s(i,im(j))%sx)+&
					s(i,j)%sy*(s(im(i),j)%sy+s(i,im(j))%sy)+&
					s(i,j)%sz*(s(im(i),j)%sz+s(i,im(j))%sz))-D*s(i,j)%sz**2	
		end do
	end do
	M=sqrt(mx**2+my**2+mz**2)
	write(*,*) E,M,D
	
	Tcycle=5
	do while(Tcycle<=155)
		T=Tcycle/100.0
		count=0
		num=0
		E1=0
		E2=0
		M1=0
		M2=0
		M4=0
		do isweep=1,totalsteps
			if(mod(isweep,2)==0) then
				call chooseR()
				i=random(seed)*l+1
				j=random(seed)*l+1
				call rotate(i,j)
				n=1
				member(n)%p=>s(i,j)
				first=n
				last=n
				last2=last

		100		do n=first,last
					i=member(n)%p%positioni
					j=member(n)%p%positionj
					if(judge(i,j,im(i),j)==1) then
						call rotate(im(i),j)
						last2=last2+1
						member(last2)%p=>s(im(i),j)
					end if
					if(judge(i,j,ip(i),j)==1) then
						call rotate(ip(i),j)
						last2=last2+1
						member(last2)%p=>s(ip(i),j)
					end if
					if(judge(i,j,i,im(j))==1) then
						call rotate(i,im(j))
						last2=last2+1
						member(last2)%p=>s(i,im(j))
					end if
					if(judge(i,j,i,ip(j))==1) then
						call rotate(i,ip(j))
						last2=last2+1
						member(last2)%p=>s(i,ip(j))
					end if
				end do

		!		write(*,*) first,last,last2
				if(last2>last) then
					first=last+1
					last=last2
					goto 100
				else
					do n=1,last2
						member(n)%p%belong=0
						member(n)%p=>pstart
					end do
				end if
			else
				call metripolisMCS()
			end if

			if(isweep>balancesteps) then
				count=count+1
				if(mod(count,5)==0) then
					num=num+1		
					E1=E1+E
					E2=E2+E**2
					M=sqrt(mx**2+my**2+mz**2)
					M1=M1+M
					M2=M2+M**2
					M4=M4+M**4
				end if
			end if
		end do
!		write(*,*) count,num		

		c=((E2/num)-(E1/num)**2)/(Kb*(T**2)*(L**2))
		kai=((M2/num)-(M1/num)**2)/(Kb*T*(L**2))
		Mav=M1/(num*L**2)
		U=1-M4*num/(3*M2**2)
		write(10,"(e15.6,e15.6,e15.6,e15.6,e15.6)") T,c,kai,Mav,U
		write(*,"(e15.6,e15.6,e15.6,e15.6,e15.6)") T,c,kai,Mav,U
		
		if(-10<=Tcycle-60.and.Tcycle-60<10) then
			Tcycle=Tcycle+1
		else
			Tcycle=Tcycle+5
		end if
	end do
stop
end

real*8 function random(ibm)
	integer*8 mult,modulo,ibm
	real*8 rmodulo
	mult=16807
	modulo=2147483647
	rmodulo=modulo
	ibm=ibm*mult
	ibm=mod(ibm,modulo)
	random=ibm/rmodulo
return
end

