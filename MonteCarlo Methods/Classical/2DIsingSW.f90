module typedef
implicit none
	type::spin
		integer::s
		integer::bond(2)
		type(spin),pointer::pnext
		type(cluster),pointer::pbelong
	end type spin

	type::cluster
		integer::nc,sc
		type(spin),pointer::phead,ptail
	end type cluster
end module typedef

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

subroutine create(newcluster,startposition)
use typedef
implicit none
	type(cluster),target::newcluster
	type(spin),target::startposition
	newcluster%sc=startposition%s
	newcluster%phead=>startposition
	newcluster%ptail=>startposition
	startposition%pbelong=>newcluster
	return
end 

subroutine addpoint(existcluster,point)
use typedef
implicit none
	type(cluster),target::existcluster
	type(spin),target::point
	existcluster%ptail%pnext=>point
	point%pbelong=>existcluster
	existcluster%ptail=>point
	return
end 

!下面函数用来连接两个本质集团，不过主程序要改动一下
!subroutine combine(cluster1,cluster2)
!use typedef
!implicit none
!	type(cluster)::cluster1,cluster2
!	type(spin),pointer::temp
!	cluster1%ptail%pnext=>cluster2%phead
!	cluster1%ptail=>cluster2%ptail
!	temp=>cluster2%phead
!	do while(.not.associated(temp,cluster2%ptail))
!		temp%pbelong=>cluster1
!		temp=>temp%pnext
!	end do
!	temp%pbelong=>cluster1
!	cluster2%phead=>pheadstart
!	cluster2%ptail=>ptailstart
!	cluster%sc=0
!	return
!end 

program main
use typedef
implicit none
		
    integer,parameter::L=50,totalsteps=50000,balancesteps=10000
	real,parameter::Jc=1,Kb=1
	integer*8::i,j,n,seed,lnc,unc,rnc,dnc,step,count,nsample
	integer::minus(L),plus(L),dt(8)
	real*8::T,E,M,E1,E2,M1,M2,M4,sheat,kai,Mav,U
	type(spin),pointer::temp,temp2
	type(spin),target::sp(L,L),pheadstart,ptailstart
	type(cluster),target::c(L**2),pbelongstart
	real,external::random
	open(unit=10,file='L=50.txt')

	call date_and_time(values=dt)
	seed=dt(1)+70*(dt(2)+12*(dt(3)+31*(dt(5)+23*(dt(6)+59*dt(7)))))

	do i=1,L
		plus(i)=i+1
		minus(i)=i-1
	end do
	plus(L)=1
	minus(1)=L
	
	do i=1,L
		do j=1,L
			sp(i,j)%s=-1
			sp(i,j)%bond(1)=0
			sp(i,j)%bond(2)=0
			sp(i,j)%pnext=>sp(i,j)
			sp(i,j)%pbelong=>pbelongstart
		end do
	end do

	do n=1,L**2
		c(n)%sc=0
		c(n)%phead=>pheadstart
		c(n)%ptail=>ptailstart
		c(n)%nc=n
	end do

	do T=2.20,2.35,0.01
		count=0
		nsample=0
		E1=0
		E2=0
		M1=0
		M2=0
		M4=0
		
		do step=1,totalsteps	
			E=0
			M=0		
			do i=1,L
				do j=1,L
					E=E-Jc*sp(i,j)%s*(sp(i,minus(j))%s+sp(minus(i),j)%s+sp(i,plus(j))%s+sp(plus(i),j)%s)
					M=M+sp(i,j)%s
					if(sp(i,j)%s==sp(i,minus(j))%s) then
						if(random(seed)<1-exp(-2*Jc/(Kb*T))) then
							sp(i,j)%bond(1)=1
						end if
					end if
					if(sp(i,j)%s==sp(minus(i),j)%s) then
						if(random(seed)<1-exp(-2*Jc/(Kb*T))) then
							sp(i,j)%bond(2)=1
						end if
					end if
				end do
			end do
			E=E/2
		
			if(step>balancesteps) then
				count=count+1
				if(mod(count,5)==0) then
					nsample=nsample+1
					E1=E1+E
					E2=E2+E**2
					M1=M1+abs(M)
					M2=M2+M**2
					M4=M4+M**4
				end if
			end if

		!	write(*,*) E,M
		!	do i=1,L
		!		do j=1,L
		!			write(*,*) associated(sp(i,j)%pnext,sp(i,j)),associated(sp(i,j)%pbelong,pbelongstart)
		!			write(*,*) sp(i,j)%s,sp(i,j)%bond(1),sp(i,j)%bond(2)
		!		end do
		!	end do

! HK算法变种
			n=0
			do i=1,L
				do j=1,L
					if(sp(i,j)%bond(1)==0.and.sp(i,j)%bond(2)==0) then	
						n=n+1
						call create(c(n),sp(i,j))
					else if(sp(i,j)%bond(1)==1.and.sp(i,j)%bond(2)==0) then
						if(associated(sp(i,minus(j))%pbelong,pbelongstart)) then
							n=n+1
							call create(c(n),sp(i,j))
						else
							call addpoint(sp(i,minus(j))%pbelong,sp(i,j))
						end if
					else if(sp(i,j)%bond(1)==0.and.sp(i,j)%bond(2)==1) then
						if(associated(sp(minus(i),j)%pbelong,pbelongstart)) then
							n=n+1
							call create(c(n),sp(i,j))
						else
							call addpoint(sp(minus(i),j)%pbelong,sp(i,j))
						end if
					else
						if(associated(sp(i,minus(j))%pbelong,pbelongstart)) then
							if(associated(sp(minus(i),j)%pbelong,pbelongstart)) then
								n=n+1
								call create(c(n),sp(i,j))
							else
								call addpoint(sp(minus(i),j)%pbelong,sp(i,j))
							end if
						else
							if(associated(sp(minus(i),j)%pbelong,pbelongstart)) then
								call addpoint(sp(i,minus(j))%pbelong,sp(i,j))
							else
								lnc=sp(i,minus(j))%pbelong%nc
								do while(lnc<0)
									lnc=c(-lnc)%nc
								end do
								unc=sp(minus(i),j)%pbelong%nc
								do while(unc<0)
									unc=c(-unc)%nc
								end do
								if(lnc<unc) then
									c(unc)%nc=-lnc
									call addpoint(c(lnc),sp(i,j))
								else if(unc<lnc) then
									c(lnc)%nc=-unc
									call addpoint(c(unc),sp(i,j))
								else
									call addpoint(c(lnc),sp(i,j))
								end if
							end if
						end if
					end if
				end do
			end do
			
!加上下面两段，实现周期边界条件
			do i=1,L
				if(sp(i,1)%bond(1)==1) then
					lnc=sp(i,1)%pbelong%nc
					do while(lnc<0)
						lnc=c(-lnc)%nc
					end do
					rnc=sp(i,L)%pbelong%nc
					do while(rnc<0)
						rnc=c(-rnc)%nc
					end do
					if(lnc<rnc) then
						c(rnc)%nc=-lnc
					else if(rnc<lnc) then
						c(lnc)%nc=-rnc
					end if
				end if
			end do

			do j=1,L
				if(sp(1,j)%bond(2)==1) then
					unc=sp(1,j)%pbelong%nc
					do while(unc<0)
						unc=c(-unc)%nc
					end do
					dnc=sp(L,j)%pbelong%nc
					do while(dnc<0)
						dnc=c(-dnc)%nc
					end do
					if(unc<dnc) then
						c(dnc)%nc=-unc
					else if(dnc<unc) then
						c(unc)%nc=-dnc
					end if
				end if
			end do


		!	write(*,*) 'test'
		!	do i=1,L
		!		do j=1,L
		!			write(*,*)	sp(i,j)%s,sp(i,j)%pbelong%sc,sp(i,j)%pbelong%nc,'shit'	
		!		end do
		!	end do

		!	write(*,*) 'test'
		!	do n=1,L**2
		!		write(*,*) n,c(n)%sc,c(n)%nc
		!	end do

			do n=1,L**2
				if(c(n)%sc==0) then
					exit
				else
					if(c(n)%nc>0) then
						if(random(seed)<0.5) then
							c(n)%sc=1
						else 
							c(n)%sc=-1
						end if
					else
						c(n)%sc=c(-c(n)%nc)%sc
						c(n)%nc=n
					end if
					temp=>c(n)%phead
					do while(.not.associated(temp,c(n)%ptail))
						temp%s=c(n)%sc
						temp%bond(1)=0
						temp%bond(2)=0
						temp%pbelong=>pbelongstart
						temp2=>temp%pnext
						temp%pnext=>temp
						temp=>temp2
					end do
					temp%s=c(n)%sc
					temp%bond(1)=0
					temp%bond(2)=0
					temp%pbelong=>pbelongstart
				end if
			end do

			do n=1,L**2
				if(c(n)%sc==0) then
					exit
				else
					c(n)%sc=0
					c(n)%phead=>pheadstart
					c(n)%ptail=>ptailstart	
				end if
			end do

		!	do i=1,L
		!		do j=1,L
		!			write(*,*) sp(i,j)%s,sp(i,j)%bond(1),sp(i,j)%bond(2),&
		!			associated(sp(i,j)%pnext,sp(i,j)),associated(sp(i,j)%pbelong,pbelongstart)
		!		end do
		!	end do
					
		!	do n=1,L**2
		!		write(*,*) c(n)%sc,c(n)%nc,associated(c(n)%phead,pheadstart),&
		!		associated(c(n)%ptail,ptailstart)		
		!	end do	
		end do

		sheat=((E2/nsample)-(E1/nsample)**2)/(Kb*(T**2)*(L**2))
		kai=((M2/nsample)-(M1/nsample)**2)/(Kb*T*(L**2))
		Mav=M1/(nsample*L**2)
		U=1-M4*nsample/(3*M2**2)
		write(10,"(e15.6,e15.6,e15.6,e15.6,e15.6)") T,sheat,kai,Mav,U
		write(*,"(e15.6,e15.6,e15.6,e15.6,e15.6)") T,sheat,kai,Mav,U
	end do


close(10)
stop
end
