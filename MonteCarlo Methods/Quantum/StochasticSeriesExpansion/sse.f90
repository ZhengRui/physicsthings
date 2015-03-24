module pre
implicit none
	type::spin
		integer::position(3)
		integer::sz,positioninsN
	end type spin

	type::bond
		type(spin),pointer::b1spin,b2spin
	end type bond

	type::operator
		type(bond),pointer::linkb
		integer::opttype
	end type operator

	type::extendedvertex
		type(spin)::leg(4)
		integer::evtxtype
	end type extendedvertex

	type::heatbath
		real*8,allocatable::cumuprob(:)
		integer,allocatable::linkedleg(:)
		integer::nexit		
	end type heatbath
	
	type(spin),target,allocatable::s(:)
	type(bond),target,allocatable::b(:)
	type(operator),allocatable::opt(:),vertex(:)
	type(extendedvertex),allocatable::evtx(:)
	type(bond),target::b0
	type(heatbath),allocatable::hbprob(:,:)
	integer,allocatable::link(:),first(:),last(:)
	integer*8::seed
	integer::dt(8),d,N,Nb,Lx,Ly,Lz,optM,optN,optMnew,optMold,Nloop,tloopl,new(6,4,4),grow
	real*8::J,T,aniso,H,e,hb,C0,C,beta,evtxp(6),temp(6),devtxptable(6,4,4),legchooseprob(6,4,4),&
			e1a,e1b,e1c,e1b1,e1b2,e1b3,e2a,e2b,e2c,e2b1,e2b2,e2b3,&
			energyav,specificheatav,uniformedabsMav,staggeredabsMav,&
			uniformedstructuref,staggeredstructuref,uniformedsusav,staggeredsusav	
	real*8,external::random

contains
	
	subroutine seedchoice()
		call date_and_time(values=dt)
		seed=dt(1)+70*(dt(2)+12*(dt(3)+31*(dt(5)+23*(dt(6)+59*dt(7)))))
		return
	end subroutine

	subroutine makelattice()
		integer::i,j,k
		write(*,*) "please input the dimension:"
		read(*,*) d
		if(d==1) then
			write(*,*) "please input the lattice size Lx:"
			read(*,*) Lx
			N=Lx
		else if(d==2) then
			write(*,*) "please input the lattice size Lx,Ly:"
			read(*,*) Lx,Ly
			N=Lx*Ly
		else if(d==3) then
			write(*,*) "please input the lattice size Lx,Ly,Lz:"
			read(*,*) Lx,Ly,Lz
			N=Lx*Ly*Lz
		else 
			write(*,*) "dimension is larger than 3, similar codes like d<=3 should be added"
			stop
		end if		
			Nb=d*N
			allocate(s(N))
			allocate(b(Nb))
			
		if(d==1) then
			do i=1,Lx
				s(i)%position(1)=i
				s(i)%position(2)=0
				s(i)%position(3)=0
				b(i)%b1spin=>s(i)
				b(i)%b2spin=>s(mod(i,Lx)+1)
			end do
		else if(d==2) then
			do i=1,Lx
				do j=1,Ly
					s((j-1)*Lx+i)%position(1)=i
					s((j-1)*Lx+i)%position(2)=j
					s((j-1)*Lx+i)%position(3)=0
					b((j-1)*Lx+i)%b1spin=>s((j-1)*Lx+i)
					b((j-1)*Lx+i)%b2spin=>s((j-1)*Lx+mod(i,Lx)+1)
					b(N+(j-1)*Lx+i)%b1spin=>s((j-1)*Lx+i)
					b(N+(j-1)*Lx+i)%b2spin=>s(mod(j,Ly)*Lx+i)
				end do
			end do
		else 
			do i=1,Lx
				do j=1,Ly
					do k=1,Lz
						s((k-1)*Lx*Ly+(j-1)*Lx+i)%position(1)=i
						s((k-1)*Lx*Ly+(j-1)*Lx+i)%position(2)=j
						s((k-1)*Lx*Ly+(j-1)*Lx+i)%position(3)=k
						b((k-1)*Lx*Ly+(j-1)*Lx+i)%b1spin=>s((k-1)*Lx*Ly+(j-1)*Lx+i)
						b((k-1)*Lx*Ly+(j-1)*Lx+i)%b2spin=>s((k-1)*Lx*Ly+(j-1)*Lx+mod(i,Lx)+1)
						b(N+(k-1)*Lx*Ly+(j-1)*Lx+i)%b1spin=>s((k-1)*Lx*Ly+(j-1)*Lx+i)
						b(N+(k-1)*Lx*Ly+(j-1)*Lx+i)%b2spin=>s((k-1)*Lx*Ly+mod(j,Ly)*Lx+i)
						b(2*N+(k-1)*Lx*Ly+(j-1)*Lx+i)%b1spin=>s((k-1)*Lx*Ly+(j-1)*Lx+i)
						b(2*N+(k-1)*Lx*Ly+(j-1)*Lx+i)%b2spin=>s(mod(k,Lz)*Lx*Ly+(j-1)*Lx+i)
					end do
				end do
			end do		
		end if
		return
	end subroutine	

	subroutine parametersinitialize()
		write(*,*) "please input J, T, uniaxial anisotropic parameter, and external positive magnetic field H:"
		read(*,*) J,T,aniso,H
		hb=H/(2*d*J)
		if(2*hb<=(1-aniso).and.2*hb<=(1+aniso)) then
			write(*,*) "hb=",hb,"aniso=",aniso,"the simulation is in the region 1 of the algorithmic phase diagram, emin=", (1-aniso)/4-hb/2
			write(*,*) "please input the value of e:"
			read(*,*) e
			e1b1=0; e1b2=0; e1b3=0; e2b1=0;	e2b2=0; e2b3=0
		else if(2*hb<(aniso-1)) then
			write(*,*) "hb=",hb,"aniso=",aniso,"the simulation is in the region 2 of the algorithmic phase diagram, emin=", 0
			write(*,*) "please input the value of e:"
			read(*,*) e
			e1b1=0; e1b2=hb-(1-aniso)/2; e1b3=0; e2b1=0; e2b2=-hb-(1-aniso)/2; e2b3=0
		else if(2*hb>(1-aniso).and.2*hb<=(1+aniso)) then
			write(*,*) "hb=",hb,"aniso=",aniso,"the simulation is in the region 3 of the algorithmic phase diagram, emin=", 0
			write(*,*) "please input the value of e:"
			read(*,*) e
			e1b1=0; e1b2=hb-(1-aniso)/2; e1b3=0; e2b1=0; e2b2=0; e2b3=0
		else if(2*hb>(1-aniso).and.2*hb>(1+aniso)) then
			write(*,*) "hb=",hb,"aniso=",aniso,"the simulation is in the region 4 of the algorithmic phase diagram, emin=", 0
			write(*,*) "please input the value of e:"
			read(*,*) e
			e1b1=0; e1b2=hb-(1-aniso)/2; e1b3=0; e2b1=0; e2b2=0; e2b3=hb-(1+aniso)/2
		else if(2*hb<-(1+aniso)) then
			write(*,*) "hb=",hb,"aniso=",aniso,"the simulation is in the region 6 of the algorithmic phase diagram, emin=", -hb-aniso/2
			write(*,*) "please input the value of e:"
			read(*,*) e
			e1b1=0; e1b2=0; e1b3=-hb-(1+aniso)/2; e2b1=0; e2b2=0; e2b3=hb-(1+aniso)/2
		else 
			write(*,*) "hb=",hb,"aniso=",aniso,"the simulation is in the region 5 of the algorithmic phase diagram, emin=", (1-aniso)/4-hb/2
			write(*,*) "please input the value of e:"
			read(*,*) e
			e1b1=0; e1b2=0; e1b3=0; e2b1=0; e2b2=0; e2b3=hb-(1+aniso)/2
		end if
		e1a=(1+aniso+2*(hb-e1b1-e1b2+e1b3))/4
		e1b=(1-aniso-2*(hb+e1b1-e1b2+e1b3))/4
		e1c=(aniso-1+2*(hb+e1b1-e1b2-e1b3)+4*e)/4
		e2a=(1+aniso-2*(hb+e2b1+e2b2-e2b3))/4
		e2b=(1-aniso+2*(hb-e2b1+e2b2-e2b3))/4
		e2c=(aniso-1+6*hb+4*e+2*(e2b1-e2b2-e2b3))/4
		C0=aniso/4+hb
		C=C0+e
		beta=J/T
		evtxp(1)=e
		evtxp(2)=aniso/2+hb+e
		evtxp(3)=evtxp(2)
		evtxp(4)=5d-1
		evtxp(5)=evtxp(4)
		evtxp(6)=e+2*hb
		temp=Nb*beta*evtxp
		return
	end subroutine

	subroutine makeprobtable()
		integer::i,j,k,l(4),table(3,12),nonzero
		real*8::sum,temp(4),cuprob
		new=0
		do i=1,6
			do j=1,4
				do k=1,4
					if(i==1) then 
						l=-1
					else if(i==6) then 
						l=1
					else if(i==2) then 
						l(1)=-1; l(2)=1; l(3)=1; l(4)=-1
					else if(i==3) then
						l(1)=1; l(2)=-1; l(3)=-1; l(4)=1
					else if(i==4) then
						l(1)=-1; l(2)=1; l(3)=-1; l(4)=1
					else
						l(1)=1; l(2)=-1; l(3)=1; l(4)=-1
					end if
					l(j)=-l(j)
					l(k)=-l(k)
					if(l(1)==-1.and.l(2)==-1.and.l(3)==-1.and.l(4)==-1) then 
						new(i,j,k)=1
					else if(l(1)==1.and.l(2)==1.and.l(3)==1.and.l(4)==1) then 
						new(i,j,k)=6
					else if(l(1)==-1.and.l(2)==1.and.l(3)==1.and.l(4)==-1) then
						new(i,j,k)=2
					else if(l(1)==1.and.l(2)==-1.and.l(3)==-1.and.l(4)==1) then
						new(i,j,k)=3
					else if(l(1)==-1.and.l(2)==1.and.l(3)==-1.and.l(4)==1) then
                		                new(i,j,k)=4
					else if(l(1)==1.and.l(2)==-1.and.l(3)==1.and.l(4)==-1) then
                       			        new(i,j,k)=5
					else
					end if		
!					write(*,*) "new(",i,j,k,")=",new(i,j,k)	
				end do
			end do
		end do
		devtxptable=0
		legchooseprob=0
		open(unit=10, file='graphtable1.in', status='old')
		read(10,*) table
		do i=1,12
			if(mod(i,3)==1) then
				devtxptable(table(1,i)/100,mod(table(1,i),100)/10,mod(table(1,i),10))=e1b1
				devtxptable(table(2,i)/100,mod(table(2,i),100)/10,mod(table(2,i),10))=e1a
				devtxptable(table(3,i)/100,mod(table(3,i),100)/10,mod(table(3,i),10))=e1b
			else if(mod(i,3)==2) then
				devtxptable(table(1,i)/100,mod(table(1,i),100)/10,mod(table(1,i),10))=e1a
				devtxptable(table(2,i)/100,mod(table(2,i),100)/10,mod(table(2,i),10))=e1b2
				devtxptable(table(3,i)/100,mod(table(3,i),100)/10,mod(table(3,i),10))=e1c
			else
				devtxptable(table(1,i)/100,mod(table(1,i),100)/10,mod(table(1,i),10))=e1b
				devtxptable(table(2,i)/100,mod(table(2,i),100)/10,mod(table(2,i),10))=e1c
				devtxptable(table(3,i)/100,mod(table(3,i),100)/10,mod(table(3,i),10))=e1b3
			end if			
		end do
		open(unit=20, file='graphtable2.in', status='old')
		read(20,*) table
		do i=1,12
			if(mod(i,3)==1) then
				devtxptable(table(1,i)/100,mod(table(1,i),100)/10,mod(table(1,i),10))=e2b1
				devtxptable(table(2,i)/100,mod(table(2,i),100)/10,mod(table(2,i),10))=e2a
				devtxptable(table(3,i)/100,mod(table(3,i),100)/10,mod(table(3,i),10))=e2b
			else if(mod(i,3)==2) then
				devtxptable(table(1,i)/100,mod(table(1,i),100)/10,mod(table(1,i),10))=e2a
				devtxptable(table(2,i)/100,mod(table(2,i),100)/10,mod(table(2,i),10))=e2b2
				devtxptable(table(3,i)/100,mod(table(3,i),100)/10,mod(table(3,i),10))=e2c
			else
				devtxptable(table(1,i)/100,mod(table(1,i),100)/10,mod(table(1,i),10))=e2b
				devtxptable(table(2,i)/100,mod(table(2,i),100)/10,mod(table(2,i),10))=e2c
				devtxptable(table(3,i)/100,mod(table(3,i),100)/10,mod(table(3,i),10))=e2b3
			end if			
		end do
				
		do i=1,6
			if(evtxp(i)/=0) then
				do j=1,4
					temp=0
					sum=0
					do k=1,4
						if(abs(devtxptable(i,j,k))>1d-15) then
							temp(k)=devtxptable(i,j,k)
							sum=sum+temp(k)
						end if
					end do
				legchooseprob(i,j,:)=temp/sum
				end do
			end if							
		end do
		
		allocate(hbprob(6,4))
		do i=1,6
			if(evtxp(i)/=0) then
				do j=1,4
					nonzero=0
					do k=1,4
						if(legchooseprob(i,j,k)/=0) then
							nonzero=nonzero+1
						end if
					end do
					allocate(hbprob(i,j)%cumuprob(nonzero))
					allocate(hbprob(i,j)%linkedleg(nonzero))
					hbprob(i,j)%nexit=nonzero
					nonzero=0
					cuprob=0
					do k=1,4
						if(legchooseprob(i,j,k)/=0) then
							cuprob=cuprob+legchooseprob(i,j,k)
							nonzero=nonzero+1
							hbprob(i,j)%cumuprob(nonzero)=cuprob
							hbprob(i,j)%linkedleg(nonzero)=k
						end if
					end do
					hbprob(i,j)%cumuprob(nonzero)=1
!					write(*,*) "vertex=",i,"entrance=",j,"nexit=",hbprob(i,j)%nexit,"exit=",hbprob(i,j)%linkedleg,"cumuprob=",hbprob(i,j)%cumuprob
				end do
			end if
		end do

		return
	end subroutine

	subroutine configinitialize()
		integer::i
		do i=1,N
			s(i)%sz=2*int(2.*random(seed))-1
			s(i)%positioninsN=i
		end do				
		write(*,*) "please input the initial operator string length:"
		read(*,*) optM
		allocate(opt(optM))
		do i=1,optM
			opt(i)%opttype=0
			opt(i)%linkb=>b0
		end do	
		optN=0
		return
	end subroutine

	subroutine Dupdate()
		integer::i,j,k
		do i=1,optM
			if(opt(i)%opttype==0) then
				j=Nb*random(seed)+1
				if(b(j)%b1spin%sz/=b(j)%b2spin%sz) then 
					k=2
				else if(b(j)%b1spin%sz==-1) then 
					k=1
				else 
					k=6
				end if
				if(temp(k)>=dfloat(optM-optN).or.temp(k)>=random(seed)*(optM-optN)) then
						opt(i)%linkb=>b(j)
						opt(i)%opttype=1
						optN=optN+1
				end if
			else if(opt(i)%opttype==1) then
				if(opt(i)%linkb%b1spin%sz/=opt(i)%linkb%b2spin%sz) then 
					k=2
				else if(opt(i)%linkb%b1spin%sz==-1) then 
					k=1
				else 
					k=6
				end if
				if(temp(k)<=dfloat(optM-optN+1).or.temp(k)<=(optM-optN+1)/random(seed)) then
						opt(i)%linkb=>b0						
						opt(i)%opttype=0
						optN=optN-1
				end if
			else
				opt(i)%linkb%b1spin%sz=-opt(i)%linkb%b1spin%sz
				opt(i)%linkb%b2spin%sz=-opt(i)%linkb%b2spin%sz
			end if
		end do
		return
	end subroutine

	subroutine directedloopupdate()
		integer::i,j,k,k2,leg0,l,loop,loopl
		real*8::r
!		write(*,*) "++after Dupdate+++++++++++++++++++++++++++++, optN=",optN
		if(optN==0) return		
		allocate(vertex(optN))
		allocate(evtx(optN))
		allocate(link(4*optN))
		allocate(first(N))
		allocate(last(N))	

		j=1
		do i=1,optM
			if(opt(i)%opttype/=0) then
				vertex(j)=opt(i)
				j=j+1		
			end if	
		end do
		do i=1,N
			first(i)=0
			last(i)=0
		end do
		do i=1,4*optN
			link(i)=0
		end do

		do i=1,optN
			j=vertex(i)%linkb%b1spin%positioninsN
			evtx(i)%leg(1)=s(j)
			if(last(j)==0) then
				first(j)=4*i-3
				last(j)=4*i
			else
				link(4*i-3)=last(j)
				link(last(j))=4*i-3
				last(j)=4*i
			end if
			j=vertex(i)%linkb%b2spin%positioninsN
			evtx(i)%leg(2)=s(j)
			if(last(j)==0) then
				first(j)=4*i-2
				last(j)=4*i-1
			else
				link(4*i-2)=last(j)
				link(last(j))=4*i-2
				last(j)=4*i-1
			end if	
			if(vertex(i)%opttype==2) then
				vertex(i)%linkb%b1spin%sz=-vertex(i)%linkb%b1spin%sz
				vertex(i)%linkb%b2spin%sz=-vertex(i)%linkb%b2spin%sz
				if(evtx(i)%leg(1)%sz==-1) then
					evtx(i)%evtxtype=4
				else
					evtx(i)%evtxtype=5
				end if
			else
				if(evtx(i)%leg(1)%sz==evtx(i)%leg(2)%sz) then
					if(evtx(i)%leg(1)%sz==-1) then
						evtx(i)%evtxtype=1
					else 
						evtx(i)%evtxtype=6
					end if
				else 
					if(evtx(i)%leg(1)%sz==-1) then
						evtx(i)%evtxtype=2
					else
						evtx(i)%evtxtype=3
					end if
				end if
			end if
			evtx(i)%leg(3)=s(vertex(i)%linkb%b2spin%positioninsN)
			evtx(i)%leg(4)=s(vertex(i)%linkb%b1spin%positioninsN)
		
		end do		
		do i=1,N
			if(last(i)/=0) then
				link(last(i))=first(i)
				link(first(i))=last(i)
			end if
		end do	
		
!		do i=1, optN
!		write(*,*) "spin1(",vertex(i)%linkb%b1spin%position,")------connected with------spin2(",vertex(i)%linkb%b2spin%position,")"
!		write(*,*) "evtxleg 1,2,3,4",evtx(i)%leg(:)%sz,"evtxtype",evtx(i)%evtxtype
!		end do
!		do i=1, 4*optN
!			write(*,*) "leg",i,"is linked to",link(i)
!		end do

		tloopl=0	
		do loop=1,Nloop
!			write(*,*) "*********************************************************"
!			write(*,*) "loop",loop,"starts to be constructed"	
			leg0=4*optN*random(seed)+1
			l=leg0
			loopl=0		
100			i=(l+3)/4
			j=mod(l-1,4)+1
			r=random(seed)
			k=1
!			write(*,*) "########vertex",i,"beforeupdated,evtxtype",evtx(i)%evtxtype,"leg1,2,3,4",evtx(i)%leg(:)%sz
			do while(r>=hbprob(evtx(i)%evtxtype,j)%cumuprob(k))
				k=k+1
			end do
!			write(*,*) "   entranceleg",j,"               exitleg",hbprob(evtx(i)%evtxtype,j)%linkedleg(k)
			k2=hbprob(evtx(i)%evtxtype,j)%linkedleg(k)
			l=4*(i-1)+k2
			evtx(i)%leg(j)%sz=-evtx(i)%leg(j)%sz
			evtx(i)%leg(k2)%sz=-evtx(i)%leg(k2)%sz
			evtx(i)%evtxtype=new(evtx(i)%evtxtype,j,k2)
!			write(*,*) "                           afterupdated,evtxtype ",evtx(i)%evtxtype,"leg1,2,3,4",evtx(i)%leg(:)%sz	
			if(j/=k2) then 
				loopl=loopl+1
			end if
			if(l/=leg0.and.link(l)/=leg0) then
				l=link(l)
				goto 100
			end if
!			write(*,*) "loop",loop,"has been construced, its length=",loopl
!			write(*,*) "---------------------------------------------------------"
			tloopl=tloopl+loopl
!			write(*,*) "total length till now",tloopl
		end do 
		
		j=1
		do i=1,optM
			if(opt(i)%opttype/=0) then
				if(evtx(j)%evtxtype==4.or.evtx(j)%evtxtype==5) then
					opt(i)%opttype=2
				else 
					opt(i)%opttype=1
				end if
				vertex(j)%opttype=opt(i)%opttype
				j=j+1
			end if
		end do
		do i=1,N
			if(first(i)==0) then
				s(i)%sz=2*int(2.*random(seed))-1
			else
				j=first(i)
				s(i)%sz=evtx((j+3)/4)%leg(mod(j-1,4)+1)%sz
			end if	
		end do
	
		do i=1,N
			if(first(i)==0) then 
				cycle
			else
				j=first(i)
				k=last(i)
				if(evtx((j+3)/4)%leg(mod(j-1,4)+1)%sz/=evtx((k+3)/4)%leg(mod(k-1,4)+1)%sz) then
					write(*,*) "error in directedloopupdate"
				end if
			end if
		end do

		deallocate(vertex)
		deallocate(evtx)
		deallocate(link)
		deallocate(first)
		deallocate(last)	
		return
	end subroutine

	subroutine adjustcutoff()
		integer::i
		type(operator),allocatable::optcopy(:)	
		optMnew=optN+optN/4
		grow=0
		if(optMnew>optM) then 
			allocate(optcopy(optM))
			optcopy(:)=opt(:)
			deallocate(opt)
			allocate(opt(optMnew))
			opt(1:optM)=optcopy(:)
			do i=optM+1,optMnew
				opt(i)%opttype=0
				opt(i)%linkb=>b0
			end do
			deallocate(optcopy)
			optMold=optM
			optM=optMnew
			grow=1
		end if
!		write(*,*) "optM,optN,optMnew",optM,optN,optMnew
		return
	end subroutine

	subroutine deallocateall()
		deallocate(s)
		deallocate(b)
		deallocate(opt)
		deallocate(hbprob)
		return
	end subroutine

end module

program main
use pre
implicit none
	integer::i,i1,i2,tmp1,tmp2,nbalance,nbins,nsimulation,time(6)
	real*8::Mu,Ms,Muabs,Mu2,Mux,Msabs,Ms1,Ms2,Msx,optNbox,optN2box,Muabsbox,Msabsbox,Mu2box,Ms2box,Muxbox,Msxbox
	open(unit=30,file='adjustoff.txt',status='replace')
	open(unit=40,file='parameters.txt',status='replace')
	open(unit=50,file='result.txt',status='replace')
	nbalance=20000
	nbins=5
	nsimulation=20000
	call seedchoice()
	call makelattice()
	call parametersinitialize()
	call makeprobtable()
	call configinitialize()
	write(40,*) "===========+++++++++++================"
	write(40,*) "Dimension:             ",d
	if(d==1) then
		write(40,*) "Lattice size L:        ",Lx
	else if(d==2) then
		write(40,*) "Lattice size Lx*Ly:    ",Lx,Ly
	else
		write(40,*) "Lattice size Lx*Ly*Lz: ",Lx,Ly,Lz
	end if
	write(40,'(a,5f15.8)') " J,T,aniso,H,e:",J,T,aniso,H,e
	write(40,'(a,4i7)') " balancesteps,bin numbers,steps in each bin:",nbalance,nbins,nsimulation
	call itime(time(1:3))
	call idate(time(4:6))
	write(40,'(a,i2,a,i2,a,i2,a,i2,a,i2,a,i4)') " simulation start time:",&
						time(1),":",time(2),":",time(3)," ",time(4),"/",time(5),"/",time(6)		
	tmp1=0
	tmp2=0
	do i=1,nbalance
!		write(*,*) "++balancestep++",i,"optM=",optM,"optN=",optN
		Nloop=5
		call Dupdate()
		call directedloopupdate()
		if(i>4*nbalance/5) then
			tmp1=tmp1+tloopl
			tmp2=tmp2+optN	
		end if		
		call adjustcutoff()
		if(grow==1) then
			write(30,*) i,optMold,optMnew,optN 	
		end if
	end do
	Nloop=2*Nloop*tmp2/tmp1
!	write(*,*) "best Nloop",Nloop
	write(50,*) &
"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~results~~in~~each bin~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	write(50,*) &
"bins       E              C            <|uM|>	       <|sM|>	        Su             Ss             uX             sX"
	do i1=1,nbins
		optNbox=0
		optN2box=0
		Muabsbox=0
		Mu2box=0
		Muxbox=0
		Msabsbox=0
		Ms2box=0
		Msxbox=0
		do i2=1,nsimulation
			call Dupdate()
			call directedloopupdate()			
			Mu=0
			Ms=0
			do i=1,N
				Mu=Mu+s(i)%sz
				Ms=Ms+s(i)%sz*(-1)**(sum(s(i)%position))
			end do	
			Mu=Mu/2
			Muabs=abs(Mu)
			Mu2=Mu**2 
			Mux=Mu2
			Ms=Ms/2
			Msabs=0
			Ms1=0
			Ms2=0
			do i=1,optM
				if(opt(i)%opttype==0) then 
					cycle
				else if(opt(i)%opttype==2) then
					opt(i)%linkb%b1spin%sz=-opt(i)%linkb%b1spin%sz
					opt(i)%linkb%b2spin%sz=-opt(i)%linkb%b2spin%sz
					Ms=Ms+2*opt(i)%linkb%b1spin%sz*(-1)**(sum(opt(i)%linkb%b1spin%position))
				end if		
				Msabs=Msabs+abs(Ms)
				Ms1=Ms1+Ms
				Ms2=Ms2+Ms**2
			end do
			if(optN==0) then
				Msabs=abs(Ms)
				Ms2=Ms**2
				Msx=Ms2
			else 
				Msx=(Ms1**2+Ms2)/(optN*(optN+1))
				Msabs=Msabs/optN
				Ms2=Ms2/optN	
			end if
			optNbox=optNbox+optN
			optN2box=optN2box+optN**2
			Muabsbox=Muabsbox+Muabs
			Msabsbox=Msabsbox+Msabs
			Mu2box=Mu2box+Mu2
			Ms2box=Ms2box+Ms2
			Muxbox=Muxbox+Mux
			Msxbox=Msxbox+Msx			
		end do
		optNbox=optNbox/nsimulation
		optN2box=optN2box/nsimulation
		Muabsbox=Muabsbox/nsimulation
		Msabsbox=Msabsbox/nsimulation			
		Mu2box=Mu2box/nsimulation
		Ms2box=Ms2box/nsimulation
		Muxbox=Muxbox/nsimulation
		Msxbox=Msxbox/nsimulation
		energyav=-optNbox/(beta*N)+J*d*(aniso/4+hb+e)
		specificheatav=(optN2box-optNbox**2-optNbox)/N
		uniformedabsMav=Muabsbox/N
		staggeredabsMav=Msabsbox/N
		uniformedstructuref=Mu2box/N
		staggeredstructuref=Ms2box/N
		uniformedsusav=beta*Muxbox/N
		staggeredsusav=beta*Msxbox/N
		write(*,*) "E,C,<|uM|>,<|sM|>,Su,Ss,uX,sX",energyav,specificheatav,uniformedabsMav,staggeredabsMav,&
			uniformedstructuref,staggeredstructuref,uniformedsusav,staggeredsusav
		write(50,'(i3,8d15.4)') i1,energyav,specificheatav,uniformedabsMav,staggeredabsMav,uniformedstructuref,&
			staggeredstructuref,uniformedsusav,staggeredsusav
	end do
	call deallocateall()
	call itime(time(1:3))
        call idate(time(4:6))
        write(40,'(a,i2,a,i2,a,i2,a,i2,a,i2,a,i4)') " simulation end   time:",&
						time(1),":",time(2),":",time(3)," ",time(4),"/",time(5),"/",time(6)
	write(40,*) "===========-----------================"
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

