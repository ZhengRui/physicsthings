#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
using namespace std;

struct segment
{
	segment *sdn,*sup,*ldn,*lup;
	struct vtx *vdn,*vup;
	double tdn,tup;
	int s,n,kkdn,kkup,c,ldnuod,lupuod;
};

struct vtx
{
	struct segment *p1,*p2,*p3,*p4;
	int type;
};

struct vtxlkcell
{
	struct vtx *pvtx;
	vtxlkcell *phead,*ptail;
};

struct loop
{
	int num,flip,tmp;
	double len,umag,smag;
};

struct rootloop
{
	int rooti;
	double len,umag,smag;
};

struct lpheadlkcell
{
	struct segment *lphead;
	lpheadlkcell *phead, *ptail;
};

struct kink
{
	int n1,n2;
	int s1,s2;
	double t;
	kink *khead, *ktail;
};

int xyz_to_n(int r[],int L,int d)	//由笛卡尔坐标xyz求出序号N
{
	int i,n,temp;
	n=0;
	temp=0;
	for(i=d-1;i>=0;i--)
	{
		n=n+(temp%2?(L-1-r[i]):r[i])*int(pow(L,i));
		temp=temp+r[i];
	}
	return n;
}

void formloop(struct segment *pstart, int uod, int c)
{
	struct segment *pnext;
	int uodnew;
	(*pstart).c=c;
//	cout<<"i:   "<<(*pstart).n<<"   headtime:"<<(*pstart).tdn<<"   tailtime:"<<(*pstart).tup<<endl;
	if(uod) pnext=(*pstart).lup, uodnew=1-(*pstart).lupuod; 
	else pnext=(*pstart).ldn, uodnew=1-(*pstart).ldnuod;
//	cout<<"i:   "<<(*pnext).n<<"   headtime:"<<(*pnext).tdn<<"   tailtime:"<<(*pnext).tup<<endl;
	if(!(*pnext).c) formloop(pnext, uodnew, c);
}

int main(void)
{
	const int L=100,d=1,N=int(pow(L,d));
	const long int mcsteps=50000,blsteps=10000,binsteps=10000;
	const double beta=10,Jxy=1,Jz=2;
	int i,j,k,bond[d*N][2],pos[N][d],nb[N][d][2],temp1,temp2,temp[d],s1,s2,lnum,lnum2,istep,nkink;
	double Jx,C,tdn,tup,tdntrial,dt,ltime,lum,lsm,Ediag,umag,smag,usus,ssus,data[8],avE,SHEAT,avUM,US,avSM,SS,smtemp;
	segment *wlhead[N],*p1,*p2,*p3,*p4,*ip1,*ip2,*p;
	vtxlkcell *vtx1,*vtx2,*vtx1_tail,*vtx2_tail,*vtemp;
	loop *lp;
	rootloop *lp2;
	lpheadlkcell *lpstart,*lpstart_tail,*ltemp;
	kink *kk1,*kk2,*ktemp,*kloc;
	time_t t;
	srand(static_cast<unsigned int>(time(0))); //srand设置随机数种子,rand()调用及RAND_MAX需要stdlib.h文件
	cout<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(8); //设置输出格式需要iomanip头文件
	ofstream ofile("exe.log",ios::out); //需要fstream头文件
	if(!ofile)
	{
		cerr<<"open error"<<endl;
		exit(1);
	}
	ofile<<"*************+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*****************"<<endl<<endl;
	ofile<<"L="<<L<<",   d="<<d<<",   beta="<<beta<<",   Jxy="<<Jxy<<",   Jz="<<Jz<<endl;	
	ofile<<"Total MCS="<<mcsteps<<",   Thermalization MCS="<<blsteps<<",   MCS in each bin="<<binsteps<<endl;
	time(&t);
	ofile<<endl<<"Start time: "<<ctime(&t)<<endl;
	ofile<<setiosflags(ios::scientific)<<setiosflags(ios::right)<<setprecision(5);

	for(i=0;i<N;i++)
	{
		wlhead[i]=new segment;
		(*wlhead[i]).sdn=wlhead[i];		//优先级 . 高于 * 高于 %
		(*wlhead[i]).sup=wlhead[i];
		(*wlhead[i]).tdn=0;
		(*wlhead[i]).tup=beta;
		(*wlhead[i]).s=2*(rand()%2)-1;
		(*wlhead[i]).n=i;
		(*wlhead[i]).c=0;
		(*wlhead[i]).kkdn=-1;
		(*wlhead[i]).kkup=-1;
		(*wlhead[i]).ldn=wlhead[i];
		(*wlhead[i]).lup=wlhead[i];
		(*wlhead[i]).vdn=NULL;	
		(*wlhead[i]).vup=NULL;
		(*wlhead[i]).lupuod=0;
		(*wlhead[i]).ldnuod=1;	
	}

//将第i个格点的笛卡尔坐标放到pos[i][0~d-1]数组中		
	for(i=0;i<N;i++)
	{		
		temp1=i/int(pow(L,d));
		temp2=i%int(pow(L,d));
		for(j=0;j<d;j++)
		{
			pos[i][d-1-j]=(temp1%2?int(pow(L,d-j))-temp2-1:temp2)/int(pow(L,d-1-j));
			temp2=(temp1%2?int(pow(L,d-j))-temp2-1:temp2)%int(pow(L,d-1-j));
			temp1=pos[i][d-1-j];
		}
//测试是否正确建立格点序号与格点坐标关系
/*
		cout<<"point "<<i<<" 's coordinates: ";
		for(j=0;j<d;j++)
		{
			cout<<pos[i][j]<<"	";
		}
		cout<<endl;
*/
//建立每个格点的最近邻nb[N][d][0-1]以及每个键所连接的两个格点的序号bond[N*d][0-1]
		for(j=0;j<d;j++)
		{
			for(k=0;k<d;k++) temp[k]=pos[i][k];
			temp[j]=(temp[j]+1)%L;
			nb[i][j][0]=xyz_to_n(temp,L,d);
			nb[xyz_to_n(temp,L,d)][j][1]=i;
			bond[i*d+j][0]=i;
			bond[i*d+j][1]=xyz_to_n(temp,L,d);
		}
	}	
//测试是否正确建立最近邻以及每个键所连接的格点序号
/*	for(i=0;i<N;i++)
	{
		cout<<"point "<<i<<" 's nearest neighbour: ";
		for(j=0;j<d;j++) cout<<nb[i][j][0]<<" "<<nb[i][j][1]<<" ";
		cout<<endl;
	}
	for(i=0;i<N*d;i++)
	{
		cout<<"bond "<<i<<" connect point "<<bond[i][0]<<" and point "<<bond[i][1]<<endl;
	}
*/	
	Jx=abs(Jxy);
	C=(Jz>=Jx?(-Jz)/4:(Jz>-Jx?(-Jx)/4:Jz/4));
	vtx1=new vtxlkcell, (*vtx1).phead=vtx1, (*vtx1).ptail=vtx1, (*vtx1).pvtx=NULL, vtx1_tail=vtx1;
	vtx2=new vtxlkcell, (*vtx2).phead=vtx2, (*vtx2).ptail=vtx2, (*vtx2).pvtx=NULL, vtx2_tail=vtx2;
	lpstart=new lpheadlkcell, (*lpstart).phead=lpstart, (*lpstart).ptail=lpstart, (*lpstart).lphead=NULL, lpstart_tail=lpstart;
	kk1=new kink, kk2=new kink;
	(*kk1).n1=-1, (*kk1).n2=-1, (*kk1).s1=0, (*kk1).s2=0, (*kk1).t=0, (*kk1).khead=kk1, (*kk1).ktail=kk2;
	(*kk2).n1=-1, (*kk2).n2=-1, (*kk2).s1=0, (*kk2).s2=0, (*kk2).t=beta, (*kk2).khead=kk1, (*kk2).ktail=kk2;
	for(i=0;i<8;i++) data[i]=0;
	cout<<"Bin      Energy            Spec_Heat         Uni_Mag        Uni_Sus        Stag_Mag       Stag_Sus"<<endl;
	ofile<<"Bin       Energy        Spec_Heat        Uni_Mag         Uni_Sus        Stag_Mag        Stag_Sus"<<endl;

	for(istep=1;istep<=mcsteps;istep++)
	{	
		if(Jz>=Jx)
		{
//			cout<<"bond   "<<"s1   "<<"s2   "<<"segment_start_time   "<<"segment_end_time   "<<"end_break_type"<<endl;
			for(i=0;i<d*N;i++)
			{
				p1=wlhead[bond[i][0]];
				p2=wlhead[bond[i][1]];
				tdn=0;
				tup=min((*p1).tup,(*p2).tup);
				p3=(*p1).sup;
				p4=(*p2).sup;
				do
				{
					if((*p1).s!=(*p2).s) tdntrial=tup;
					else tdntrial=tdn-2*log(static_cast<double>(rand())/RAND_MAX)/Jx;
//					cout<<i<<"      "<<(*p1).s<<"    "<<(*p2).s<<"          ";
					if(tdntrial>=tup) 
					{
//						cout<<tdn<<"        "<<tup<<"        ";
						tdn=tup;
						if(tup!=beta)
						{
							if((*p1).kkup==bond[i][1]&&(*p2).kkup==bond[i][0]) 
							{
//								cout<<"kink===>1"<<endl;
								(*p1).lup=p4, (*p1).lupuod=0;
								(*p4).ldn=p1, (*p4).ldnuod=1;
								(*p2).lup=p3, (*p2).lupuod=0;
								(*p3).ldn=p2, (*p3).ldnuod=1;
								(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*p3).vdn=(*p1).vup, (*p4).vdn=(*p1).vup;
								(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=p3, (*(*p1).vup).p4=p4, (*(*p1).vup).type=1;
						
								vtemp=new vtxlkcell;
								(*vtemp).phead=vtx1_tail, (*vtemp).ptail=vtemp, (*vtemp).pvtx=(*p1).vup;
								(*vtx1_tail).ptail=vtemp, vtx1_tail=vtemp;
						
								p1=p3;
								p2=p4;
								tup=min((*p1).tup,(*p2).tup);
								p3=(*p1).sup;
								p4=(*p2).sup;
							}
							else
							{
//								cout<<"trivial"<<endl;
								if((*p1).tup==tup) p1=p3, p3=(*p1).sup;
								if((*p2).tup==tup) p2=p4, p4=(*p2).sup;
								tup=min((*p1).tup,(*p2).tup);
							}
						}
//						else cout<<"boundary"<<endl;
					}
					else
					{
//						cout<<tdn<<"        "<<tdntrial<<"        "<<"1"<<endl;
						tdn=tdntrial;
						ip1=new segment, ip2=new segment;
						(*ip1).c=0, (*ip2).c=0;
	
						(*p3).sdn=ip1, (*(*p1).lup).ldn=ip1;
						(*ip1).s=(*p1).s, (*ip1).n=(*p1).n, (*ip1).sup=p3, (*ip1).sdn=p1, (*ip1).tup=(*p1).tup, (*ip1).tdn=tdn;
						(*ip1).kkup=(*p1).kkup, (*ip1).kkdn=-1; 
						(*ip1).lup=(*p1).lup, (*ip1).lupuod=0;
						(*p1).tup=tdn, (*p1).sup=ip1, (*p1).kkup=-1;
	
						(*p4).sdn=ip2, (*(*p2).lup).ldn=ip2;
						(*ip2).s=(*p2).s, (*ip2).n=(*p2).n, (*ip2).sup=p4, (*ip2).sdn=p2, (*ip2).tup=(*p2).tup, (*ip2).tdn=tdn;
						(*ip2).kkup=(*p2).kkup, (*ip2).kkdn=-1;
						(*ip2).lup=(*p2).lup, (*ip2).lupuod=0;
						(*p2).tup=tdn, (*p2).sup=ip2, (*p2).kkup=-1;
						
						if((*p1).vup!=NULL) ((*(*p1).vup).p1==p1?(*(*p1).vup).p1:(*(*p1).vup).p2)=ip1; //!!!!!!!判断条件?a:b=c实际上是先将c赋值给b后再与前面构成三目条件表达式，而(判断条件?a:b)=c是先进行条件表达式运算，再进行赋值，即=优先级高于三目运算符，切记切记!!!!!!!
						if((*p2).vup!=NULL) ((*(*p2).vup).p2==p2?(*(*p2).vup).p2:(*(*p2).vup).p1)=ip2;
						(*ip1).vup=(*p1).vup, (*ip2).vup=(*p2).vup;
						(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*ip1).vdn=(*p1).vup, (*ip2).vdn=(*p1).vup;
						(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=ip1, (*(*p1).vup).p4=ip2, (*(*p1).vup).type=1;
					
						vtemp=new vtxlkcell;
						(*vtemp).phead=vtx1_tail, (*vtemp).ptail=vtemp, (*vtemp).pvtx=(*p1).vup;
	 					(*vtx1_tail).ptail=vtemp, vtx1_tail=vtemp;
					
						(*p1).lup=ip2, (*p2).lup=ip1;
						(*ip2).ldn=p1, (*ip2).ldnuod=1;
						(*ip1).ldn=p2, (*ip1).ldnuod=1;
					
						p1=ip1;
						p2=ip2;					
					}
				}while(tdn!=beta);
	
				if(Jz>Jx)
				{
					p1=wlhead[bond[i][0]];
					p2=wlhead[bond[i][1]];
					tdn=0;
					tup=min((*p1).tup,(*p2).tup);
					p3=(*p1).sup;
					p4=(*p2).sup;
					do
					{	
						if((*p1).s!=(*p2).s) tdntrial=tup;
						else tdntrial=tdn-2*log(static_cast<double>(rand())/RAND_MAX)/(Jz-Jx);
//						cout<<i<<"      "<<(*p1).s<<"    "<<(*p2).s<<"          ";
						if(tdntrial>=tup)
						{
//							cout<<tdn<<"        "<<tup<<"        "<<"types_except_3"<<endl;
							tdn=tup;
							if(tup!=beta)
							{
								if((*p1).tup==tup) p1=(*p1).sup, p3=(*p1).sup;
								if((*p2).tup==tup) p2=(*p2).sup, p4=(*p2).sup;
								tup=min((*p1).tup,(*p2).tup);
							}
						}
						else
						{
//							cout<<tdn<<"        "<<tdntrial<<"        "<<"3"<<endl;
							tdn=tdntrial;
							ip1=new segment;ip2=new segment;
							(*ip1).c=0, (*ip2).c=0;
							
							(*p3).sdn=ip1, (*(*p1).lup).ldn=ip1;
							(*ip1).s=(*p1).s, (*ip1).n=(*p1).n, (*ip1).sup=p3, (*ip1).sdn=p1, (*ip1).tup=(*p1).tup, (*ip1).tdn=tdn;
							(*ip1).kkup=(*p1).kkup, (*ip1).kkdn=-1;
							(*ip1).lup=(*p1).lup, (*ip1).lupuod=0;
							(*p1).tup=tdn, (*p1).sup=ip1, (*p1).kkup=-1;
	
							(*p4).sdn=ip2, (*(*p2).lup).ldn=ip2;
							(*ip2).s=(*p2).s, (*ip2).n=(*p2).n, (*ip2).sup=p4, (*ip2).sdn=p2, (*ip2).tup=(*p2).tup, (*ip2).tdn=tdn;
							(*ip2).kkup=(*p2).kkup, (*ip2).kkdn=-1;
							(*ip2).lup=(*p2).lup, (*ip2).lupuod=0;
							(*p2).tup=tdn, (*p2).sup=ip2, (*p2).kkup=-1;
	
							if((*p1).vup!=NULL) ((*(*p1).vup).p1==p1?(*(*p1).vup).p1:(*(*p1).vup).p2)=ip1;
							if((*p2).vup!=NULL) ((*(*p2).vup).p2==p2?(*(*p2).vup).p2:(*(*p2).vup).p1)=ip2;
							(*ip1).vup=(*p1).vup, (*ip2).vup=(*p2).vup;
							(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*ip1).vdn=(*p1).vup, (*ip2).vdn=(*p1).vup;
							(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=ip1, (*(*p1).vup).p4=ip2, (*(*p1).vup).type=3;
	
							vtemp=new vtxlkcell;
							(*vtemp).phead=vtx2_tail, (*vtemp).ptail=vtemp, (*vtemp).pvtx=(*p1).vup;
							(*vtx2_tail).ptail=vtemp, vtx2_tail=vtemp;
	
							(*p1).lup=ip2, (*p2).lup=ip1;
							(*ip2).ldn=p1, (*ip2).ldnuod=1;
							(*ip1).ldn=p2, (*ip1).ldnuod=1;
	
							p1=ip1;
							p2=ip2;
						}
					}while(tdn!=beta);
				}
			}
		}

		else if(Jz<Jx && Jz>-Jx)
		{
//			cout<<"bond   "<<"s1   "<<"s2   "<<"segment_start_time   "<<"segment_end_time   "<<"end_break_type"<<endl;
			for(i=0;i<d*N;i++)
			{
				p1=wlhead[bond[i][0]];
				p2=wlhead[bond[i][1]];
				tdn=0;
				tup=min((*p1).tup,(*p2).tup);
				p3=(*p1).sup;
				p4=(*p2).sup;
				do
				{
					if((*p1).s==(*p2).s) tdntrial=tdn-4*log(static_cast<double>(rand())/RAND_MAX)/(Jx+Jz);
					else tdntrial=tdn-4*log(static_cast<double>(rand())/RAND_MAX)/(Jx-Jz);
//					cout<<i<<"      "<<(*p1).s<<"    "<<(*p2).s<<"          ";
					if(tdntrial>=tup)
					{
//						cout<<tdn<<"        "<<tup<<"        ";
						tdn=tup;
						if(tup!=beta)
						{
							if((*p1).kkup==bond[i][1]&&(*p2).kkup==bond[i][0])
							{
								if((static_cast<double>(rand())/RAND_MAX)<=(Jx+Jz)/(2*Jx))
								{
//									cout<<"kink===>1"<<endl;
									(*p1).lup=p4, (*p1).lupuod=0;
									(*p4).ldn=p1, (*p4).ldnuod=1;
									(*p2).lup=p3, (*p2).lupuod=0;
									(*p3).ldn=p2, (*p3).ldnuod=1;
									(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*p3).vdn=(*p1).vup, (*p4).vdn=(*p1).vup;
									(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=p3, (*(*p1).vup).p4=p4, (*(*p1).vup).type=1;
								}
								else
								{
//									cout<<"kink===>2"<<endl;
									(*p1).lup=p2, (*p1).lupuod=1;
									(*p2).lup=p1, (*p2).lupuod=1;
									(*p3).ldn=p4, (*p3).ldnuod=0;
									(*p4).ldn=p3, (*p4).ldnuod=0;
									(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*p3).vdn=(*p1).vup, (*p4).vdn=(*p1).vup;
									(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=p3, (*(*p1).vup).p4=p4, (*(*p1).vup).type=2;
								}

								vtemp=new vtxlkcell;
								(*vtemp).phead=vtx1_tail, (*vtemp).ptail=vtemp, (*vtemp).pvtx=(*p1).vup;
								(*vtx1_tail).ptail=vtemp, vtx1_tail=vtemp;

								p1=p3;
								p2=p4;
								tup=min((*p1).tup,(*p2).tup);
								p3=(*p1).sup;
								p4=(*p2).sup;
							}
							else
							{
//								cout<<"trivial"<<endl;
								if((*p1).tup==tup) p1=p3, p3=(*p1).sup;
								if((*p2).tup==tup) p2=p4, p4=(*p2).sup;
								tup=min((*p1).tup,(*p2).tup);
							}
						}
//						else cout<<"boundary"<<endl;
					}
					else
					{
						tdn=tdntrial;
						ip1=new segment, ip2=new segment;
						(*ip1).c=0, (*ip2).c=0;

						((*p1).lupuod?(*(*p1).lup).lup:(*(*p1).lup).ldn)=ip1;
						((*p2).lupuod?(*(*p2).lup).lup:(*(*p2).lup).ldn)=ip2;
						(*ip1).lup=(*p1).lup;
						(*ip2).lup=(*p2).lup;

						(*p3).sdn=ip1;
						(*ip1).s=(*p1).s, (*ip1).n=(*p1).n, (*ip1).sup=p3, (*ip1).sdn=p1, (*ip1).tup=(*p1).tup, (*ip1).tdn=tdn;
						(*ip1).kkup=(*p1).kkup, (*ip1).kkdn=-1;
						(*ip1).lupuod=(*p1).lupuod;
						(*p1).tup=tdn, (*p1).sup=ip1, (*p1).kkup=-1;

						(*p4).sdn=ip2;
						(*ip2).s=(*p2).s, (*ip2).n=(*p2).n, (*ip2).sup=p4, (*ip2).sdn=p2, (*ip2).tup=(*p2).tup, (*ip2).tdn=tdn;
						(*ip2).kkup=(*p2).kkup, (*ip2).kkdn=-1;
						(*ip2).lupuod=(*p2).lupuod;
						(*p2).tup=tdn, (*p2).sup=ip2, (*p2).kkup=-1;
						
						if((*p1).vup!=NULL) ((*(*p1).vup).p1==p1?(*(*p1).vup).p1:(*(*p1).vup).p2)=ip1;
						if((*p2).vup!=NULL) ((*(*p2).vup).p2==p2?(*(*p2).vup).p2:(*(*p2).vup).p1)=ip2;
						(*ip1).vup=(*p1).vup, (*ip2).vup=(*p2).vup;
						(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*ip1).vdn=(*p1).vup, (*ip2).vdn=(*p1).vup;
						(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=ip1, (*(*p1).vup).p4=ip2;

						vtemp=new vtxlkcell;
						(*vtemp).phead=vtx1_tail, (*vtemp).ptail=vtemp, (*vtemp).pvtx=(*p1).vup;
						(*vtx1_tail).ptail=vtemp, vtx1_tail=vtemp;

						if((*p1).s==(*p2).s)
						{
//							cout<<tdn<<"        "<<tdntrial<<"        "<<"1"<<endl;
							(*p1).lupuod=0;
							(*p2).lupuod=0;
							(*(*p1).vup).type=1;
							(*p1).lup=ip2, (*p2).lup=ip1;
							(*ip2).ldn=p1, (*ip2).ldnuod=1;
							(*ip1).ldn=p2, (*ip1).ldnuod=1;
						}
						else
						{
//							cout<<tdn<<"        "<<tdntrial<<"        "<<"2"<<endl;
							(*p1).lupuod=1;
							(*p2).lupuod=1;
							(*(*p1).vup).type=2;
							(*p1).lup=p2, (*p2).lup=p1;
							(*ip1).ldn=ip2, (*ip1).ldnuod=0;
							(*ip2).ldn=ip1, (*ip2).ldnuod=0;
						}

						p1=ip1;
						p2=ip2;
					}
				}while(tdn!=beta);
			}	
		}

		else
		{
//			cout<<"bond   "<<"s1   "<<"s2   "<<"segment_start_time   "<<"segment_end_time   "<<"end_break_type"<<endl;
			for(i=0;i<d*N;i++)
			{
				p1=wlhead[bond[i][0]];
				p2=wlhead[bond[i][1]];
				tdn=0;
				tup=min((*p1).tup,(*p2).tup);
				p3=(*p1).sup;
				p4=(*p2).sup;
				do
				{
					if((*p1).s==(*p2).s) tdntrial=tup;
					else tdntrial=tdn-2*log(static_cast<double>(rand())/RAND_MAX)/Jx;
//					cout<<i<<"      "<<(*p1).s<<"    "<<(*p2).s<<"          ";
					if(tdntrial>=tup)
					{
//						cout<<tdn<<"        "<<tup<<"        ";
						tdn=tup;
						if(tup!=beta)
						{
							if((*p1).kkup==bond[i][1]&&(*p2).kkup==bond[i][0])
							{
//								cout<<"kink===>2"<<endl;
								(*p1).lup=p2, (*p1).lupuod=1;
								(*p2).lup=p1, (*p2).lupuod=1;
								(*p3).ldn=p4, (*p3).ldnuod=0;
								(*p4).ldn=p3, (*p4).ldnuod=0;
								(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*p3).vdn=(*p1).vup, (*p4).vdn=(*p1).vup;
								(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=p3, (*(*p1).vup).p4=p4, (*(*p1).vup).type=2;

								vtemp=new vtxlkcell;
								(*vtemp).phead=vtx1_tail, (*vtemp).ptail=vtemp, (*vtemp).pvtx=(*p1).vup;
								(*vtx1_tail).ptail=vtemp, vtx1_tail=vtemp;

								p1=p3;
								p2=p4;
								tup=min((*p1).tup,(*p2).tup);
								p3=(*p1).sup;
								p4=(*p2).sup;
							}
							else
							{
//								cout<<"trivial"<<endl;
								if((*p1).tup==tup) p1=p3, p3=(*p1).sup;
								if((*p2).tup==tup) p2=p4, p4=(*p2).sup;
								tup=min((*p1).tup,(*p2).tup);
							}
						}
//						else cout<<"boundary"<<endl;
					}
					else
					{
//						cout<<tdn<<"        "<<tdntrial<<"        "<<"2"<<endl;
						tdn=tdntrial;
						ip1=new segment, ip2=new segment;
						(*ip1).c=0, (*ip2).c=0;

						((*p1).lupuod?(*(*p1).lup).lup:(*(*p1).lup).ldn)=ip1;
						((*p2).lupuod?(*(*p2).lup).lup:(*(*p2).lup).ldn)=ip2;
						(*ip1).lup=(*p1).lup;
						(*ip2).lup=(*p2).lup;

						//当L=2的时候，这个唯一的bond会考虑两次，所以就算是放2算符，也会出现和放4算符一样，在上面已经有跨越这个bond的算符的情况，需要特许考虑，而当L>2时，放2算符的时候，上方是不会有跨越这个bond的算符的。L=2时，注销掉上面一段，用下面一段来判断是否当前bond上面是否已经有算符。

/*						if((*p1).lup==p2)
						{
							(*ip1).lup=ip2;
							(*ip2).lup=ip1;
						}
						else
						{
							((*p1).lupuod?(*(*p1).lup).lup:(*(*p1).lup).ldn)=ip1;
							((*p2).lupuod?(*(*p2).lup).lup:(*(*p2).lup).ldn)=ip2;
							(*ip1).lup=(*p1).lup;
							(*ip2).lup=(*p2).lup;
						}
*/

						(*p3).sdn=ip1;
						(*ip1).s=(*p1).s, (*ip1).n=(*p1).n, (*ip1).sup=p3, (*ip1).sdn=p1, (*ip1).tup=(*p1).tup, (*ip1).tdn=tdn;
						(*ip1).kkup=(*p1).kkup, (*ip1).kkdn=-1;
					   	(*ip1).lupuod=(*p1).lupuod;
						(*p1).lupuod=1;
						(*p1).tup=tdn, (*p1).sup=ip1, (*p1).kkup=-1;

						(*p4).sdn=ip2;
						(*ip2).s=(*p2).s, (*ip2).n=(*p2).n, (*ip2).sup=p4, (*ip2).sdn=p2, (*ip2).tup=(*p2).tup, (*ip2).tdn=tdn;
						(*ip2).kkup=(*p2).kkup, (*ip2).kkdn=-1;
						(*ip2).lupuod=(*p2).lupuod;
						(*p2).lupuod=1;
						(*p2).tup=tdn, (*p2).sup=ip2, (*p2).kkup=-1;

						if((*p1).vup!=NULL) ((*(*p1).vup).p1==p1?(*(*p1).vup).p1:(*(*p1).vup).p2)=ip1;
						if((*p2).vup!=NULL) ((*(*p2).vup).p2==p2?(*(*p2).vup).p2:(*(*p2).vup).p1)=ip2;
						(*ip1).vup=(*p1).vup, (*ip2).vup=(*p2).vup;
						(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*ip1).vdn=(*p1).vup, (*ip2).vdn=(*p1).vup;
						(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=ip1, (*(*p1).vup).p4=ip2, (*(*p1).vup).type=2;

						vtemp=new vtxlkcell;
						(*vtemp).phead=vtx1_tail, (*vtemp).ptail=vtemp, (*vtemp).pvtx=(*p1).vup;
						(*vtx1_tail).ptail=vtemp, vtx1_tail=vtemp;

						(*p1).lup=p2, (*p2).lup=p1;
						(*ip1).ldn=ip2, (*ip1).ldnuod=0;
						(*ip2).ldn=ip1, (*ip2).ldnuod=0;

						p1=ip1;
						p2=ip2;
//						cout<<(*(*p1).lup).tdn<<"   "<<(*(*p1).lup).tup<<"   "<<(*(*p1).lup).n<<endl;
//						cout<<(*(*p2).lup).tdn<<"   "<<(*(*p2).lup).tup<<"   "<<(*(*p2).lup).n<<endl;
					}
				}while(tdn!=beta);
				
				if(Jz<-Jx)
				{
					p1=wlhead[bond[i][0]];
					p2=wlhead[bond[i][1]];
					tdn=0;
					tup=min((*p1).tup,(*p2).tup);
					p3=(*p1).sup;
					p4=(*p2).sup;
					do
					{
						if((*p1).s==(*p2).s) tdntrial=tup;
						else tdntrial=tdn-2*log(static_cast<double>(rand())/RAND_MAX)/(-Jz-Jx);
//						cout<<i<<"      "<<(*p1).s<<"    "<<(*p2).s<<"          ";
						if(tdntrial>=tup)
						{
//							cout<<tdn<<"        "<<tup<<"        "<<"types_except_4"<<endl;
							tdn=tup;
							if(tup!=beta)
							{
								if((*p1).tup==tup) p1=(*p1).sup, p3=(*p1).sup;
								if((*p2).tup==tup) p2=(*p2).sup, p4=(*p2).sup;
								tup=min((*p1).tup,(*p2).tup);
							}
						}
						else
						{
//								cout<<tdn<<"        "<<tdntrial<<"        "<<"4"<<endl;
								tdn=tdntrial;
								ip1=new segment;ip2=new segment;
								(*ip1).c=0, (*ip2).c=0;
								
								if((*p1).lup==p2)
								{
									(*ip1).lup=ip2;
									(*ip2).lup=ip1;
								}
								else
								{
									((*p1).lupuod?(*(*p1).lup).lup:(*(*p1).lup).ldn)=ip1;
									((*p2).lupuod?(*(*p2).lup).lup:(*(*p2).lup).ldn)=ip2;
									(*ip1).lup=(*p1).lup;
									(*ip2).lup=(*p2).lup;
								}
								
								(*p3).sdn=ip1;
								(*ip1).s=(*p1).s, (*ip1).n=(*p1).n, (*ip1).sup=p3, (*ip1).sdn=p1, (*ip1).tup=(*p1).tup, (*ip1).tdn=tdn;
								(*ip1).kkup=(*p1).kkup, (*ip1).kkdn=-1;
								(*ip1).lupuod=(*p1).lupuod;
								(*p1).lupuod=1;
								(*p1).tup=tdn, (*p1).sup=ip1, (*p1).kkup=-1;

								(*p4).sdn=ip2;
								(*ip2).s=(*p2).s, (*ip2).n=(*p2).n, (*ip2).sup=p4, (*ip2).sdn=p2, (*ip2).tup=(*p2).tup, (*ip2).tdn=tdn;
								(*ip2).kkup=(*p2).kkup, (*ip2).kkdn=-1;
								(*ip2).lupuod=(*p2).lupuod;
								(*p2).lupuod=1;
								(*p2).tup=tdn, (*p2).sup=ip2, (*p2).kkup=-1;

								if((*p1).vup!=NULL) ((*(*p1).vup).p1==p1?(*(*p1).vup).p1:(*(*p1).vup).p2)=ip1;
								if((*p2).vup!=NULL) ((*(*p2).vup).p2==p2?(*(*p2).vup).p2:(*(*p2).vup).p1)=ip2;
								(*ip1).vup=(*p1).vup, (*ip2).vup=(*p2).vup;
								(*p1).vup=new vtx, (*p2).vup=(*p1).vup, (*ip1).vdn=(*p1).vup, (*ip2).vdn=(*p1).vup;
								(*(*p1).vup).p1=p1, (*(*p1).vup).p2=p2, (*(*p1).vup).p3=ip1, (*(*p1).vup).p4=ip2, (*(*p1).vup).type=4;

								vtemp=new vtxlkcell;
								(*vtemp).phead=vtx2_tail, (*vtemp).ptail=vtemp, (*vtemp).pvtx=(*p1).vup;
								(*vtx2_tail).ptail=vtemp, vtx2_tail=vtemp;

								(*p1).lup=p2, (*p2).lup=p1;
								(*ip1).ldn=ip2, (*ip1).ldnuod=0;
								(*ip2).ldn=ip1, (*ip2).ldnuod=0;

								p1=ip1;
								p2=ip2;

//								cout<<(*(*p1).lup).tdn<<"   "<<(*(*p1).lup).tup<<"   "<<(*(*p1).lup).n<<endl;
//								cout<<(*(*p2).lup).tdn<<"   "<<(*(*p2).lup).tup<<"   "<<(*(*p2).lup).n<<endl;
						}
					}while(tdn!=beta);
				}		
			}
		}
/*		
		if((*vtx1).ptail!=vtx1)
		{		
			cout<<"vortices_type   "<<"n1   "<<"n2   "<<"time"<<endl;
			vtemp=vtx1;
			do
			{
				vtemp=(*vtemp).ptail;
				cout<<(*(*vtemp).pvtx).type<<"               "<<(*(*(*vtemp).pvtx).p1).n<<"    "<<\
(*(*(*vtemp).pvtx).p2).n<<"    "<<(*(*(*vtemp).pvtx).p1).tup<<endl;
			}while((*vtemp).ptail!=vtemp);
		}
	
		if((*vtx2).ptail!=vtx2)
		{
			cout<<"vortices_type   "<<"n1   "<<"n2   "<<"time"<<endl;
			vtemp=vtx2;
			do
			{
				vtemp=(*vtemp).ptail;
				cout<<(*(*vtemp).pvtx).type<<"               "<<(*(*(*vtemp).pvtx).p1).n<<"    "<<\
(*(*(*vtemp).pvtx).p2).n<<"    "<<(*(*(*vtemp).pvtx).p1).tup<<endl;
			}while((*vtemp).ptail!=vtemp);
		}
*/	
		lnum=1;
		for(i=0;i<N;i++)
		{
//			cout<<"worldline "<<i<<" 's segments' loop indices: "<<endl;
			p=wlhead[i];		
			do
			{
				if(!(*p).c)
				{
					ltemp=new lpheadlkcell;
					(*ltemp).phead=lpstart_tail, (*ltemp).ptail=ltemp, (*ltemp).lphead=p;
					(*lpstart_tail).ptail=ltemp,lpstart_tail=ltemp; 	
					formloop(p,1,lnum);
					lnum++;	
				}
//				cout<<(*p).c<<"   ";
				p=(*p).sup;	
			}while(p!=wlhead[i]);
//			cout<<endl;	
		}	
		
		lp=new loop[lnum];
		for(i=1;i<lnum;i++)
		{
			lp[i].num=i;
			lp[i].len=0;
			lp[i].umag=0;
			lp[i].smag=0;
			lp[i].tmp=i;
		}
		if((*vtx2).ptail!=vtx2)
		{
			vtemp=vtx2;
			do
			{
				vtemp=(*vtemp).ptail;
				i=(*(*(*vtemp).pvtx).p1).c;
				if((*(*vtemp).pvtx).type==3)	j=(*(*(*vtemp).pvtx).p2).c;
				else j=(*(*(*vtemp).pvtx).p3).c;
//				cout<<"loop1_index="<<i<<"   loop2_index="<<j<<endl;
				if(i!=j)
				{
					while(i!=lp[i].num) i=lp[i].num;
					while(j!=lp[j].num) j=lp[j].num;
					lp[max(i,j)].num=min(i,j);
				} 
			}while((*vtemp).ptail!=vtemp);		
		}
		lnum2=1;
		for(i=1;i<lnum;i++)
		{
			if(lp[i].num==i) 
			{
				lp[i].flip=2*(rand()%2)-1;
				lp[i].tmp=lnum2;
				lnum2++;
			}
			else 
			{
				lp[i].num=lp[lp[i].num].num;
				lp[i].flip=lp[lp[i].num].flip;
				lp[i].tmp=lp[lp[i].num].tmp;
			}
		}	
		lp2=new rootloop[lnum2];
		j=1;
		for(i=1;i<lnum;i++) 
		{
				if(lp[i].num==i) lp2[j].rooti=i, j++;
		}
		for(i=1;i<lnum2;i++) lp2[i].len=0,lp2[i].umag=0,lp2[i].smag=0;
	
		ltemp=lpstart;
		i=0;
		do
		{
			ltemp=(*ltemp).ptail;
			i++;
			p1=(*ltemp).lphead, p2=p1;
			ltime=0, lum=0, lsm=0, j=1;
			do
			{
				dt=(*p2).tup-(*p2).tdn;
				ltime=ltime+dt;
				lum=lum+dt*(*p2).s;
				if((*p2).n%2) lsm=lsm-dt*(*p2).s;
				else lsm=lsm+dt*(*p2).s;
				(*p2).s=lp[i].flip*(*p2).s;
				if(j) j=1-(*p2).lupuod, p2=(*p2).lup;
				else j=1-(*p2).ldnuod, p2=(*p2).ldn; 
			}while(p2!=p1);
			lp[i].len=ltime, lp[i].umag=lum, lp[i].smag=lsm;
			lp2[lp[i].tmp].len=lp2[lp[i].tmp].len+ltime;
			lp2[lp[i].tmp].umag=lp2[lp[i].tmp].umag+lum;
			lp2[lp[i].tmp].smag=lp2[lp[i].tmp].smag+lsm;			
		}while((*ltemp).ptail!=ltemp);	
	/*	
		cout<<"loop_indices   root_indices   flip    length            umag           smag"<<endl;
		for(i=1;i<lnum;i++)
		{
			cout<<i<<"               "<<lp[i].num<<"              "<<lp[i].flip<<"     "<<lp[i].len<<"     "<<lp[i].umag<<"     "<<lp[i].smag<<endl;
		}
		cout<<"root_indices   "<<"total_length     total_umag      total_smag"<<endl;
		for(i=1;i<lnum2;i++) cout<<lp2[i].rooti<<"              "<<lp2[i].len<<"     "<<lp2[i].umag<<"     "<<lp2[i].smag<<endl;
		for(i=0;i<N;i++)
		{
			cout<<"worldline "<<i<<" 's segments':"<<endl<<"indices   spinvalue   length"<<endl;
			p=wlhead[i];
			do
			{
				cout<<(*p).c<<"        "<<(*p).s<<"           "<<(*p).tup-(*p).tdn<<endl;
				p=(*p).sup;
			}while(p!=wlhead[i]);		
		}*/

		usus=0, ssus=0;
		for(i=1;i<lnum2;i++)
		{
			usus=usus+lp2[i].umag*lp2[i].umag;
			ssus=ssus+lp2[i].smag*lp2[i].smag;
		}
		usus=usus/4;
		ssus=ssus/4;

		nkink=0;
		if((*vtx1).ptail!=vtx1)
		{
			vtemp=(*vtx1).ptail;
			while(vtemp!=vtx1_tail)
			{
				if((*(*(*vtemp).pvtx).p1).s!=((*(*(*vtemp).pvtx).p3).s)) 
				{
					(*(*(*vtemp).pvtx).p1).kkup=(*(*(*vtemp).pvtx).p2).n;
					(*(*(*vtemp).pvtx).p2).kkup=(*(*(*vtemp).pvtx).p1).n;
					(*(*(*vtemp).pvtx).p3).kkdn=(*(*(*vtemp).pvtx).p4).n;
					(*(*(*vtemp).pvtx).p4).kkdn=(*(*(*vtemp).pvtx).p3).n;
					nkink++;
					kloc=kk1;
					while((*kloc).t<(*(*(*vtemp).pvtx).p1).tup) kloc=(*kloc).ktail;
					ktemp=new kink;
					(*ktemp).n1=(*(*(*vtemp).pvtx).p1).n;
					(*ktemp).s1=(*(*(*vtemp).pvtx).p1).s;
					(*ktemp).n2=(*(*(*vtemp).pvtx).p2).n;
					(*ktemp).s2=(*(*(*vtemp).pvtx).p2).s;
					(*ktemp).t=(*(*(*vtemp).pvtx).p1).tup;
					(*ktemp).khead=(*kloc).khead;
					(*ktemp).ktail=kloc;
					(*(*kloc).khead).ktail=ktemp;
					(*kloc).khead=ktemp;
				}
			/* 下面else的情况屏蔽掉也无所谓，因为在后面segments合并的过程中，会将这些中间段的信息抹掉	
				else 
				{
					(*(*(*vtemp).pvtx).p1).kkup=-1;
					(*(*(*vtemp).pvtx).p2).kkup=-1;
					(*(*(*vtemp).pvtx).p3).kkdn=-1;
					(*(*(*vtemp).pvtx).p4).kkdn=-1;
				}
			*/
				delete (*vtemp).pvtx;
				vtemp=(*vtemp).ptail;
				delete (*vtemp).phead;
			}
			if((*(*(*vtemp).pvtx).p1).s!=((*(*(*vtemp).pvtx).p3).s))
			{
				(*(*(*vtemp).pvtx).p1).kkup=(*(*(*vtemp).pvtx).p2).n;
				(*(*(*vtemp).pvtx).p2).kkup=(*(*(*vtemp).pvtx).p1).n;
				(*(*(*vtemp).pvtx).p3).kkdn=(*(*(*vtemp).pvtx).p4).n;
				(*(*(*vtemp).pvtx).p4).kkdn=(*(*(*vtemp).pvtx).p3).n;
				nkink++;
				kloc=kk1;
				while((*kloc).t<(*(*(*vtemp).pvtx).p1).tup) kloc=(*kloc).ktail;
				ktemp=new kink;
				(*ktemp).n1=(*(*(*vtemp).pvtx).p1).n;
				(*ktemp).s1=(*(*(*vtemp).pvtx).p1).s;
				(*ktemp).n2=(*(*(*vtemp).pvtx).p2).n;
				(*ktemp).s2=(*(*(*vtemp).pvtx).p2).s;
				(*ktemp).t=(*(*(*vtemp).pvtx).p1).tup;
				(*ktemp).khead=(*kloc).khead;
				(*ktemp).ktail=kloc;
				(*(*kloc).khead).ktail=ktemp;
				(*kloc).khead=ktemp;
			}
		/*
			else
			{
				(*(*(*vtemp).pvtx).p1).kkup=-1;
				(*(*(*vtemp).pvtx).p2).kkup=-1;
				(*(*(*vtemp).pvtx).p3).kkdn=-1;
				(*(*(*vtemp).pvtx).p4).kkdn=-1;
			}
		*/
			delete (*vtemp).pvtx;
			delete vtemp;
			(*vtx1).ptail=vtx1;
			vtx1_tail=vtx1;
		}
	
		if((*vtx2).ptail!=vtx2)
		{
			vtemp=(*vtx2).ptail;
			while(vtemp!=vtx2_tail)
			{
				delete (*vtemp).pvtx;
				vtemp=(*vtemp).ptail;
				delete (*vtemp).phead;
			}
			delete (*vtemp).pvtx;
			delete vtemp;
			(*vtx2).ptail=vtx2;
			vtx2_tail=vtx2;
		}
	
		ltemp=(*lpstart).ptail;
		while(ltemp!=lpstart_tail)
		{
			ltemp=(*ltemp).ptail;
			delete (*ltemp).phead;
		}
		delete ltemp;
		(*lpstart).ptail=lpstart;
		lpstart_tail=lpstart;
	
		delete [] lp;
		delete [] lp2;
	
		umag=0;
//		smag=0;	
		for(i=0;i<N;i++)
		{
			p1=wlhead[i], p2=(*p1).sup;
			umag=umag+(*p1).s;
//			if((*p1).n%2) smag=smag-(*p1).s;
//			else smag=smag+(*p1).s;
			j=1;
			while(j)
			{
				k=0;
				while((*p2).s==(*p1).s)
				{
					if(p2==wlhead[i])
					{
						j=0;
						break;
					}
					k=1;
					p2=(*p2).sup;
					delete (*p2).sdn;
				}
				if(k)
				{
					(*p1).sup=p2;
					(*p2).sdn=p1;
					(*p1).kkup=(*p2).kkdn;
					if(j) (*p1).tup=(*p2).tdn;
					else (*p1).tup=beta;
				}
				(*p1).lup=p2, (*p1).vup=NULL, (*p1).c=0;
				(*p2).ldn=p1, (*p2).vdn=NULL, (*p2).c=0;
				(*p1).lupuod=0, (*p2).ldnuod=1;
				p1=p2, p2=(*p1).sup;
			}
		}
		umag=abs(umag)/2;
//		smag=abs(smag)/2;
/*
		for(i=0;i<N;i++)
		{
			p1=wlhead[i], p2=(*p1).lup;
			while(p2!=wlhead[i])
			{
				cout<<(*p1).tdn<<"   "<<(*p1).tup<<"   "<<(*p1).lupuod<<"   "<<(*p1).ldnuod<<"   "<<(*p1).n<<endl;
				cout<<(*p2).tdn<<"   "<<(*p2).tup<<"   "<<(*p2).lupuod<<"   "<<(*p2).ldnuod<<"   "<<(*p2).n<<endl;
				cout<<endl;
				p1=p2;
				p2=(*p1).lup;
			}
		}
*/

		smag=0, smtemp=0;
		for(i=0;i<N;i++)
		{
			p=wlhead[i];
			if((*p).n%2) smtemp=smtemp-(*p).s;
			else smtemp=smtemp+(*p).s;
		}
		ktemp=(*kk1).ktail;
		tdn=0, tup=(*ktemp).t;
		while(ktemp!=kk2)
		{
			smag=smag+abs(smtemp)*(tup-tdn);
			if((*ktemp).n1%2) smtemp=smtemp+4*(*ktemp).s1;
			else smtemp=smtemp-4*(*ktemp).s1;
			ktemp=(*ktemp).ktail;
			tdn=tup, tup=(*ktemp).t;
			delete (*ktemp).khead;
		}
		smag=smag+abs(smtemp)*(tup-tdn);
		smag=smag/(2*beta);
		(*kk1).ktail=kk2, (*kk2).khead=kk1;

		Ediag=0;
		for(i=0;i<d*N;i++)
		{
			p1=wlhead[bond[i][0]];
			p2=wlhead[bond[i][1]];
			tdn=0;
			tup=min((*p1).tup,(*p2).tup);
			while(tup!=beta)
			{
				Ediag=Ediag+(tup-tdn)*(-Jz)*(*p1).s*(*p2).s;
				tdn=tup;
				if((*p1).tup==tup) p1=(*p1).sup;
				if((*p2).tup==tup) p2=(*p2).sup;
				tup=min((*p1).tup,(*p2).tup);
			}
			Ediag=Ediag+(tup-tdn)*(-Jz)*(*p1).s*(*p2).s;
		}
		Ediag=Ediag/4;		

		if(istep>blsteps)
		{
			data[0]=data[0]+Ediag;
			data[1]=data[1]+Ediag*Ediag;
			data[2]=data[2]+nkink;
			data[3]=data[3]+nkink*nkink;
			data[4]=data[4]+umag;
			data[5]=data[5]+usus;
			data[6]=data[6]+smag;
			data[7]=data[7]+ssus;
			if(!((istep-blsteps)%binsteps))
			{
				avE=(data[0]-data[2])/(binsteps*beta*N);
				SHEAT=(data[1]+data[3]-data[2])/(N*binsteps)-(data[0]*data[0]+data[2]*data[2])/(N*binsteps*binsteps);
				avUM=data[4]/(binsteps*N);
				US=data[5]/(binsteps*beta*N);
				avSM=data[6]/(binsteps*N);
				SS=data[7]/(binsteps*beta*N);
				cout<<(istep-blsteps)/binsteps<<"     "<<avE<<"     "<<SHEAT<<"     "<<avUM<<"     "<<US<<"     "<<avSM<<"     "<<SS<<endl;
				ofile<<(istep-blsteps)/binsteps<<"     "<<avE<<"     "<<SHEAT<<"     "<<avUM<<"     "<<US<<"     "<<avSM<<"     "<<SS<<endl;
				for(i=0;i<8;i++) data[i]=0;
			}
		}
	}
	time(&t);
	ofile<<endl<<"End   time: "<<ctime(&t)<<endl;
	ofile<<"*************-----------------------------------------------------------------------*****************"<<endl<<endl;
		
	return 0;
}
