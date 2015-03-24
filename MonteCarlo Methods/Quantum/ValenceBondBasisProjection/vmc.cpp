#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
using namespace std;

struct Plist	//Plist为h(x,y),考虑格子对称性，x,y均为正值，其他的情况可以由对称变换得到
{
	int x,y;
	double p;
};

int sign(double x)	//自定义符号函数，自洽迭代时要用到
{
	return x>0?1:x<0?-1:0;	
}

void reform(int l[],int n)		//对l[n]序列的重组函数
{
    int i,j,k,m,cp[n];
	for(i=0;i<n;i++) cp[i]=l[i];
	for(i=0;i<n;i++)
	{
		j=rand()%(n-i);
		m=-1;
		for(k=0;k<n;k++)
		{
			if(cp[k]!=-1) m=m+1;
			if(m==j) break;
		}
		l[i]=cp[k];
		cp[k]=-1;		
	}	
}

int bisearch(double h[],double r,int n)	//在递增序列h中找到第一个元素大于r的位置
{
    int i,j,k;
    i=0;
    j=n-1;
    do
    {
        k=(i+j)/2;
        if(h[k]>r)  j=k;
        else    i=k;
    } while((i+1)!=j);
//    cout<<"h["<<i<<"]="<<h[i]<<"    h["<<j<<"]="<<h[j]<<"   r="<<r<<endl;
    return i+(h[i]<=r);
}

int trialpoint(int j,int x,int y,int L)	//bondloopupdate中，由当前点j,和根据二分法查找heatbath概率表得到的x,y，寻找试探点i
{
	int x0,y0,i;
	y0=j/L;
	if(y0%2)	x0=L-1-j%L;
	else	x0=j%L;
	x0=(x0+x+L)%L;
	y0=(y0+y+L)%L;
	if(y0%2)	i=y0*L+L-1-x0;
	else i=y0*L+x0;
	return i;
}

int effectivepos(int i,int j,int L,int *pmap)	//由i,j两点返回它们的有效相对位置在h[nPlist]中的标号
{
	int xi,yi,xj,yj,x,y;
	yi=i/L;
	if(yi%2)	xi=L-1-i%L;
	else	xi=i%L;
	yj=j/L;
	if(yj%2)	xj=L-1-j%L;
	else	xj=j%L;
	x=abs(xi-xj);
	y=abs(yi-yj);
	x=min(x,L-x);
	y=min(y,L-y);
	return	*(pmap+x*(L/2+1)+y);
}


int main(void)
{
	const int L=6,vmcsteps=10000,ministep=5,scsteps=10;
	int i,j,k,x,y,N,nPlist,i0,i1,j1,len,rchoose,length,istep,iscstep;
	N=L*L;
	int rfm[N/2],vl[N],vr[N],s[N],cpv[N],*v,map[L/2+1][L/2+1],*pmap,sig[N],sigl,nnb[N][4],nstat;
	double r,Ms,Msavtotal,E,Eavtotal,Msav,Eav;
	srand(static_cast<unsigned int>(time(0)));

//由trialpoint函数构建所有点的最近邻表nnb
	for(i=0;i<N;i++)
	{
		nnb[i][0]=trialpoint(i,1,0,L);
		nnb[i][1]=trialpoint(i,0,1,L);
		nnb[i][2]=trialpoint(i,-1,0,L);
		nnb[i][3]=trialpoint(i,0,-1,L);
	}

//找出独立的h(x,y)的数目nPlist,考虑正方格子情形
	nPlist=0;
	for(i=0;i<=L/2;i++)
	{
		for(j=i+1;j<=L/2;j=j+2)
		{
			nPlist=nPlist+1;
		}
	}
//map的作用是由知道相对位置为x,y的成键的两格点返回这个有效相对位置在h[nPlist]中的位置，先对其全部初始化为-1
	for(i=0;i<=L/2;i++)
	{
		for(j=0;j<=L/2;j++)	map[i][j]=-1;
	}

//下面创建概率表h(x,y)并初始化为所有成键概率相等，即Neel态
	double htable[nPlist],hnav[nPlist],hnavtotal[nPlist],nEav[nPlist],nEavtotal[nPlist],part[nPlist];	//htable为heatbath概率表，后面vmc中bond loop update要利用这个表
	int hn[nPlist];			//hn记录在一次统计的<vl|vr>环pattern中各种独立键h(x,y)出现的总次数
	Plist h[nPlist];
	nPlist=-1;	
	for(i=0;i<=L/2;i++)
	{
		for(j=i+1;j<=L/2;j=j+2)
		{
			nPlist=nPlist+1;
			h[nPlist].x=i;
			h[nPlist].y=j;
			h[nPlist].p=1;
			map[i][j]=nPlist;	//建立有效位置的map值
			map[j][i]=nPlist;
//			cout<<"h["<<nPlist<<"].x="<<i<<"	h["<<nPlist<<"].y="<<j<<endl;  测试哪些h(x,y)为独立的 
		}
	}
	nPlist++;

//检验map表，并测试effectivepos是否可以由两点i,j正确返回它们的有效相对位置标号	
/*	for(j=L/2;j>=0;j--)
	{
		for(i=0;i<=L/2;i++)	cout<<map[i][j]<<"	";
		cout<<endl;
	}	
*/
	pmap=map[0];	//pmap指向a[0][0]的地址，直接输入pmap=map会出现错误，除非map为一维数组
//	cout<<effectivepos(0,2,L,pmap)<<endl;

//下面建立VMC初始态，<Vl|Vr>,以及某个与之兼容的|Sz>组态
	for(i=0;i<N/2;i++)	rfm[i]=i;
	reform(rfm,N/2);
	for(i=0;i<N/2;i++)	
	{
		vl[2*i]=2*rfm[i]+1;
		vl[2*rfm[i]+1]=2*i;
	}
	reform(rfm,N/2);
	for(i=0;i<N/2;i++)
	{
		vr[2*i]=2*rfm[i]+1;
		vr[2*rfm[i]+1]=2*i;
	}
	for(i=0;i<N;i++)	s[i]=0;
	for(i=0;i<N;i=i+2)
	{
		if(s[i]==0)
		{	s[i]=2*(rand()%2)-1;
			j=vl[i];
			s[j]=-s[i];
			while(vr[j]!=i)
			{
				s[vr[j]]=-s[j];
				s[vl[vr[j]]]=s[j];
				j=vl[vr[j]];
			}	
//			cout<<"i="<<i<<endl; 此处可以修改以记录环的数目
		}
	}
//	for(i=0;i<N;i++) cout<<"vl["<<i<<"]="<<vl[i]<<"	vr["<<i<<"]="<<vr[i]<<"	s["<<i<<"]="<<s[i]<<endl; 	\
此处用来检验<vl|vr>是否正确形成环pattern


	for(iscstep=1;iscstep<=scsteps;iscstep++)
	{
//生成heatbath概率表
		htable[0]=h[0].p;
		for(i=1;i<nPlist;i++)
		{
			htable[i]=h[i].p+htable[i-1];	
		}

/*		for(i=0;i<nPlist;i++)	//对htable和h[].p进行归一化
		{
			h[i].p=h[i].p/htable[nPlist-1];
			htable[i]=htable[i]/htable[nPlist-1];
		}	
*/
		nstat=0;
		Msavtotal=0;
		Eavtotal=0;
		for(i=0;i<nPlist;i++)	
		{	
			hnavtotal[i]=0;
			nEavtotal[i]=0;
		}
		
		for(istep=0;istep<vmcsteps;istep++)
		{
/************下面为vmc中采用键环反转方式(bond loop updates)的一个Mcs START*******************/
			length=0;
			do
			{		
				if(rand()%2) v=vl;
				else v=vr;
				for(i=0;i<N;i++)	cpv[i]=*(v+i);	//cpv拷贝一份链接状态，在bondloopupdate过程中通过与之比较判断是否是bounce过程
				i0=rand()%N;	//i0为一个bondloopupdate过程的起始点
				i=i0;
				j=*(v+i);
				len=0;			//len为一个bondloopupdate过程中被改变的键数目
				do
				{
					do
					{
						r=(static_cast<double>(rand())/RAND_MAX)*htable[nPlist-1];	
		//				r=static_cast<double>(rand())/RAND_MAX;	如果对htable和h[].p归一化了就用这一行
						k=bisearch(htable,r,nPlist);				
						rchoose=rand()%8;
						switch(rchoose)
						{
							case 0: x=h[k].x;y=h[k].y;break;
							case 1: x=h[k].x;y=-h[k].y;break;
							case 2: x=-h[k].x;y=h[k].y;break;
							case 3: x=-h[k].x;y=-h[k].y;break;
							case 4: x=h[k].y;y=h[k].x;break;
							case 5: x=-h[k].y;y=h[k].x;break;
							case 6: x=h[k].y;y=-h[k].x;break;
							default : x=-h[k].y;y=-h[k].x;break;
						}	
						i1=trialpoint(j,x,y,L);
					} while(s[i1]==s[j]);
					if(i1==cpv[j]) len=--len+(i1==i0);
					else len++;		
					*(v+j)=i1;
					j1=*(v+i1);
					*(v+i1)=j;	
					i=i1;
					j=j1;
				} while(i!=i0);	
				length=length+len;
			} while(length<N);
//在进行了总长度为N左右的键环反转后，用产生的<vl|vr>环组态产生一个新的兼容|sz>态
			for(i=0;i<N;i++)    s[i]=0;
			for(i=0;i<N;i=i+2)
			{   
		        if(s[i]==0)
		        {   s[i]=2*(rand()%2)-1;
		            j=vl[i];
		            s[j]=-s[i];
		            while(vr[j]!=i)
		            {   
		                s[vr[j]]=-s[j];
	                	s[vl[vr[j]]]=s[j];
	            	    j=vl[vr[j]];
	        	    }   
	    	    }   
		    }   
/*********************************END Of 1MCS**************************************/
//			if(istep%ministep==0)	//是否想隔ministep个mcs抽一次样
//			{
				nstat++;
				for(i=0;i<nPlist;i++)	hn[i]=0;	//在统计之前先让所有独立h(x,y)的统计数目为0
				sigl=0;
				for(i=0;i<N;i++)	sig[i]=0;
				for(i=0;i<N;i=i+2)
				{
					if(sig[i]==0)	
					{
						sigl++;
						sig[i]=sigl;
						j=vl[i];
						sig[j]=sigl;
						hn[effectivepos(i,j,L,pmap)]++;
						while(vr[j]!=i)
						{
							hn[effectivepos(j,vr[j],L,pmap)]++;
							sig[vr[j]]=sigl;
							hn[effectivepos(vr[j],vl[vr[j]],L,pmap)]++;
							sig[vl[vr[j]]]=sigl;
							j=vl[vr[j]];
						}
						hn[effectivepos(j,i,L,pmap)]++;
					}	
				}
				Ms=0;
				E=0;			
				for(i=0;i<N;i++)
				{
					Ms=Ms+pow(-1,i%2)*s[i];
					E=E-(3.0/4)*((sig[i]==sig[nnb[i][0]])+(sig[i]==sig[nnb[i][1]]));
				}
				Msavtotal=Msavtotal+fabs(Ms);	
				Eavtotal=Eavtotal+E;
				for(i=0;i<nPlist;i++)
				{
					hnavtotal[i]=hnavtotal[i]+hn[i];
					nEavtotal[i]=nEavtotal[i]+hn[i]*E;
				}
//			}
		}
		Msav=Msavtotal/nstat;
		Eav=Eavtotal/nstat;
		for(i=0;i<nPlist;i++)
		{
			hnav[i]=hnavtotal[i]/(nstat*h[i].p);
			nEav[i]=nEavtotal[i]/(nstat*h[i].p);
		}
	
		cout<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(8);
		for(i=0;i<nPlist;i++)
		{
			part[i]=nEav[i]-hnav[i]*Eav;
//			cout<<"partial E to h"<<i<<"="<<setw(10)<<part[i]<<endl;
			h[i].p=h[i].p/exp((static_cast<double>(rand())/RAND_MAX)*sign(part[i])/pow(iscstep,0.75));
		}
		if(iscstep%5==0)	
		{
/*			cout<<iscstep<<" Msav="<<Msav/N<<" Eav="<<Eav/N<<endl;
			for(i=0;i<nPlist;i++)
			{
				cout<<"partial E to h"<<i<<"="<<setw(10)<<part[i]<<"	h["<<i<<"].p="<<h[i].p<<endl;
			}*/
			cout<<iscstep<<"	"<<Msav/N<<"	"<<Eav/N<<endl;
		}
	}
	for(i=0;i<nPlist;i++)
	{
		cout<<"partial E to h"<<i<<"="<<setw(10)<<part[i]<<"    h["<<i<<"].p="<<h[i].p<<endl;
	}
	return 0;
}
