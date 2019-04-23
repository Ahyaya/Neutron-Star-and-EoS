#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#define pi 3.14159265358979
#define Mscale 2.033931261665867e5
#define c 2.99792458e8

FILE *Xout;
double dM,dMa,dMp,dp,dy,dI,dAg00,temp_rho,temp_E,temp_cs2;
double temp_R,temp_M,temp_I,temp_Ma,temp_Mp,temp_k2,temp_Lambda;
double p_surf,E_core,EOS_lgDen[3000],EOS_lgE[3000],EOS_lgp[3000];
double Ecentral[3000];
int length_EOS=0,length_DenFile=0,loopTOV=1,neos=1;
char sT[30];
struct EoS
	{
		char path[60],filename[36],name[32];
	};
struct EoS list[30001];

/*pchip interpolation*/
int interp_p2rhoE(cp)
double cp;
{
if(cp<p_surf)
{loopTOV=0;return 0;}

double logp=log10(cp),logE,cE,a0,a1,b0,b1,h0,h1,h2,d0,d1,d2,k0,k1,dydx;
double logrho,drho0,drho1,drho2,krho0,krho1;
int n=length_EOS-1;
int np=n;
for(;n>0;n--)
{
	if(logp>EOS_lgp[n-1])
	{
		if((n>1)&&(n<np))
		{
			h0=EOS_lgp[n-1]-EOS_lgp[n-2];
			h1=EOS_lgp[n]-EOS_lgp[n-1];
			h2=EOS_lgp[n+1]-EOS_lgp[n];
			d0=(EOS_lgE[n-1]-EOS_lgE[n-2])/h0;
			drho0=(EOS_lgDen[n-1]-EOS_lgDen[n-2])/h0;
			d1=(EOS_lgE[n]-EOS_lgE[n-1])/h1;
			drho1=(EOS_lgDen[n]-EOS_lgDen[n-1])/h1;
			d2=(EOS_lgE[n+1]-EOS_lgE[n])/h2;
			drho2=(EOS_lgDen[n+1]-EOS_lgDen[n])/h2;
			k0=3*d0*d1/(d0+d1+(h0*d0+h1*d1)/(h0+h1));
			krho0=3*drho0*drho1/(drho0+drho1+(h0*drho0+h1*drho1)/(h0+h1));
			k1=3*d2*d1/(d2+d1+(h2*d2+h1*d1)/(h2+h1));
			krho1=3*drho2*drho1/(drho2+drho1+(h2*drho2+h1*drho1)/(h2+h1));
		}else if(n==1)
		{
			h1=EOS_lgp[n]-EOS_lgp[n-1];
			h2=EOS_lgp[n+1]-EOS_lgp[n];
			d1=(EOS_lgE[n]-EOS_lgE[n-1])/h1;
			drho1=(EOS_lgDen[n]-EOS_lgDen[n-1])/h1;
			d2=(EOS_lgE[n+1]-EOS_lgE[n])/h2;
			drho2=(EOS_lgDen[n+1]-EOS_lgDen[n])/h2;
			k0=((2*h1+h2)*d1-h1*d2)/(h1+h2);
			krho0=((2*h1+h2)*drho1-h1*drho2)/(h1+h2);
			if(k0*d1<=0)
			{
				k0=0;
			}else if((d1*d2<=0)&&(k0*k0>9*d1*d1))
			{
				k0=3*d1;
			}
			k1=3*d2*d1/(d2+d1+(h2*d2+h1*d1)/(h2+h1));
			if(krho0*drho1<=0)
			{
				krho0=0;
			}else if((drho1*drho2<=0)&&(krho0*krho0>9*drho1*drho1))
			{
				krho0=3*drho1;
			}
			krho1=3*drho2*drho1/(drho2+drho1+(h2*drho2+h1*drho1)/(h2+h1));
		}else if(n==np)
		{
			h0=EOS_lgp[n-1]-EOS_lgp[n-2];
			h1=EOS_lgp[n]-EOS_lgp[n-1];
			d0=(EOS_lgE[n-1]-EOS_lgE[n-2])/h0;
			drho0=(EOS_lgDen[n-1]-EOS_lgDen[n-2])/h0;
			d1=(EOS_lgE[n]-EOS_lgE[n-1])/h1;
			drho1=(EOS_lgDen[n]-EOS_lgDen[n-1])/h1;
			k0=3*d0*d1/(d0+d1+(h0*d0+h1*d1)/(h0+h1));
			krho0=3*drho0*drho1/(drho0+drho1+(h0*drho0+h1*drho1)/(h0+h1));
			k1=((2*h1+h0)*d1-h1*d0)/(h1+h0);
			krho1=((2*h1+h0)*drho1-h1*drho0)/(h1+h0);
			if(k1*d1<=0)
			{
				k1=0;
			}else if((d0*d1<=0)&&(k1*k1>9*d1*d1))
			{
				k1=3*d1;
			}
			if(krho1*drho1<=0)
			{
				krho1=0;
			}else if((drho0*drho1<=0)&&(krho1*krho1>9*drho1*drho1))
			{
				krho1=3*drho1;
			}
		}else
		{
			printf("Warning: pchip error 1\n");
			return 0;
		}
		a0=(1+2/h1*(logp-EOS_lgp[n-1]))*(logp-EOS_lgp[n])*(logp-EOS_lgp[n])/h1/h1;
		a1=(1-2/h1*(logp-EOS_lgp[n]))*(logp-EOS_lgp[n-1])*(logp-EOS_lgp[n-1])/h1/h1;
		b0=(logp-EOS_lgp[n-1])*(logp-EOS_lgp[n])*(logp-EOS_lgp[n])/h1/h1;
		b1=(logp-EOS_lgp[n])*(logp-EOS_lgp[n-1])*(logp-EOS_lgp[n-1])/h1/h1;
		logE=a0*EOS_lgE[n-1]+a1*EOS_lgE[n]+k0*b0+k1*b1;
		logrho=a0*EOS_lgDen[n-1]+a1*EOS_lgDen[n]+krho0*b0+krho1*b1;
		cE=pow(10,logE);
		temp_rho=pow(10,logrho);
		temp_E=cE;
		dydx=6/h1/h1/h1*(logp-EOS_lgp[n])*(logp-EOS_lgp[n-1])*(EOS_lgE[n-1]-EOS_lgE[n])+k0/h1/h1*(logp-EOS_lgp[n])*(3*logp-EOS_lgp[n]-EOS_lgp[n-1]-EOS_lgp[n-1])+k1/h1/h1*(logp-EOS_lgp[n-1])*(3*logp-EOS_lgp[n-1]-EOS_lgp[n]-EOS_lgp[n]);
		temp_cs2=cp/cE/dydx;
		return 1;
	}
}
printf("Warning: pchip error 2\n");
return 0;
}


/*interpolation from EnergyDensity to Pressure*/
double interp_E2p(cE)
double cE;
{
double logp,logE=log10(cE),cp;
int n=length_EOS-1;
for(;n>0;n--)
{
if(logE>EOS_lgE[n-1])
{
logp=(logE-EOS_lgE[n-1])*(EOS_lgp[n]-EOS_lgp[n-1])/(EOS_lgE[n]-EOS_lgE[n-1])+EOS_lgp[n-1];
cp=pow(10,logp);
return(cp);
}
}
printf("Warning:E2p interpolation error\n");
return(0);
}


/*Euler step*/
int fEul(fr,fm,fp,fy,fAg00,h)
double fr,fm,fp,fy,fAg00,h;
{
	double fpir=4.0*pi*fr*fr,hbt=1.0/sqrt(1-2*fm/fr),fE,frho;
	interp_p2rhoE(fp);
	fE=temp_E;
	frho=temp_rho;
	dM=fpir*fE*h;
	dMa=fpir*h*frho*hbt;
	dMp=fpir*h*fE*hbt;
	dp=h*(fE+fp)*(fm+4*pi*fr*fr*fr*fp)/(2*fm*fr-fr*fr);
	dI=8.0/3*h*pi*fr*fr*fr*fr*(fE+fp)*hbt/sqrt(fAg00);
	dAg00=-2.0*fAg00*dp/(fE+fp);
	dy=h*(-fy*fy/fr-fy/fr*(1+fpir*(fp-fE))/(1-2*fm/fr)-fr*(4*pi*(5*fE+9*fp+(fp+fE)/temp_cs2)/(1-2*fm/fr)-6/fr/fr/(1-2*fm/fr)-4/fr/fr/fr/fr/(1-2*fm/fr)/(1-2*fm/fr)*(fm+4*pi*fr*fr*fr*fp)*(fm+4*pi*fr*fr*fr*fp)));
	return 1;
}


/*Main iteration process*/
int getRM(Ec)
double Ec;
{
	loopTOV=1;
	double fEc=Ec*6.6741e-11;
	double h=2e-8,localRes,r,y;
	double m,ma,mp,p=interp_E2p(fEc),E,rho,I,Ag00,Bg00;
	double m11,p11,y11,m21,p21,y21,m22,p22,y22,m31,p31,y31,m61,p61,y61,m62,p62,y62;
	double ma21,ma22,ma31,ma61,ma62,mp21,mp22,mp31,mp61,mp62;
	double Ag0021,Ag0022,Ag0031,Ag0061,Ag0062,I21,I22,I31,I61,I62;
	double TR3m,TR3ma,TR3mp,TR3p,TR3I,TR3Ag00,TR3y,TR2m,TR2p,TR2y;
	double res1,res2,res3;
	int repeat=0;
	r=h;
	interp_p2rhoE(p);
	E=temp_E;
	rho=temp_rho;
	m=4.0/3*pi*r*r*r*E;
	ma=4.0/3*pi*r*r*r*rho;
	mp=m;
	y=2.0;
	I=0.0;
	Ag00=1.0;
	for(repeat=0;repeat<5e4;repeat++)
	{
		fEul(r,m,p,y,Ag00,h);
		m11=m+dM;
		p11=p+dp;
		y11=y+dy;
		
		m21=m+0.5*dM;
		p21=p+0.5*dp;
		y21=y+0.5*dy;
		ma21=ma+0.5*dMa;
		mp21=mp+0.5*dMp;
		I21=I+0.5*dI;
		Ag0021=Ag00+0.5*dAg00;
		
		m31=m+dM/3;
		p31=p+dp/3;
		y31=y+dy/3;
		ma31=ma+dMa/3;
		mp31=mp+dMp/3;
		I31=I+dI/3;
		Ag0031=Ag00+dAg00/3;
		
		m61=m+dM/6;
		p61=p+dp/6;
		y61=y+dy/6;
		ma61=ma+dMa/6;
		mp61=mp+dMp/6;
		I61=I+dI/6;
		Ag0061=Ag00+dAg00/6;
		
		fEul(r+0.5*h,m21,p21,y21,Ag0021,0.5*h);
		m22=m21+dM;
		p22=p21+dp;
		y22=y21+dy;
		ma22=ma21+dMa;
		mp22=mp21+dMp;
		I22=I21+dI;
		Ag0022=Ag0021+dAg00;
		
		fEul(r+h/6,m61,p61,y61,Ag0061,h/6);
		m62=m61+dM;
		p62=p61+dp;
		y62=y61+dy;
		ma62=ma61+dMa;
		mp62=mp61+dMp;
		I62=I61+dI;
		Ag0062=Ag0061+dAg00;
		
		if(loopTOV<1) break;
		
		TR3m=m22+9*(m62-m31);
		TR3p=p22+9*(p62-p31);
		TR3y=y22+9*(y62-y31);
		TR3ma=ma22+9*(ma62-ma31);
		TR3mp=mp22+9*(mp62-mp31);
		TR3I=I22+9*(I62-I31);
		TR3Ag00=Ag0022+9*(Ag0062-Ag0031);
		TR2m=2*m22-m11;
		TR2p=2*p22-p11;
		TR2y=2*y22-y11;
		
		res1=(TR3m-TR2m)/(TR3m-m);
		res1=res1>0?res1:-res1;
		res2=(TR3p-TR2p)/(TR3p-p);
		res2=res2>0?res2:-res2;
		res3=(TR3y-TR2y)/(TR3y-y);
		res3=res3>0?res3:-res3;
		res2=(res2>res3)?res2:res3;
		localRes=(res1>res2)?res1:res2;
		
		if(localRes<1e-6)
		{
			r=r+h;
			m=TR3m;
			p=TR3p;
			y=TR3y;
			ma=TR3ma;
			mp=TR3mp;
			I=TR3I;
			Ag00=TR3Ag00;
			h=h*2;
		}else if(localRes>2e-5 && h>1e-8){
			h=h*0.5;
		}else{
			r=r+h;
			m=TR3m;
			p=TR3p;
			y=TR3y;
			ma=TR3ma;
			mp=TR3mp;
			I=TR3I;
			Ag00=TR3Ag00;
		}
		
	}
	
	double bt=m/r;
	Bg00=(1-2.0*m/r)/Ag00;
	temp_I=I/(m*r*r*sqrt(Bg00));
	temp_R=r*c*1e-3;
	temp_M=m*Mscale;
	temp_k2=1.6*bt*bt*bt*bt*bt*(1-2*bt)*(1-2*bt)*(2-y+2*bt*(y-1))/(2*bt*(6-3*y+3*bt*(5*y-8))+4*bt*bt*bt*(13-11*y+bt*(3*y-2)+2*bt*bt*(1+y))+3*(1-2*bt)*(1-2*bt)*(2-y+2*bt*(y-1))*log(1-2*bt));
	temp_Ma=ma*Mscale;
	temp_Mp=mp*Mscale;
	temp_Lambda=9495*temp_k2*pow(temp_R/10,5)/pow(temp_M,5);
	return 1;
}


/*EoS files import*/
int loadEoS(pathname)
char pathname[60];
{
	FILE *inf;
	int n;
	if((inf=fopen(pathname,"r"))==NULL)
	{
		printf("\ncannot open %s\n",pathname);
		return -1;
	}
	length_EOS=0;
	while(fscanf(inf,"%lf",&EOS_lgE[length_EOS])==1)
	{
		fscanf(inf,"%lf",&EOS_lgp[length_EOS]);
		fscanf(inf,"%lf%*[^\n]",&EOS_lgDen[length_EOS]);
		length_EOS++;
	}
	fclose(inf);
	E_core=pow(10,EOS_lgE[length_EOS-1]);
	
	for(n=0;n<length_EOS;n++)
	{
		EOS_lgE[n]=EOS_lgE[n]-10.175607290470733;
		EOS_lgDen[n]=EOS_lgDen[n]-10.175607290470733;
		EOS_lgp[n]=EOS_lgp[n]-27.129849799910058;
	}
	p_surf=pow(10,EOS_lgp[0]);
	return 1;
}


/*load density file*/
int loadDenFile()
{
	FILE *DenFile;
	double Ec;
	length_DenFile=0;
	if((DenFile=fopen("./Central_Density/DenSeq.txt","r"))==NULL)
	{
		printf("Cannot load density file: ./Central_Density/DenSeq.txt");
		return 0;
	}
	while(fscanf(DenFile,"%lf",&Ec)==1)
	{
		Ecentral[length_DenFile]=Ec;
		length_DenFile++;
	}
	fclose(DenFile);
	return 1;
}


/*Working Mode: Density File Dependent*/
int DenMode()
{
	int pf;
	for(pf=0;pf<length_DenFile;pf++)
	{
		getRM(Ecentral[pf]);
		fprintf(Xout,"%e\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Ecentral[pf],temp_M,temp_R,temp_I,temp_Lambda,temp_k2,temp_Ma,temp_Mp);
	}
	return 1;
}


/*Main UI*/
int main()
{
	FILE *temp;
	time_t timep;
	struct tm *p;
	char pathName[50];
	int n=0,mode,len1,len2,pf;
	time(&timep);
	p=gmtime(&timep);

	printf("\nWhite Dwarfs computation code ver 1.0 for Linux\n========================================\n\n");
	temp=fopen("temp_dir","w");
	fclose(temp);
	system("dir ./EoS_lib >> temp_dir");
	if((temp=fopen("temp_dir","r"))==NULL)
	{
		printf("Unknown error, automatically quit.");fclose(temp);return 0;
	}else{
		printf("Equation of State files detected:\n\n");
	}
	while(fscanf(temp,"%s",list[n].filename)>0)
	{
		len2=strlen(list[n].filename);
		for(pf=0;pf<=len2-2;pf++)
		{
			if(list[n].filename[pf]=='.'){break;}
			list[n].name[pf]=list[n].filename[pf];
		}
		printf("%d\t%s\n",n+1,list[n].name);
		sprintf(list[n].path,"./EoS_lib/%s",list[n].filename);
		n++;
	}
	len1=n;
	fclose(temp);
	remove("temp_dir");
	printf("\nCaution: Please ensure all files above are available equation of states\n\n");
	Xout=fopen("BE_WD_report.xls","a");
	loadDenFile();
	for(pf=0;pf<len1;pf++)
	{
		loadEoS(list[pf].path);
		DenMode();
	}
	fclose(Xout);
	printf("Computation Done.\n");
	return 0;
}
