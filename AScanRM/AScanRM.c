#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define pi 3.14159265358979
#define Mscale 2.033931261665867e5
#define c 2.99792458e8

FILE *Pout,*Dataout,*strucpfout;
double dM,dp,dy,dI,dAg00,temp_E,temp_cs2;
double temp_R,temp_M,temp_I,temp_k2,temp_Lambda;
double p_surf,E_core,EOS_lgE[3000],EOS_lgp[3000];
double Ecentral[3000];
int length_EOS=0,length_DenFile=0,loopTOV=1,neos=1;
long Rec=0,strucpf=1;
char sT[30];
struct EoS
	{
		char path[60],filename[36],name[32];
	};
struct EoS list[100001];

struct Results
{
	int dataLength;
	double Mmax,R14,Lambda14,I14;
	double Ec[500],R[500],M[500],I[500],Lambda[500];
};
struct Results Data;

//linear interpolation
int interp_p2rhoE(cp)
double cp;
{
	if(cp<p_surf)
	{
		loopTOV=0;
		return 0;
	}
	double logp=log10(cp),logE;
	int n=length_EOS-1,np=length_EOS-1;
	for(;n>0;n--)
	{
		if(logp>EOS_lgp[n-1])
		{
			logE=(logp-EOS_lgp[n-1])*(EOS_lgE[n]-EOS_lgE[n-1])/(EOS_lgp[n]-EOS_lgp[n-1])+EOS_lgE[n-1];
			temp_E=pow(10,logE);
			temp_cs2=cp*(EOS_lgp[n]-EOS_lgp[n-1])/(temp_E*(EOS_lgE[n]-EOS_lgE[n-1]));
			return 1;
		}
	}
	return 0;
}

//interpolation from EnergyDensity to Pressure
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
return 0;
}

//Euler step
int fEul(fr,fm,fp,fy,fAg00,h)
double fr,fm,fp,fy,fAg00,h;
{
	double fpir=4.0*pi*fr*fr,hbt=1.0/sqrt(1-2*fm/fr),fE;
	interp_p2rhoE(fp);
	fE=temp_E;
	dM=fpir*fE*h;
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
	double h=1e-8,localRes,r,y;
	double m,p=interp_E2p(fEc),E,I,Ag00,Bg00;
	double m11,p11,y11,m21,p21,y21,m22,p22,y22,m31,p31,y31,m61,p61,y61,m62,p62,y62;
	double Ag0021,Ag0022,Ag0031,Ag0061,Ag0062,I21,I22,I31,I61,I62;
	double TR3m,TR3p,TR3I,TR3Ag00,TR3y,TR2m,TR2p,TR2y;
	double res1,res2,res3;
	r=h;
	interp_p2rhoE(p);
	E=temp_E;
	m=4.0/3*pi*r*r*r*E;
	y=2.0;
	I=0.0;
	Ag00=1.0;
	while(r<1e-4)
	{
		fEul(r,m,p,y,Ag00,h);
		m11=m+dM;
		p11=p+dp;
		y11=y+dy;
		
		m21=m+0.5*dM;
		p21=p+0.5*dp;
		y21=y+0.5*dy;
		I21=I+0.5*dI;
		Ag0021=Ag00+0.5*dAg00;
		
		m31=m+dM/3;
		p31=p+dp/3;
		y31=y+dy/3;
		I31=I+dI/3;
		Ag0031=Ag00+dAg00/3;
		
		m61=m+dM/6;
		p61=p+dp/6;
		y61=y+dy/6;
		I61=I+dI/6;
		Ag0061=Ag00+dAg00/6;
		
		fEul(r+0.5*h,m21,p21,y21,Ag0021,0.5*h);
		m22=m21+dM;
		p22=p21+dp;
		y22=y21+dy;
		I22=I21+dI;
		Ag0022=Ag0021+dAg00;
		
		fEul(r+h/6,m61,p61,y61,Ag0061,h/6);
		m62=m61+dM;
		p62=p61+dp;
		y62=y61+dy;
		I62=I61+dI;
		Ag0062=Ag0061+dAg00;
		
		if(loopTOV<1) break;
		
		TR3m=m22+9*(m62-m31);
		TR3p=p22+9*(p62-p31);
		TR3y=y22+9*(y62-y31);
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
			I=TR3I;
			Ag00=TR3Ag00;
			h=h*2;
		}else if(localRes>2e-5 && h>1.25e-9){
			h=h*0.5;
		}else{
			r=r+h;
			m=TR3m;
			p=TR3p;
			y=TR3y;
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
		fscanf(inf,"%lf%*[^\n]",&EOS_lgp[length_EOS]);
		length_EOS++;
	}
	fclose(inf);
	E_core=pow(10,EOS_lgE[length_EOS-1]);
	
	for(n=0;n<length_EOS;n++)
	{
		EOS_lgE[n]=EOS_lgE[n]-10.175607290470733;
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
		if(Ec<1e13)
		{
			printf("warning: incorrect density setting.");
			break;
		}
		Ecentral[length_DenFile]=Ec;
		length_DenFile++;
	}
	fclose(DenFile);
	return 1;
}


/*Working Mode: Automatically scan*/
int ScanMode()
{
	double Ec=1e18,dE=5e16,lastM,lastR,lastLambda,lastI;
	int n=0,pf=0;
	Data.dataLength=0;
	Data.Mmax=0;
	getRM(Ec);
	while(temp_M>0.75)
	{
		Ec=Ec*0.87;
		getRM(Ec);
	}
	Data.Ec[0]=Ec;Data.M[0]=temp_M;Data.R[0]=temp_R;Data.I[0]=temp_I;Data.Lambda[0]=temp_Lambda;(Data.dataLength)++;
	n++;
	lastM=temp_M;lastR=temp_R,lastLambda=temp_Lambda,lastI=temp_I;

	while(Ec<E_core)
	{
		Ec=Ec+dE;
		getRM(Ec);
		if(temp_M<lastM)
		{
			Data.Mmax=lastM;
			break;
		}else if((temp_M-lastM>0.02)||fabs(temp_R-lastR)>0.10)
		{
			dE=dE*0.5;
		}else if((temp_M-lastM<0.01)&&fabs(temp_R-lastR)<0.05)
		{
			dE=dE*2;
		}
		Data.Ec[n]=Ec;Data.M[n]=temp_M;Data.R[n]=temp_R;Data.I[n]=temp_I;Data.Lambda[n]=temp_Lambda;(Data.dataLength)++;
		n++;
		if(temp_M>=1.4&&lastM<=1.4)
		{
			Data.R14=lastR+(1.4-lastM)/(temp_M-lastM)*(temp_R-lastR);
			Data.Lambda14=lastLambda+(1.4-lastM)/(temp_M-lastM)*(temp_Lambda-lastLambda);
			Data.I14=lastI+(1.4-lastM)/(temp_M-lastM)*(temp_I-lastI);
		}
		lastM=temp_M;lastR=temp_R,lastLambda=temp_Lambda,lastI=temp_I;
	}
	Data.Mmax=lastM;
	if(Data.Mmax>1.0)
	{
		fprintf(Pout,"%lf\t%lf\t%lf\t%lf\n",Data.Mmax,Data.R14,Data.I14,Data.Lambda14);
		fprintf(strucpfout,"%ld\n",strucpf);
		for(pf=0;pf<Data.dataLength;pf++)
		{
			fprintf(Dataout,"%lf\t%lf\t%lf\t%lf\n",Data.R[pf],Data.M[pf],Data.Lambda[pf],Data.I[pf]);
		}
		Rec++;
		strucpf=strucpf+Data.dataLength;
	}

	return 1;
}


/*Working Mode: Density File Dependent*/
int DenMode()
{
	int pf;
	for(pf=0;pf<length_DenFile;pf++)
	{
		getRM(Ecentral[pf]);
		fprintf(Dataout,"%lf\t%lf\t%lf\t%lf\n",temp_R,temp_M,temp_Lambda,temp_I);
	}
	return 1;
}


/*Main UI*/
int main()
{
	FILE *temp;
	char pathName[50];
	int n=0,mode,len1,len2=60,pf;

	printf("Welcome to use Auto-Scan Radii and Masses --- version 9.05\n========================================\n\n");
	Pout=fopen("Results.txt","w");
	Dataout=fopen("InteralData.txt","w");
	strucpfout=fopen("DataStrucPf","w");
	system("cd EoS_lib && dir *.txt > EoS_Path.dir");
	if((temp=fopen("./EoS_lib/EoS_Path.dir","r"))==NULL)
	{
		printf("Unknown error, automatically quit.");fclose(temp);return 0;
	}else{
		printf("Equation of State files detected:\n\n");
	}
	while(fscanf(temp,"%s",list[n].filename)>0)
	{
		for(pf=0;pf<=len2;pf++)
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
	printf("\nCaution: Please ensure all *.txt files above are available equation of states!\n\nPress any key to continue\n");
	getchar();
	printf("\n\nSelect a work mode:\t\n[0]\tAuto Scan\t\n[1]\tDensity File\n\nEnter a Mode serial number: ");
	fflush(stdin);scanf("%d",&mode);
	if(mode==0)
	{
		for(pf=0;pf<len1;pf++)
		{
			loadEoS(list[pf].path);
			ScanMode();
		}
	}else{
		loadDenFile();
		for(pf=0;pf<len1;pf++)
		{
			loadEoS(list[pf].path);
			DenMode();
		}
	}
	fclose(Pout);
	fclose(Dataout);
	fclose(strucpfout);
	printf("\n\nComputation Done!\n\nValid results: %ld\n",Rec);
	return 0;
}
