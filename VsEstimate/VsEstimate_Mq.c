#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define pi 3.14159265358979
#define Mscale 2.033931261665867e5
#define c 2.99792458e8

//NV+BPS EoS log10 data in SI units
double E[]={7.0187,	7.418633,	7.818688,	8.218534999999999,	8.618676000000001,	9.018700000000001,	9.418633,	9.818754,	9.918711999999999,	10.218798,	10.518777,	10.81882,	11.118926,	11.418964,	11.51904,	11.719083,	11.91913,	12.019116,	12.119256,	12.219323,	12.419295,	12.619511,	12.81961,	12.919706,	13.019532,	13.119915,	13.219846,	13.320146,	13.420121,	13.520221,	13.620344,	13.72049,	13.820661,	13.920749,	14.020775,	14.121231,	14.221153,	14.265761,	14.321391,	14.421604,	14.521792,	14.622007,	14.633367,	14.825516,	15.001863,	15.223989,	15.826723,	16.128141,	16.304512,	16.429704,	16.52683,	16.606209,	16.673335,	16.731499,	16.782817,	16.828725,	16.870263,	16.908195,	16.943094,	16.975413,	17.005524};
double Pre[]={17.988737,	18.696182,	19.385785,	20.061075,	20.721481,	21.365113,	21.989227,	22.592288,	22.720903,	23.156852,	23.583539,	24.002598,	24.415641,	24.824516,	24.941412,	25.211921,	25.481299,	25.615845,	25.702086,	25.836324,	26.104487,	26.372175,	26.639686,	26.75297,	26.886604,	27.020361,	27.153815,	27.287354,	27.398461,	27.53199,	27.665393,	27.774444,	27.907895,	28.041393,	28.174641,	28.308137,	28.414472,	28.461198,	28.517196,	28.650599,	28.764624,	28.877256,	28.892373,	28.941472,	29.100646,	29.270166,	29.829799,	30.2512,	30.517791,	30.709075,	30.85755,	30.978838,	31.081527,	31.170731,	31.24981,	31.320999,	31.38589,	31.445682,	31.501211,	31.553179,	31.602091};
int length_NVBPS=sizeof(E)/sizeof(E[0]);

//Sound Speed relevant EOS setting
double VsMIN[]={0.15,	0.185,	0.2,	0.23,	0.3,	0.4,	0.47,	0.5,	0.52,	0.55,	0.5,	0.4,	0.35,	0.32,	0.3},	VsMAX[]={0.3,	0.42,	0.6,	0.7,	0.84,	0.98,	1,	1,	1,	1,	1,	1,	1,	1,	1};
double EOSgrid[]={1,	1.5,	2,	2.5,	3,	3.5,	4,	4.5,	5,	5.5,	6,	6.5,	7,	8,	10};
int length_EOSgrid=sizeof(EOSgrid)/sizeof(EOSgrid[0]);

double Emax,p_surf,E_core,EOS_lgE[3000],EOS_lgp[3000];
double dM,dp,dy,dI,dAg00,temp_E,temp_cs2;
double temp_R,temp_M,temp_I,temp_k2,temp_Lambda;

long Rec=0,strucpf=1;
int length_EOS=0,loopTOV=1;

struct Results
{
	int dataLength;
	double VsSeq[50],Mmax,R14,Lambda14,I14;
	double Ec[500],R[500],M[500],I[500],Lambda[500];
};

struct Results Data;

double ComLambda[][2]=
{
0,	0,
20,	0,
44.303,	15,
85.192,	65,
126.08,	166,
166.97,	331,
207.86,	446,
248.75,	402,
289.63,	317,
330.52,	215,
371.41,	140,
412.3,	127,
453.19,	167,
494.08,	185,
534.97,	215,
575.85,	244,
616.74,	184,
657.63,	147,
698.52,	123,
739.41,	94,
780.3,	100,
821.18,	81,
862.07,	68,
902.96,	47,
943.85,	21,
984.74,	14,
1025.6,	8,
1066.5,	3,
1107.4,	6,
1148.3,	0,
1189.2,	2,
1230.1,	2,
1271,	2,
1311.8,	3,
1352.7,	1,
1393.6,	4,
1434.5,	2,
1475.4,	2,
1516.3,	0,
1557.2,	1,
1598.1,	1,
1639,	1,
1650,	0,
1700,	0
};
int length_ComLambda=sizeof(ComLambda)/sizeof(ComLambda[0][0])/2;

int update_VsSeq()
{
	int pf,npara=length_EOSgrid;
	double R_0_1;
	for(pf=0;pf<npara;pf++)
	{
		R_0_1= (double) rand()/RAND_MAX;
		Data.VsSeq[pf]=VsMIN[pf]+(VsMAX[pf]-VsMIN[pf])*R_0_1;
	}
	return 1;
}

double PDF_ComLambda(Lda)
double Lda;
{
	int pf,n=length_ComLambda;
	for(pf=1;pf<n;pf++)
	{
		if(Lda<ComLambda[pf][0])
		{
			return((Lda-ComLambda[pf-1][0])*(ComLambda[pf][1]-ComLambda[pf-1][1])/(ComLambda[pf][0]-ComLambda[pf-1][0])+ComLambda[pf-1][1]);
		}
	}
	return 0;
}

double m1Seq[5500],m2Seq[5500];
int length_massSeq=0;

int loadMassSeq()
{
	FILE *inf;
	int n;
	if((inf=fopen("massSeq","r"))==NULL)
	{
		printf("missing file massSeq at root path!\nCaution: this *.exe file cannot run independently.");
		return -1;
	}
	length_massSeq=0;
	while(fscanf(inf,"%lf%lf",&m1Seq[length_massSeq],&m2Seq[length_massSeq])>1){length_massSeq++;}
	fclose(inf);
	return 0;
}

//linear interpolation
int interp_p2E(cp)
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
	interp_p2E(fp);
	fE=temp_E;
	dM=fpir*fE*h;
	dp=h*(fE+fp)*(fm+4*pi*fr*fr*fr*fp)/(2*fm*fr-fr*fr);
	dI=8.0/3*h*pi*fr*fr*fr*fr*(fE+fp)*hbt/sqrt(fAg00);
	dAg00=-2.0*fAg00*dp/(fE+fp);
	dy=h*(-fy*fy/fr-fy/fr*(1+fpir*(fp-fE))/(1-2*fm/fr)-fr*(4*pi*(5*fE+9*fp+(fp+fE)/temp_cs2)/(1-2*fm/fr)-6/fr/fr/(1-2*fm/fr)-4/fr/fr/fr/fr/(1-2*fm/fr)/(1-2*fm/fr)*(fm+4*pi*fr*fr*fr*fp)*(fm+4*pi*fr*fr*fr*fp)));
	return 1;
}

//Main iteration process
int getRM(Ec)
double Ec;
{
	loopTOV=1;
	double fEc=Ec*6.6741e-11;
	double h=1e-8,localRes,r,y;
	double m,p=interp_E2p(fEc),E,I,Ag00,Bg00;
	double m11,p11,y11,m21,p21,y21,m22,p22,y22;
	double Ag0011,I11,Ag0021,Ag0022,I21,I22;
	double TR2m,TR2p,TR2I,TR2Ag00,TR2y;
	double res1,res2,res3;
	r=h;
	interp_p2E(p);
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
		I11=I+dI;
		Ag0011=Ag00+dAg00;

		m21=m+0.5*dM;
		p21=p+0.5*dp;
		y21=y+0.5*dy;
		I21=I+0.5*dI;
		Ag0021=Ag00+0.5*dAg00;

		fEul(r+0.5*h,m21,p21,y21,Ag0021,0.5*h);
		m22=m21+dM;
		p22=p21+dp;
		y22=y21+dy;
		I22=I21+dI;
		Ag0022=Ag0021+dAg00;

		if(loopTOV<1) break;

		TR2m=2*m22-m11;
		TR2p=2*p22-p11;
		TR2y=2*y22-y11;
		TR2I=2*I22-I11;
		TR2Ag00=2*Ag0022-Ag0011;

		res1=(TR2m-m11)/(TR2m-m);
		res1=res1>0?res1:-res1;
		res2=(TR2p-p11)/(TR2p-p);
		res2=res2>0?res2:-res2;
		res3=(TR2y-y11)/(TR2y-y);
		res3=res3>0?res3:-res3;
		res2=(res2>res3)?res2:res3;
		localRes=(res1>res2)?res1:res2;

		if(localRes<1e-5)
		{
			r=r+h;
			m=TR2m;
			p=TR2p;
			y=TR2y;
			I=TR2I;
			Ag00=TR2Ag00;
			h=h*2;
		}else if(localRes>5e-4 && h>2e-10){
			h=h*0.5;
		}else{
			r=r+h;
			m=TR2m;
			p=TR2p;
			y=TR2y;
			I=TR2I;
			Ag00=TR2Ag00;
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


//interpolate from m to Lambda
double interp_L(m)
double m;
{
	int pf,n=Data.dataLength-1;
	for(pf=n;pf>0;pf--)
	{
		if(m>Data.M[pf-1])
		{
			return((m-Data.M[pf-1])*(Data.Lambda[pf]-Data.Lambda[pf-1])/(Data.M[pf]-Data.M[pf-1])+Data.Lambda[pf-1]);
		}
	}
	return 4000;
}

double getLih()
{
	double Lih=0,m1,m2,L1,L2;
	int pf,n=length_massSeq;
	for(pf=0;pf<n;pf++)
	{
		m1=m1Seq[pf];
		m2=m2Seq[pf];
		L1=interp_L(m1);
		L2=interp_L(m2);
		Lih=Lih+PDF_ComLambda(16.0/13*((m1+12*m2)*(m1*m1*m1*m1)*L1+(m2+12*m1)*(m2*m2*m2*m2)*L2)/pow((m1+m2),5));
	}
	return Lih;
}


//Working Mode: Automatically scan
int ScanMode()
{
	FILE *Vout,*Pout;
	double Ec=7e17,dE=5e16,lastM,lastR,lastLambda,lastI,Lih;
	int n=0,pf=0,npara=length_EOSgrid;

	Data.dataLength=0;
	Data.Mmax=0;
	getRM(Ec);
	while(temp_M>1.0)
	{
		Ec=Ec*0.87;
		getRM(Ec);
		pf++;
		if(pf>10) {return -1;}
	}
	Data.Ec[0]=Ec;Data.M[0]=temp_M;Data.R[0]=temp_R;Data.I[0]=temp_I;Data.Lambda[0]=temp_Lambda;(Data.dataLength)++;
	n++;
	lastM=temp_M;lastR=temp_R,lastLambda=temp_Lambda,lastI=temp_I;
	while(Ec<E_core)
	{
		Ec=Ec+dE;
		getRM(Ec);
		if((temp_M<lastM)||(temp_M>2.5))
		{
			Data.Mmax=lastM;
			break;
		}else if(temp_M-lastM>0.02)
		{
			dE=dE*0.5;
		}else if(temp_M-lastM<0.01)
		{
			dE=dE*2;
		}
		Data.Ec[n]=Ec;Data.M[n]=temp_M;Data.R[n]=temp_R;Data.I[n]=temp_I;Data.Lambda[n]=temp_Lambda;(Data.dataLength)++;
		n++;
		if((temp_M>=1.4)&&(lastM<=1.4))
		{
			Data.R14=lastR+(1.4-lastM)/(temp_M-lastM)*(temp_R-lastR);
			Data.Lambda14=lastLambda+(1.4-lastM)/(temp_M-lastM)*(temp_Lambda-lastLambda);
			Data.I14=lastI+(1.4-lastM)/(temp_M-lastM)*(temp_I-lastI);
		}
		lastM=temp_M;lastR=temp_R,lastLambda=temp_Lambda,lastI=temp_I;
	}
	Data.Mmax=lastM;
	if(Data.Mmax>2.01)
	{
		Lih=getLih();
		Vout=fopen("VsPara","a");
		for(pf=0;pf<npara;pf++)
		{
			fprintf(Vout,"%lf\t",Data.VsSeq[pf]);
		}
		fprintf(Vout,"\n");
		fclose(Vout);
		Pout=fopen("Results","a");
		fprintf(Pout,"%lf\t%lf\t%lf\t%lf\t%lf\n",Lih,Data.Mmax,Data.R14,Data.I14,Data.Lambda14);
		fclose(Pout);
		Rec++;
	}
	return 1;
}


//load ANM EOS with paras
int loadEoS()
{
	length_EOS=length_NVBPS;
	double VsSeqSI[50],ESeqSI[50],PSeqSI[50],dESI,diffVs,dvE;
	int pf,n;

	for(pf=0;pf<length_NVBPS;pf++)
	{
			EOS_lgE[pf]=E[pf];
			EOS_lgp[pf]=Pre[pf];
	}
	ESeqSI[0]=pow(10,E[length_NVBPS-1]);
	PSeqSI[0]=pow(10,Pre[length_NVBPS-1]);
	VsSeqSI[0]=sqrt(PSeqSI[0]*(Pre[length_NVBPS-1]-Pre[length_NVBPS-2])/(ESeqSI[0]*(E[length_NVBPS-1]-E[length_NVBPS-2])));

	for(pf=0;pf<length_EOSgrid;pf++)
	{
		ESeqSI[pf+1]=EOSgrid[pf]*2.68e17;
		VsSeqSI[pf+1]=c*Data.VsSeq[pf];
		PSeqSI[pf+1]=PSeqSI[pf]+(ESeqSI[pf+1]-ESeqSI[pf])*(VsSeqSI[pf+1]*VsSeqSI[pf+1]+VsSeqSI[pf]*VsSeqSI[pf+1]+VsSeqSI[pf]*VsSeqSI[pf])/3;
		dESI=(ESeqSI[pf+1]-ESeqSI[pf])/16;
		diffVs=(VsSeqSI[pf+1]-VsSeqSI[pf])/(ESeqSI[pf+1]-ESeqSI[pf]);
		for(n=0;n<16;n++)
		{
			dvE=(n+1)*dESI;
			EOS_lgE[length_EOS]=log10(ESeqSI[pf]+dvE);
			EOS_lgp[length_EOS]=log10(diffVs*diffVs*dvE*dvE*dvE/3+VsSeqSI[pf]*diffVs*dvE*dvE+VsSeqSI[pf]*VsSeqSI[pf]*dvE+PSeqSI[pf]);
			length_EOS++;
		}
	}
E_core=pow(10,EOS_lgE[length_EOS-1]);
for(n=0;n<length_EOS;n++)
	{
		EOS_lgE[n]=EOS_lgE[n]-10.175607290470733;
		EOS_lgp[n]=EOS_lgp[n]-27.129849799910058;
	}
p_surf=pow(10,EOS_lgp[0]);
return 0;
}


//Main process
int main()
{
	long n=0,pf;

	printf("Welcome to use Vs Estimation --- version Mq\n========================================\n\n");

	srand(time(NULL));

	printf("Loading mass frames...\n");
	loadMassSeq();
	printf("\nComputing...\n\n");
	for(pf=0;pf<100;pf++)
	{
		update_VsSeq();
		loadEoS();
		printf("EoS counts: %8ld\n",++n);//\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b
		ScanMode();
	}
	printf("\n\nComputation Done!\n\nValid results: %ld\n",Rec);
	return 0;
}
