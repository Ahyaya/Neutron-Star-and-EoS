#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define pi 3.14159265358979
#define XM 939.0
#define XME 0.511
#define XMMU 105.658
#define Rho0 0.16
#define H 1239.8522
#define HBAR 197.3286
#define Mscale 2.033931261665867e5
#define c 2.99792458e8

FILE *Pout;
/*double E[]={1.044000e+004	,2.622000e+004	,6.587000e+004	,1.654000e+005	,4.156000e+005	,1.044000e+006	,2.622000e+006	,6.588000e+006	,8.293000e+006	,1.655000e+007	,3.302000e+007	,6.589000e+007	,1.315000e+008	,2.624000e+008	,3.304000e+008	,5.237000e+008	,8.301000e+008	,1.045000e+009	,1.316000e+009	,1.657000e+009	,2.626000e+009	,4.164000e+009	,6.601000e+009	,8.312000e+009	,1.046000e+010	,1.318000e+010	,1.659000e+010	,2.090000e+010	,2.631000e+010	,3.313000e+010	,4.172000e+010	,5.254000e+010	,6.617000e+010	,8.332000e+010	,1.049000e+011	,1.322000e+011	,1.664000e+011	,1.844000e+011	,2.096000e+011	,2.640000e+011	,3.325000e+011	,4.188000e+011	,4.299000e+011	,6.691390e+011	,1.004300e+012	,1.674900e+012	,6.710000e+012	,1.343200e+013	,2.016100e+013	,2.689700e+013	,3.363800e+013	,4.038400e+013	,4.713410e+013	,5.388890e+013	,6.064810e+013	,6.741010e+013	,7.417600e+013	,8.094600e+013	,8.771900e+013	,9.449600e+013	,1.012800e+014	,1.080600e+014	,1.148400e+014	,1.216300e+014	,1.284300e+014	,1.352200e+014	,1.518900e+014	,1.689300e+014	,1.860100e+014	,2.031390e+014	,2.202900e+014	,2.374910e+014	,2.547500e+014};
double Pre[]={9.744000e+018	,4.968000e+019	,2.431000e+020	,1.151000e+021	,5.266000e+021	,2.318000e+022	,9.755000e+022	,3.911000e+023	,5.259000e+023	,1.435000e+024	,3.833000e+024	,1.006000e+025	,2.604000e+025	,6.676000e+025	,8.738000e+025	,1.629000e+026	,3.029000e+026	,4.129000e+026	,5.036000e+026	,6.860000e+026	,1.272000e+027	,2.356000e+027	,4.362000e+027	,5.662000e+027	,7.702000e+027	,1.048000e+028	,1.425000e+028	,1.938000e+028	,2.503000e+028	,3.404000e+028	,4.628000e+028	,5.949000e+028	,8.089000e+028	,1.100000e+029	,1.495000e+029	,2.033000e+029	,2.597000e+029	,2.892000e+029	,3.290000e+029	,4.473000e+029	,5.816000e+029	,7.538000e+029	,7.805000e+029	,8.739200e+029	,1.260800e+030	,1.862800e+030	,6.757700e+030	,1.783200e+031	,3.294510e+031	,5.117700e+031	,7.203600e+031	,9.524410e+031	,1.206500e+032	,1.481600e+032	,1.777500e+032	,2.094110e+032	,2.431590e+032	,2.790500e+032	,3.171110e+032	,3.574200e+032	,4.000290e+032	,4.452010e+032	,4.924600e+032	,5.424390e+032	,5.950390e+032	,6.503300e+032	,1.165200e+033	,1.486500e+033	,1.869990e+033	,2.289900e+033	,2.780600e+033	,3.335110e+033	,3.950400e+033};
double Den[]={6.295000e-012	,1.581000e-011	,3.972000e-011	,9.976000e-011	,2.506000e-010	,6.294000e-010	,1.581000e-009	,3.972000e-009	,5.000000e-009	,9.976000e-009	,1.990000e-008	,3.972000e-008	,7.924000e-008	,1.581000e-007	,1.990000e-007	,3.155000e-007	,5.000000e-007	,6.294000e-007	,7.924000e-007	,9.976000e-007	,1.581000e-006	,2.506000e-006	,3.972000e-006	,5.000000e-006	,6.294000e-006	,7.924000e-006	,9.976000e-006	,1.256000e-005	,1.581000e-005	,1.990000e-005	,2.506000e-005	,3.155000e-005	,3.972000e-005	,5.000000e-005	,6.294000e-005	,7.924000e-005	,9.976000e-005	,1.105000e-004	,1.256000e-004	,1.581000e-004	,1.990000e-004	,2.506000e-004	,2.572000e-004	,4.000000e-004	,6.000000e-004	,1.000000e-003	,4.000000e-003	,8.000000e-003	,1.200000e-002	,1.600000e-002	,2.000000e-002	,2.400000e-002	,2.800000e-002	,3.200000e-002	,3.600000e-002	,4.000000e-002	,4.400000e-002	,4.800000e-002	,5.200000e-002	,5.600000e-002	,6.000000e-002	,6.400000e-002	,6.800000e-002	,7.200000e-002	,7.600000e-002	,8.000000e-002	,9.000000e-002	,1.000000e-001	,1.100000e-001	,1.200000e-001	,1.300000e-001	,1.400000e-001	,1.500000e-001};*/
double Emax,p_surf,E_core,EOS_lgDen[3000],EOS_lgE[3000],EOS_lgp[3000];
double KsymTotal[5500],JsymTotal[5500],J0Total[5500],RhoTTotal[5500],PTTotal[5500];
double dM,dMa,dMp,dp,dy,dI,dAg00,temp_rho,temp_E,temp_cs2;
double temp_R,temp_M,temp_I,temp_Ma,temp_Mp,temp_k2,temp_Lambda;
double temp_Ksym,temp_Jsym,temp_J0;
static long PtLen=0;
long Rec=0;
int length_EOS=0,loopTOV=1,neos;

struct Results
{
	int dataLength;
	double Ksym,Jsym,J0,Mmax,R14,Lambda14,I14;
	double Ec[400],R[400],M[400],I[400],Lambda[400],Ma[400],Mp[400];
};

struct Results Data;

struct EoS
{
	char Name[55],tName[50];
};

struct EoS list[30001];


double Data_m1[][2]=
{
	1.00000000000000,	0,
	1.36570000000000,	0,
	1.38570000000000,	338,
	1.40573000000000,	367,
	1.42576000000000,	311,
	1.44579000000000,	348,
	1.46582000000000,	360,
	1.48585000000000,	332,
	1.50588000000000,	324,
	1.52591000000000,	315,
	1.54594000000000,	278,
	1.56597000000000,	219,
	1.58600000000000,	204,
	1.60603000000000,	153,
	1.62606000000000,	140,
	1.64609000000000,	99,
	1.66612000000000,	60,
	1.68615000000000,	46,
	1.70618000000000,	25,
	1.72621000000000,	25,
	1.74624000000000,	6,
	1.76627000000000,	2,
	1.78630000000000,	0,
	1.80000000000000,	0
};

double Data_m2[][2]=
{
	1.00000000000000,	0,
	1.07110000000000,	0,
	1.08592751835859,	4,
	1.10077829068030,	16,
	1.11562906300201,	35,
	1.13047983532372,	51,
	1.14533060764543,	73,
	1.16018137996714,	132,
	1.17503215228885,	153,
	1.18988292461056,	179,
	1.20473369693227,	226,
	1.21958446925398,	239,
	1.23443524157569,	282,
	1.24928601389741,	295,
	1.26413678621912,	298,
	1.27898755854083,	315,
	1.29383833086254,	302,
	1.30868910318425,	289,
	1.32353987550596,	241,
	1.33839064782767,	284,
	1.35324142014938,	280,
	1.36809219247109,	258,
	1.38290000000000,	0,
	1.40000000000000,	0
};

double Data_Lambda1[][2]=
{
	-32.000000000000,	0,
	-30.427000000000,	0,
	30.4757916894131,	546,
	91.3786130549647,	517,
	152.281434420516,	494,
	213.184255786068,	403,
	274.087077151619,	347,
	334.989898517171,	285,
	395.892719882723,	210,
	456.795541248274,	201,
	517.698362613826,	144,
	578.601183979377,	121,
	639.504005344929,	97,
	700.406826710481,	104,
	761.309648076032,	77,
	822.212469441584,	82,
	883.115290807135,	62,
	944.018112172687,	59,
	1004.92093353824,	47,
	1065.82375490379,	29,
	1126.72657626934,	26,
	1187.62939763489,	27,
	1248.53221900044,	18,
	1309.43504036600,	15,
	1370.33786173155,	14,
	1431.24068309710,	6,
	1492.14350446265,	8,
	1553.04632582820,	3,
	1613.94914719375,	1,
	1674.85196855931,	1,
	1735.75478992486,	2,
	1796.65761129041,	3,
	1857.56043265596,	0,
	1918.46325402151,	1,
	1979.36607538706,	0,
	2040.26889675261,	0,
	2101.17171811817,	1,
	2162.07453948372,	0,
	2222.97736084927,	0,
	2283.88018221482,	0,
	2344.78300358037,	0,
	2405.68582494592,	1,
	2466.66666666666,	0,
	2500.00000000000,	0
};

double Data_Lambda2[][2]=
{
	-48.000000000000,	0,
	-47.007800000000,	0,
	47.3240916547122,	463,
	141.656018702800,	419,
	235.987945750887,	457,
	330.319872798975,	431,
	424.651799847062,	389,
	518.983726895150,	314,
	613.315653943237,	253,
	707.647580991325,	191,
	801.979508039412,	181,
	896.311435087500,	136,
	990.643362135587,	129,
	1084.97528918367,	90,
	1179.30721623176,	102,
	1273.63914327985,	97,
	1367.97107032794,	75,
	1462.30299737602,	53,
	1556.63492442411,	37,
	1650.96685147220,	28,
	1745.29877852029,	26,
	1839.63070556837,	19,
	1933.96263261646,	21,
	2028.29455966455,	7,
	2122.62648671264,	10,
	2216.95841376073,	7,
	2311.29034080881,	3,
	2405.62226785690,	1,
	2499.95419490499,	2,
	2594.28612195307,	3,
	2688.61804900116,	2,
	2782.94997604925,	0,
	2877.28190309734,	2,
	2971.61383014543,	0,
	3065.94575719351,	0,
	3160.27768424160,	1,
	3254.60961128969,	0,
	3348.94153833777,	1,
	3443.27346538586,	0,
	3537.60539243395,	0,
	3631.93731948204,	0,
	3726.26924653012,	2,
	3820.60000000000,	0,
	4000.00000000000,	0
};

int LenData_m1=sizeof(Data_m1)/sizeof(Data_m1[0][0])/2;
int LenData_m2=sizeof(Data_m2)/sizeof(Data_m2[0][0])/2;
int LenData_Lambda1=sizeof(Data_Lambda1)/sizeof(Data_Lambda1[0][0])/2;
int LenData_Lambda2=sizeof(Data_Lambda2)/sizeof(Data_Lambda2[0][0])/2;

double PDF_m1(m)
double m;
{
	int pf,n=LenData_m1-1;
	for(pf=n;pf>0;pf--)
	{
		if(m>Data_m1[pf-1][0])
		{
			return((m-Data_m1[pf-1][0])*(Data_m1[pf][1]-Data_m1[pf-1][1])/(Data_m1[pf][0]-Data_m1[pf-1][0])+Data_m1[pf-1][1]);
		}
	}
	return 0;
}

double PDF_m2(m)
double m;
{
	int pf,n=LenData_m2-1;
	for(pf=n;pf>0;pf--)
	{
		if(m>Data_m2[pf-1][0])
		{
			return((m-Data_m2[pf-1][0])*(Data_m2[pf][1]-Data_m2[pf-1][1])/(Data_m2[pf][0]-Data_m2[pf-1][0])+Data_m2[pf-1][1]);
		}
	}
	return 0;
}

double PDF_Lambda1(Lda)
double Lda;
{
	int pf,n=LenData_Lambda1-1;
	for(pf=n;pf>0;pf--)
	{
		if(Lda>Data_Lambda1[pf-1][0])
		{
			return((Lda-Data_Lambda1[pf-1][0])*(Data_Lambda1[pf][1]-Data_Lambda1[pf-1][1])/(Data_Lambda1[pf][0]-Data_Lambda1[pf-1][0])+Data_Lambda1[pf-1][1]);
		}
	}
	return 0;
}

double PDF_Lambda2(Lda)
double Lda;
{
	int pf,n=LenData_Lambda2-1;
	for(pf=n;pf>0;pf--)
	{
		if(Lda>Data_Lambda2[pf-1][0])
		{
			return((Lda-Data_Lambda2[pf-1][0])*(Data_Lambda2[pf][1]-Data_Lambda2[pf-1][1])/(Data_Lambda2[pf][0]-Data_Lambda2[pf-1][0])+Data_Lambda2[pf-1][1]);
		}
	}
	return 0;
}

//linear interpolation
int interp_p2rhoE(cp)
double cp;
{
	if(cp<p_surf)
	{
		loopTOV=0;
		return 0;
	}
	double logp=log10(cp),logDen,logE;
	int n=length_EOS-1,np=length_EOS-1;
	for(;n>0;n--)
	{
		if(logp>EOS_lgp[n-1])
		{
			logDen=(logp-EOS_lgp[n-1])*(EOS_lgDen[n]-EOS_lgDen[n-1])/(EOS_lgp[n]-EOS_lgp[n-1])+EOS_lgDen[n-1];
			logE=(logp-EOS_lgp[n-1])*(EOS_lgE[n]-EOS_lgE[n-1])/(EOS_lgp[n]-EOS_lgp[n-1])+EOS_lgE[n-1];
			temp_rho=pow(10,logDen);
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

//Main iteration process
int getRM(Ec)
double Ec;
{
	loopTOV=1;
	double fEc=Ec*6.6741e-11;
	double h=1e-8,localRes,r,y;
	double m,ma,mp,p=interp_E2p(fEc),E,rho,I,Ag00,Bg00;
	double m11,p11,y11,m21,p21,y21,m22,p22,y22,m31,p31,y31,m61,p61,y61,m62,p62,y62;
	double ma21,ma22,ma31,ma61,ma62,mp21,mp22,mp31,mp61,mp62;
	double Ag0021,Ag0022,Ag0031,Ag0061,Ag0062,I21,I22,I31,I61,I62;
	double TR3m,TR3ma,TR3mp,TR3p,TR3I,TR3Ag00,TR3y,TR2m,TR2p,TR2y;
	double res1,res2,res3;
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
	while(r<1e-4)
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
		
		if(localRes<1e-5 && h<1e-7)
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
		}else if(localRes>2e-4 && h>2e-10){
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

double getLih_1()
{
	double m,dLih11,dLih21,dLih22,dLih31,dLih61,dLih62,dm=0.003,Lih=0;
	
	for(m=1.36;m<1.8;m=m+dm)
	{
		dLih11=PDF_m1(m)*PDF_Lambda1(interp_L(m))*dm;
		dLih21=dLih11/2;
		dLih31=dLih11/3;
		dLih61=dLih11/6;
		dLih22=dLih21+PDF_m1(m+dm/2)*PDF_Lambda1(interp_L(m+dm/2))*dm/2;
		dLih62=dLih61+PDF_m1(m+dm/6)*PDF_Lambda1(interp_L(m+dm/6))*dm/6;
		Lih=Lih+dLih22+9*(dLih62-dLih31);
	}
	return Lih;
}

double getLih_2()
{
	double m,dLih11,dLih21,dLih22,dLih31,dLih61,dLih62,dm=0.003,Lih=0;
	
	for(m=1.07;m<1.4;m=m+dm)
	{
		dLih11=PDF_m2(m)*PDF_Lambda2(interp_L(m))*dm;
		dLih21=dLih11/2;
		dLih31=dLih11/3;
		dLih61=dLih11/6;
		dLih22=dLih21+PDF_m2(m+dm/2)*PDF_Lambda2(interp_L(m+dm/2))*dm/2;
		dLih62=dLih61+PDF_m2(m+dm/6)*PDF_Lambda2(interp_L(m+dm/6))*dm/6;
		Lih=Lih+dLih22+9*(dLih62-dLih31);
	}
	return Lih;
}


//Working Mode: Automatically scan
int ScanMode()
{
	double Ec=1e18,dE=5e16,lastM,lastR,lastLambda,lastI,Lih;
	int n=0,pf=0;
	
	Data.dataLength=0;
	Data.Mmax=0;
	getRM(Ec);
	while(temp_M>1.0)
	{
		Ec=Ec*0.87;
		getRM(Ec);
	}
	//Data.Ksym=temp_Ksym;Data.Jsym=temp_Jsym;Data.J0=temp_J0;
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
	//if(Data.Mmax>2.01)
	{
		Lih=getLih_1()*getLih_2();
		fprintf(Pout,"%s\t%lf\t%lf\t%lf\t%lf\t%lf\n",list[neos-1].Name,Lih,Data.Mmax,Data.R14,Data.I14,Data.Lambda14);
		Rec++;
	}
	return 1;
}


int loadEoS(EOSname)
char EOSname[65];
{
	FILE *inf;
	int pf,n;
	if((inf=fopen(EOSname,"r"))==NULL)
	{
		printf("\nfile not found!\n");
		return 0;
	}
	length_EOS=0;
	while(fscanf(inf,"%lf",&EOS_lgE[length_EOS])==1)
	{
		fscanf(inf,"%lf%*[^\n]",&EOS_lgp[length_EOS]);
		//fscanf(inf,"%lf%*[^\n]",&EOS_lgDen[length_EOS]);
		length_EOS++;
	}
	fclose(inf);
	E_core=pow(10,EOS_lgE[length_EOS-1]);

	for(n=0;n<length_EOS;n++)
	{
		//EOS_lgDen[n]=EOS_lgDen[n]-10.175607290470733;
		EOS_lgE[n]=EOS_lgE[n]-10.175607290470733;
		EOS_lgp[n]=EOS_lgp[n]-27.129849799910058;
	}
	p_surf=pow(10,EOS_lgp[0]);
	return 0;

}	

//Main process
int main()
{
	int n=0;
	double Ksym,Jsym,J0;
	FILE *temp;
	char pathName[65];
	int len1,len2,pf;
	
	printf("Welcome to use ModMCMC --- version 1.0\n========================================\n\n");
	
	Pout=fopen("ModResults.txt","w");
	system("dir ./EoS_lib >> temp_dir");
	if((temp=fopen("temp_dir","r"))==NULL)
	{
		printf("No Equation of State files detected, process has been shutdown.");fclose(temp);return 0;
	}else{
		printf("Equation of State files detected:\n\n");
	}
	while(fscanf(temp,"%s",list[n].Name)>0)
	{
		len2=strlen(list[n].Name);
		for(pf=0;pf<=len2-2;pf++)
		{
			if(list[n].Name[pf]=='.'){break;}
			list[n].tName[pf]=list[n].Name[pf];
		}
		printf("%d\t%s\n",n+1,list[n].tName);
		n++;
	}
	len1=n;
	fclose(temp);
	remove("temp_dir");
	printf("\nCaution: Please ensure all *.txt files in above are available equation of states!\n\nPress any key to continue\n");
	getchar();
	printf("\nComputing...\n\n");
	for(neos=1;neos<=len1;neos++)
	{
		sprintf(pathName,"./EoS_lib/%s",list[neos-1].Name);
		loadEoS(pathName);
		ScanMode();
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bEoS counts: %8d",neos);
	}


	fclose(Pout);
	printf("\n\nComputation Done!\n\nValid results: %ld\n",Rec);
	return 0;
}
