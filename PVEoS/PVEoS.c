#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define pi 3.14159265358979
#define Mscale 2.033931261665867e5
#define c 2.99792458e8

double E[]={7.0187,	7.418633,	7.818688,	8.218534999999999,	8.618676000000001,	9.018700000000001,	9.418633,	9.818754,	9.918711999999999,	10.218798,	10.518777,	10.81882,	11.118926,	11.418964,	11.51904,	11.719083,	11.91913,	12.019116,	12.119256,	12.219323,	12.419295,	12.619511,	12.81961,	12.919706,	13.019532,	13.119915,	13.219846,	13.320146,	13.420121,	13.520221,	13.620344,	13.72049,	13.820661,	13.920749,	14.020775,	14.121231,	14.221153,	14.265761,	14.321391,	14.421604,	14.521792,	14.622007,	14.633367,	14.825516,	15.001863,	15.223989,	15.826723,	16.128141,	16.304512,	16.429704,	16.52683,	16.606209,	16.673335,	16.731499,	16.782817,	16.828725,	16.870263,	16.908195,	16.943094,	16.975413,	17.005524};
double Pre[]={17.988737,	18.696182,	19.385785,	20.061075,	20.721481,	21.365113,	21.989227,	22.592288,	22.720903,	23.156852,	23.583539,	24.002598,	24.415641,	24.824516,	24.941412,	25.211921,	25.481299,	25.615845,	25.702086,	25.836324,	26.104487,	26.372175,	26.639686,	26.75297,	26.886604,	27.020361,	27.153815,	27.287354,	27.398461,	27.53199,	27.665393,	27.774444,	27.907895,	28.041393,	28.174641,	28.308137,	28.414472,	28.461198,	28.517196,	28.650599,	28.764624,	28.877256,	28.892373,	28.941472,	29.100646,	29.270166,	29.829799,	30.2512,	30.517791,	30.709075,	30.85755,	30.978838,	31.081527,	31.170731,	31.24981,	31.320999,	31.38589,	31.445682,	31.501211,	31.553179,	31.602091};
int length_NVBPS=sizeof(E)/sizeof(E[0]);

double Emax,p_surf,E_core,EOS_lgE[3000],EOS_lgp[3000];

/*double dM,dp,dy,dI,dAg00,temp_E,temp_cs2;
double temp_R,temp_M,temp_I,temp_k2,temp_Lambda;*/

//long Rec=0,strucpf=1;
int length_EOS=0;

double VsMIN[]={0.15,	0.185,	0.2,	0.23,	0.3,	0.4,	0.47,	0.5,	0.52,	0.55,	0.5,	0.4,	0.35,	0.32,	0.3};
double VsMAX[]={0.3,	0.42,	0.6,	0.7,	0.84,	0.98,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};
double EOSgrid[]={1.0,	1.5,	2.0,	2.5,	3.0,	3.5,	4.0,	4.5,	5.0,	5.5,	6.0,	6.5,	7.0,	8.0,	10.0};
int length_EOSgrid=sizeof(EOSgrid)/sizeof(EOSgrid[0]);

struct Parameters
{
	double VsSeq[50];
} Data;

//uniform distribution for Vs at each boundary
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
/*for(n=0;n<length_EOS;n++)
	{
		EOS_lgE[n]=EOS_lgE[n]-10.175607290470733;
		EOS_lgp[n]=EOS_lgp[n]-27.129849799910058;
	}
p_surf=pow(10,EOS_lgp[0]);*/
return 0;
}


//Main process
int main()
{
	FILE *outf;
	long pf;
	int n,eoslen;	
	char filename[30];

	printf("This is pre-generator of Piecewise Velocity Equation of States\n========================================\n\n");
	srand(time(NULL));

	printf("EoS counts: Initiate");
	for(pf=0;pf<100;pf++)
	{
		update_VsSeq();
		loadEoS();
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bEoS counts: %8ld",pf);
		eoslen=length_EOS;
		sprintf(filename,"%ld.txt",pf+1);
		outf=fopen(filename,"w");
		for(n=0;n<eoslen;n++)
		{
			fprintf(outf,"%lf\t%lf\n",EOS_lgE[n],EOS_lgp[n]);
		}
		fclose(outf);
	}
	printf("\n\nComputation Done!\n\nValid results: %ld\n\n",pf);
	printf("Attempt to make quickview...\n(need octave supported)\n");
	system("octave -q --persist quickView.m");
	return 0;
}
