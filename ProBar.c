#include <stdio.h>
#include <stdlib.h>

double *Rcano,*Lih,x[20],y[20];
long LenArr=0;

int deleteAarry(pf)
long pf;
{
	long n,fLen=--LenArr;
	for(n=pf;n<fLen;n++)
	{
		Rcano[n]=Rcano[n+1];
		Lih[n]=Lih[n+1];
	}
	return 1;
}

int loadResults()
{
	FILE *inf;
	double temp;
	inf=fopen("Results.txt","r");
	LenArr=0;
	while(fscanf(inf,"%lf%lf%lf%lf%lf%lf%*[^\n]",&temp,&temp,&temp,&Lih[LenArr],&temp,&Rcano[LenArr])>0)
	{
		LenArr++;
	}
	fclose(inf);
	return 1;
}

int initiateBar()
{
	int pf;
	for(pf=0;pf<20;pf++)
	{
		x[pf]=9.51+0.21*pf;
		y[pf]=0;
	}
	return 1;
}


int main()
{
	FILE *outf;
	double hx,Normfactor=0;
	long n,t,k;
	Rcano=(double*)malloc(4000000*sizeof(double));
	Lih=(double*)malloc(4000000*sizeof(double));
	printf("\nImporting files...\n");
	loadResults();
	initiateBar();
	hx=0.5*(x[1]-x[0]);
	n=LenArr;
	k=-1;
	printf("\nsorting...      ");
	while(n>0 && k<20)
	{
		t=0;k++;
		printf("\b\b\b\b\b\b%5ld%%",k*5);
		while(t<n)
		{
			if((Rcano[t]>x[k]?(Rcano[t]-x[k]):(x[k]-Rcano[t]))<hx)
			{
				y[k]=y[k]+Lih[t];
				deleteAarry(t);
				n--;
			}else{
				t++;
			}
		}
	}
	for(n=0;n<20;n++)
	{
		Normfactor=Normfactor+y[n];
	}
	
	outf=fopen("BarData.txt","w");	
	for(n=0;n<20;n++)
	{
		y[n]=y[n]/Normfactor;
		fprintf(outf,"%lf\t%lf\n",x[n],y[n]);
	}
	fclose(outf);
	
	
	printf("\n\nDone.");
	free(Rcano);
	free(Lih);
	return 0;
}