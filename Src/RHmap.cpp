#include<stdio.h>
#include<stdlib.h>

class R2H{
	public:
	int ndat;	// number of data points
	double *rH;	// relative humidity
	double *hz;	// basal spacing [nm]
	int load(char fname[128]);
	double rH_min, rH_max;
	double hmin,hmax;
	int h_index(double h);
	int r_index(double r);
	double get_rH(double h);
	double get_hz(double r);

	private:
	protected:
};
int R2H::r_index(double r){
	if(r<=rH_min) return(-1);
	if(r>=rH_max) return(ndat);
	int indx=0;
	while(rH[indx+1]<r) indx++;
	return(indx);
}
double R2H::get_hz(double r){	// rH -> h mapping
	int i1=R2H::r_index(r);
	if(i1<0) return(hz[0]);
	if(i1>=ndat) return(hz[ndat-1]);

	double h1=hz[i1];
	double h2=hz[i1+1];
	double xi=(r-rH[i1])/(rH[i1+1]-rH[i1]);
	double eta=1.-xi;
	return(h1*eta+h2*xi);
};

int R2H::h_index(double h){
	if(h<=hmin) return(-1);
	if(h>=hmax) return(ndat);
	int indx=0;
	while(hz[indx+1]<h) indx++;
	return(indx);
};
double R2H::get_rH(double h){	// h -> rH mapping
	int i1=R2H::h_index(h);
	if(i1<0) return(rH[0]);
	if(i1>=ndat) return(rH[ndat-1]);

	double rH1=rH[i1];
	double rH2=rH[i1+1];
	double xi=(h-hz[i1])/(hz[i1+1]-hz[i1]);
	double eta=1.-xi;
	return(rH1*eta+rH2*xi);
};
int R2H::load(char fname[128]){
	FILE *fp=fopen(fname,"r");
	if(fp==NULL){
		printf("File %s cannot be found !\n",fname);
		printf(" --abort process\n");
		exit(-1);
	};
	char cbff[128];
	fscanf(fp,"%d\n",&ndat);
	fgets(cbff,128,fp);

	hz=(double *)malloc(sizeof(double)*ndat);
	rH=(double *)malloc(sizeof(double)*ndat);
	int i;
	double tmp;
	for(i=0;i<ndat;i++){
		fscanf(fp,"%lf, %lf, %lf\n",rH+i,hz+i,&tmp);
	};

	rH_min=rH[0];
	rH_max=rH[ndat-1];
	hmin=hz[0];
	hmax=hz[ndat-1];

	//printf("rH_lim=%lf %lf\n",rH_min,rH_max);
	//printf("hz_lim=%lf %lf\n",hmin,hmax);
};

int main(){
	char fname[128]="../rH_x.dat";
	R2H rh;
	rh.load(fname);

	double r1=-0.1; 
	double r2=1.1;
	int N=51;
	double dr=(r2-r1)/(N-1);

	double h1=9.0, h2=18.0;
	double dh=(h2-h1)/(N-1);

	int i;
	double r,h,hz,rH;
	for(i=0;i<N;i++){
		r=r1+dr*i;
		h=h1+dh*i;
		hz=rh.get_hz(r);
		rH=rh.get_rH(h);
		printf("%lf, %lf, %lf, %lf\n",r,hz, h,rH);
	};
	return(0);
};
