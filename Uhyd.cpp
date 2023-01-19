#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//	N: molar number of H2O molecules
//	H: basal spacing 

//------------------------------------------------------------
class N2H{
	public:
		int ndat;	// number of nH2O data points (equi-spaced, prescribed)
		double *nH2O;   // number of H2O/unit structure 
		double *hz;	// basal spacing [nm] (N->H map)
		int load(char fname[128]); // load data from "NaMt.dat" or "CaMt.dat"
		double dn;	// nH2O increment (for MD simlation)
		double hmin,hmax; // min and max basal spacing 
		double nmin,nmax; // min and max number of H2O mole number
		double get_h(double nw); // return h(basal spacing) for given n(H2O)
		int h_index(double h); // return index for given h(basal spacing) 
		double get_nH2O(double h);// return n(H2O) for givne h

		int ndat_h;	// equi-spaced hz data points (arbitrary) 
		double dh0;	// increment in h(basal spacing) 
		double *nH2O_intpl; // interplated n(H2O) data ( H->N map: inverse function)
		void setup_H2N(int nn); // establish H->N mapping
		void fwrite(); // export N->H(fwd.) and H->N (inv.) mappings 

	private:
	protected:
};
//------------------------------------------------------------
double N2H::get_h(double nw){

	if(nw<0.0) return(hmin);
	if(nw>=nH2O[ndat-1]) return(hmax);

	int i1,i2;
	i1=int(nw/dn);
	i2=i1+1;
	double xi=nw/dn-i1;
	double eta=1.-xi;
	return(eta*hz[i1]+xi*hz[i2]);
};
int N2H::h_index(double h){
	if(h<=hmin) return(-1);
	if(h>=hmax) return(ndat);
	int indx=0;
	while(hz[indx+1]<h) indx++;
	return(indx);
};
double N2H::get_nH2O(double h){
	int i1=N2H::h_index(h);
	if(i1<0) return(nH2O[0]);
	if(i1>=ndat) return(nH2O[ndat-1]);

	double n1=nH2O[i1];
	double n2=nH2O[i1+1];
	double xi=(h-hz[i1])/(hz[i1+1]-hz[i1]);
	double eta=1.-xi;
	return(n1*eta+n2*xi);
};
void N2H::setup_H2N(int nn){
	ndat_h=nn;
	dh0=(hmax-hmin)/(ndat_h-1);
	nH2O_intpl=(double *)malloc(sizeof(double)*ndat_h);

	int i;
	double h;
	for(i=0;i<ndat_h;i++){
		h=hmin+dh0*i;
		nH2O_intpl[i]=N2H::get_nH2O(h);
	};
};
void N2H::fwrite(){
	FILE *fp1=fopen("n2h.dat","w");
	FILE *fp2=fopen("h2n.dat","w");

	int i;
	fprintf(fp1,"%d\n",ndat);
	fprintf(fp1,"# n(H2O), h(basal spacing [nm]) \n");
	for(i=0;i<ndat;i++){
		fprintf(fp1,"%lf, %lf\n",nH2O[i],hz[i]);
	};
	double h;
	fprintf(fp2,"%d\n",ndat_h);
	fprintf(fp2,"# n(H2O), h(basal spacing [nm]) \n");
	for(i=0;i<ndat_h;i++){
		h=hmin+dh0*i;
		fprintf(fp2,"%lf, %lf\n",nH2O_intpl[i],h);
	};

	fclose(fp1);
	fclose(fp2);

};

int N2H:: load(char fname[128]){
	FILE *fp=fopen(fname,"r");
	if(fp==NULL){
		printf("File %s cannot be found !\n",fname);
		printf(" --abort process\n");
		exit(-1);
	};

	char cbff[128];
	fgets(cbff,128,fp);

	// Count data lines
	ndat=0;
	while(fgets(cbff,128,fp)!=NULL){
		ndat++;
	};
	printf("#ndat=%d\n",ndat);
	fclose(fp);

	hz=(double *)malloc(sizeof(double)*ndat);
	nH2O=(double *)malloc(sizeof(double)*ndat);

	fp=fopen(fname,"r");
	fgets(cbff,128,fp);
	int i;
	double Un,Nc,dP;
	for(i=0;i<ndat;i++){
		fscanf(fp,"%lf, %lf, %lf, %lf, %lf\n",nH2O+i,&Un,hz+i,&Nc,&dP);
		hz[i]*=1.e-01;	// chage unit to [nm]
	};
	fclose(fp);

	dn=nH2O[1]-nH2O[0];
	hmin=hz[0];
	hmax=hz[ndat-1];
	nmin=nH2O[0];
	nmax=nH2O[ndat-1];
	return(ndat);
};
//-----------------------------------------------------------
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
//-----------------------------------------------------------
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
		hz[i]*=1.e-01;
	};

	rH_min=rH[0];
	rH_max=rH[ndat-1];
	hmin=hz[0];
	hmax=hz[ndat-1];

	//printf("rH_lim=%lf %lf\n",rH_min,rH_max);
	//printf("hz_lim=%lf %lf\n",hmin,hmax);
};

//-----------------------------------------------------------
int main(){

	char fname1[128]="MD_Data/NaMt.dat";
	char fname2[128]="rH_x.dat";

	N2H n2h_MD;	// MD-generated swelling curve
	n2h_MD.load(fname1);

	R2H r2h_XRD;	// XRD-meassured swelling curve
	r2h_XRD.load(fname2);

/*
	int i;
	double n1=n2h_MD.nmin;
	double n2=n2h_MD.nmax;
	int N=101;
	double dn=(n2-n1)/(N-1);
	double nw,hz_n,rh_n; 
	printf("# n(H2O), hz, rH\n");
	for(i=0;i<N;i++){
		nw=n1+dn*i;
		hz_n=n2h_MD.get_h(nw); // MD-swelling curve
		rh_n=r2h_XRD.get_rH(hz_n);
		printf("%lf, %lf, %lf\n",nw,hz_n,rh_n);
	};
*/
	// rH --> hz --> n(H2O)
	int i,N=101;
	double drH=1./(N-1);
	double hz_r,nw_r;
	printf("# rH, hz, n(H2O)\n");
	for(i=0;i<N;i++){
		hz_r=r2h_XRD.get_hz(i*drH);
		nw_r=n2h_MD.get_nH2O(hz_r);
//		printf("%lf, %lf, %lf\n",i*drH,hz_r,nw_r);
	};

	i=0;
	hz_r=r2h_XRD.get_hz(i*drH);
	nw_r=n2h_MD.get_nH2O(hz_r);
	double nmin=nw_r;
	i=N;
	hz_r=r2h_XRD.get_hz(i*drH);
	nw_r=n2h_MD.get_nH2O(hz_r);
	double nmax=nw_r;
//	printf("nmin=%lf, nmax=%lf\n",nmin,nmax);


	double dn=(nmax-nmin)/(N-1);
	double nw,hz_n,rh_n,mu; 
	printf("# n(H2O), hz, rH\n");
	double beta=1./18.;
	double Gh=0.0;
	for(i=N-1;i>=0;i--){
		nw=nmin+dn*i;
		hz_n=n2h_MD.get_h(nw); // MD-swelling curve
		rh_n=r2h_XRD.get_rH(hz_n);
		//mu=-1.+beta*log(rh_n+1.e-08);
		mu=beta*log(rh_n+1.e-08);
		printf("%lf, %lf, %lf, %lf, %lf\n",nw,hz_n,rh_n,mu,Gh);
		Gh+=(-mu);
	};

	return(0);
};
