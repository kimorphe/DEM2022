#include<stdio.h>
#include<stdlib.h>

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
	return(ndat);
};
//-----------------------------------------------------------
int main(){
	char fname[128]="NaMt.dat";

	int ndat_h=51;

	N2H swelling_curve;	// MD-generated swelling curve
	swelling_curve.load(fname);
	swelling_curve.setup_H2N(ndat_h);
	swelling_curve.fwrite(); // --> "n2h.dat" and "h2n.dat"

/*
	int i,N=46;
	double dn=0.2;
	for(i=0;i<N;i++){
//		printf("%lf, %lf\n",i*dn,hnc.get_h(i*dn));
	}
	printf("\n");

	int M=46;
	double h1=hnc.hmin;
	double h2=hnc.hmax;
	double dh=(h2-h1+1.)/(M-1);
	for(i=0;i<M;i++){
		printf("%lf, %lf\n",i*dh+h1,hnc.get_nH2O(i*dh+h1));
	};
*/

	return(0);
};
