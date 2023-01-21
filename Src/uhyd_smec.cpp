#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"dem.h"

/*
double Uhyd_sig(double U0, double sig){
	double sig0=0.9;	//[nm]
	double sigb=1.1;	//[nm] break point
	double alph=log(2)/(sigb-sig0);
	return( -U0*(1.0-exp(-alph*(sig-sig0))) );
};
class Crv{
	public:
		int np;
		void setup(int n);
		void set_xlim(double x1, double x2);
		double *x,*y;
		double xs,xe,dx;
		double xmin,xmax;
		void set_stair(double *Xk, int NXk);
		void write(char fname[128]);
		void saw_tooth(double *Xk, int NXk,double gmm);
		void smooth(int nsmp);
		void trend(int type,double amp_wv);
		double eval(double xval);
	private:
};
class CLAY{
	public:
		void load(char fname[128]);
		Crv hz; // basal spacing [nm]
		Crv rH; // relative humidity
		Crv mu_var; // chemical potential (nolineary varying part) [J/particle/face]
		Crv G_var; // hydration free energy 
		double mu_sat; // chemical potential at saturated R.H.
		int ndat;
	private:
	protected:
};
*/
void CLAY::load(char fname[128]){
	FILE *fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname);
	
	fscanf(fp,"%d\n",&ndat);
	hz.setup(ndat);
	rH.setup(ndat);
	mu_var.setup(ndat);
	G_var.setup(ndat);

	char cbff[256];
	fgets(cbff,256,fp);
	double Na=6.022e+23;

	int i;
	double nw, hn, rn, mvn,gvn, gsn; 
	for(i=0;i<ndat;i++){
		fscanf(fp,"%lf, %lf, %lf, %lf, %lf, %lf, %lf\n",&nw, &hn,&rn,&mvn,&mu_sat,&gvn,&gsn);
		hz.x[i]=nw; 
		rH.x[i]=nw; 
		mu_var.x[i]=nw;
		G_var.x[i]=nw;

		hz.y[i]=hn;	//[nm]
		rH.y[i]=rn;
		mu_var.y[i]=mvn*1.e03/Na; // [J]
		G_var.y[i]=gvn*1.e03/Na;  // [J]
	};
	mu_sat=mu_sat*1.e03/Na; //[J]


	double n1,n2;
	n1=hz.x[0];
	n2=hz.x[ndat-1];
	hz.set_xlim(n1,n2);
	rH.set_xlim(n1,n2);
	mu_var.set_xlim(n1,n2);
	G_var.set_xlim(n1,n2);
	nw_min=n1;
	nw_max=n2;

	fclose(fp);
};
void CLAY::change_unit(double unit){
	mu_sat/=unit;
	for(int i=0;i<ndat;i++){
		mu_var.y[i]/=unit;
		G_var.y[i]/=unit;
	};
};
void CLAY::print(){
	printf("-------------------------------------------------\n");
	printf("mu_sat=%le \n",mu_sat);
	printf("# n(H2O), hz[nm], R.H, mu_var, G_var\n");
	for(int i=0;i<ndat;i++){
		printf("%le, %le, %le, %le, %le\n",hz.x[i],hz.y[i],rH.y[i],mu_var.y[i],G_var.y[i]);
	};
};

//--------------------------------------------
void Crv::load(char fname[128]){
	FILE *fp=fopen(fname,"r");
	char cbff[128];
	if(fp==NULL){
		printf("File %s cannot found!",fname);
		printf(" --> abort process\n");
		exit(-1);
	};
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&np);
	printf("ndat=%d\n",np);
	Crv::setup(np);

	double tmp;
	for(int i=0;i<np;i++){
		fscanf(fp,"%lf, %lf, %lf\n",x+i,y+i,&tmp);
	};

	fclose(fp);
};
void Crv::setup(int n){
	np=n;

	x=(double *)malloc(sizeof(double)*np);
	y=(double *)malloc(sizeof(double)*np);
	for(int i=0;i<np;i++){
		x[i]=0.0;
		y[i]=0.0;
	};
};

void Crv::set_xlim(double x1,double x2){
	xmin=x1;
	xmax=x2;
	dx=(x2-x1)/(np-1);

	for(int i=0;i<np;i++) x[i]=xmin+dx*i;
};
void Crv::set_stair(double *Xk, int NXk){
	int i=0;
	double fac;
	for(int k=0;k<NXk;k++){
		if(i>=np) break;
		while(x[i]<Xk[k]){
			y[i]=k;
			i++;
			if(i>=np) break;
		};
	}
};
void Crv::saw_tooth(double *Xk, int NXk,double gmm){
	int i=0;
	int k;
	double Wk,xb;
	double fac;
	for(i=0;i<np;i++){
		k=int(y[i]);
		fac=pow(gmm,k-1);
		Wk=Xk[k]-Xk[k-1];
		xb=(x[i]-Xk[k-1])/Wk;
		y[i]=xb*fac;
	};
};
void Crv::write(char fname[128]){
	FILE *fp=fopen(fname,"w");
	int i;
	for(i=0;i<np;i++) fprintf(fp,"%lf, %lf\n",x[i],y[i]);
	fclose(fp);
};
void Crv::smooth(int nsmp){
	double sum;
	int N=nsmp/2;
	int i,j,k;
	double *tmp=(double *)malloc(sizeof(double)*np);
	for(i=0;i<np;i++) tmp[i]=y[i];

	for(i=0;i<np;i++){
		sum=0.;
		for(k=i-N;k<=i+N;k++){
			j=k;
			if(j<0) j+=np;
			if(j>=np) j%np;
			sum+=tmp[j];
		};
		y[i]=sum/(2*N+1);
	};
};
void Crv::trend(int type,double amp_wv){
	if(type==0){
		double alph=2.0;
		for(int i=0;i<np;i++){
			y[i]=5.*exp(-alph*x[i])-y[i];
		};
	}
	if(type==1){
		double U0=1.0;
		for(int i=0;i<np;i++){
			y[i]=Uhyd_sig(U0,x[i])-y[i]*amp_wv;
		};
	};

};
double Crv::eval(double xval){
	double yval;

	if(xval<=x[0]) return(y[0]);
	if(xval>=x[np-1]) return(y[np-1]);
	int indx=int((xval-x[0])/dx);
	double xi=(xval-x[0])/dx-indx;

	double y1,y2;

	y1=y[indx];
	y2=y[indx+1];
	yval=y1*(1.0-xi)+y2*xi;
	return(yval);
};
double Crv::y2x(double yval){
	if(yval<=y[0]) return(x[0]);
	if(yval>=y[np-1]) return(x[np-1]);

	int i1=0;
	while(y[i1+1]<yval) i1++;
	double x1=x[i1];
	double x2=x[i1+1];
	double xi=(yval-y[i1])/(y[i1+1]-y[i1]);
	double eta=1.-xi;
	return(x1*eta+x2*xi);

};
//-------------------------------------

