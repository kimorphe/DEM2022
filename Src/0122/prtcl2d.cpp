#include "dem.h"
#include "math.h"

void PRTCL::init(){
	x[0]=0.0; x[1]=0.0;
	v[0]=0.0; v[1]=0.0;
	irev[0]=0; irev[1]=0;
	F[0]=0.0; F[1]=0.0;
	m=1.0;	// mass
	mobile=1;	// set to "mobile" mode
}; 

void PRTCL::setX(
	double x1, double x2,	// given particle coordinate
	REV rev,
	CNTRL prms		// DEM controle parameter incl. cell info			
){
	Vec2 xf;
	xf.set(x1-rev.Xa[0],x2-rev.Xa[1]);
	double xa=iprod(xf,rev.ah);
	double xb=iprod(xf,rev.bh);
	
	irev[0]=0; irev[1]=0;
	if(prms.iprd[0]==1) irev[0]=floor(xa);
	if(prms.iprd[1]==1) irev[1]=floor(xb);

	xa-=floor(xa);
	xb-=floor(xb);
	xf= vsum(rev.av.times(xa),rev.bv.times(xb));
	x[0]=xf.x[0]+rev.Xa[0];
	x[1]=xf.x[1]+rev.Xa[1];

};

void PRTCL::setX(double x1, double x2){

	x[0]=x1;
	x[1]=x2;

};

void PRTCL::incV(CNTRL prms){
	
	double edt2=prms.eta2;
	double dt=prms.dt;

	if( mobile != 1) return;
	v[0]=(v[0]*(m-edt2)+F[0]*prms.dt)/(m+edt2);
	v[1]=(v[1]*(m-edt2)+F[1]*prms.dt)/(m+edt2);
};

void PRTCL::incX(CNTRL prms, REV rev){

	int i;

	if( mobile != 1) return;
	x[0]+= (v[0]*prms.dt);
	x[1]+= (v[1]*prms.dt);

	double xi[2];
	Vec2 xf;
	xf.set(x[0]-rev.Xa[0],x[1]-rev.Xa[1]);
	xi[0]=iprod(xf,rev.ah);
	xi[1]=iprod(xf,rev.bh);

	for(i=0;i<2;i++){
	if(prms.iprd[i] == 1){
		//if(xi[i]>1.0){
		while(xi[i]>1.0){
			xi[i]-=1.0;
			irev[i]++;
		}
		//if(xi[i] < 0.0){
		while(xi[i] < 0.0){
			xi[i]+=1.0;
			irev[i]--;
		}
	}
	}
	xf=vsum(rev.av.times(xi[0]), rev.bv.times(xi[1]));
	x[0]=xf.x[0]+rev.Xa[0];
	x[1]=xf.x[1]+rev.Xa[1];
};

void PRTCL :: scale(int it,double dt, REV rev){

	int i;
	double Wt1[2], Wt2[2],xi[2];

	Vec2 xf;
	double xa,xb;

	xf.set(x[0]-rev.Xa[0],x[1]-rev.Xa[1]);
	xa=iprod(xf,rev.ah);
	xb=iprod(xf,rev.bh);
	
	//rev.update(it+1,dt);
	rev.update(it,dt);
	xf=vsum( rev.av.times(xa), rev.bv.times(xb));
	for(i=0;i<2;i++){
		if(rev.iprd[i] == 1) x[i]=xf.x[i]+rev.Xa[i];
	}
};
void PRTCL :: x2xh(REV rev){
//	int i;
	//double Wt1[2], Wt2[2],xi[2];
	Vec2 xf;
	//double xa,xb;
	xf.set(x[0]-rev.Xa[0],x[1]-rev.Xa[1]);
	//xa=iprod(xf,rev.ah);
	//xb=iprod(xf,rev.bh);
	x[0]=iprod(xf,rev.ah);
	x[1]=iprod(xf,rev.bh);
	
	//rev.update(it+1,dt);
	/*
	rev.update(it,dt);
	xf=vsum( rev.av.times(xa), rev.bv.times(xb));
	for(i=0;i<2;i++){
		if(rev.iprd[i] == 1) x[i]=xf.x[i]+rev.Xa[i];
	}
	*/
};
void PRTCL :: xh2x(REV rev){
	Vec2 xf;
	//xf=vsum( rev.av.times(xa), rev.bv.times(xb));
	xf=vsum( rev.av.times(x[0]), rev.bv.times(x[1]));
	for(int i=0;i<2;i++){
		if(rev.iprd[i] == 1) x[i]=xf.x[i]+rev.Xa[i];
	}
};

double PRTCL :: KE(double m0){	// Kinetic energy  
	double dKE;

	dKE=0.5*(v[0]*v[0]+v[1]*v[1])*m*m0; // [Nm = kg (nm/ps)^2]	
	dKE*=1.e06;	// [kg(m/s)^2]

	return(dKE);
};
double PRTCL::get_sig_iso(double xf[2], int *sgn_xf, double dsig[2]){

	double rf[2];
	rf[0]=xf[0]-x[0];
	rf[1]=xf[1]-x[1];
	//doubler rr=sqrt(rf[0]*rf[0]+rf[1]*rf[1]);
	//rf[0]/=rr;
	//rf[1]/=rr;

	double mn=rf[0]*nx[0]+rf[1]*nx[1];
	double sigf;
	if(mn > 0.0){
		*sgn_xf=0;
		sigf=sigs[0]+dsig[0];
	}else{
		*sgn_xf=1;
		sigf=sigs[1]+dsig[1];
	};
	return(sigf);
};
double PRTCL::get_sig(double xf[2], int *sgn_xf, double dsig[2]){

	Vec2 yf1,yf0;	// use Vec2 instances

	yf0.set(x);	// register particle coordinate
	yf1.set(xf);	// register field point 

	yf0=vdiff(yf1,yf0); // get relative position vector
	double rr=yf0.len(); // vector length
	yf0=yf0.div(rr); // make vector length unity 

	double mt,mn;

	yf1.x[0]= nx[1];	// tangential vector
	yf1.x[1]=-nx[0];	// tangential vector
	mt=iprod(yf0,yf1); // tangential compoonent

	yf1.x[0]= nx[0];	// tangential vector
	yf1.x[1]= nx[1];	// tangential vector
	mn=iprod(yf0,yf1); // normal component


	//double st=mt/sig;
	double sig_zero=1.35;
	double st=mt/sig_zero;
	double sn=mn/(sigs[0]+dsig[0]);
	*sgn_xf=0;	// xf is located above this molecule 
	if(sn < 0.0){
	       sn=mn/(sigs[1]+dsig[1]);
	       *sgn_xf=1;	// xf is located below this molecule
	}
	return(1./sqrt(st*st+sn*sn));

};
