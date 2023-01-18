#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dem.h"

using namespace std;

void show_msg(char *fname){
	printf("Can't find '%s'\n",fname);
	printf(" --> process terminated.");
	exit(-1);
}

void CNTRL :: load(char *fname){
	int i;
	char cbff[128];
	FILE *fp=fopen(fname,"r");
	double Xb[2],ma;
	double PI=4.0*atan(1.0);
	if(fp==NULL) show_msg(fname);
	
	fgets(cbff,128,fp);
	fgets(Dir,128,fp);


	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&rstat);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fnrst);
	puts(fnrst);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&Nout);

	fgets(cbff,128,fp);
	fscanf(fp,"%le %le\n",&m0,&Eps0);// [kg],[kg m^2/sec^2=Nm] 
	printf("%le %le\n",m0,Eps0);// [kg],[kg m^2/sec^2=Nm] 

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %d\n",&dt,&Nt); // [ps]
	itime=0;
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&eta);	// dimensionless parameter
	eta2=eta*eta;

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",&K1, &K2); // [kg/s^2]
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&sig);	// [nm]

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&mvw);	// 0: water fixed, 1: move 
	printf("mvw=%d\n",mvw);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",&vmin,&vmax);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %lf  %lf %lf\n",&T_cntrl, &Tset, &time_T_start,&time_T_end);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xa,Xa+1); // [nm]
	fscanf(fp,"%lf %lf\n",Xb,Xb+1); // [nm]
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",iprd,iprd+1);
		prd=iprd[0]+iprd[1];

		if(prd > 0){
			FILE *ftmp=fopen("bbox.dat","r");
			strcpy(cbff,"bbox.dat");
			if(ftmp == NULL) show_msg(cbff);
			fscanf(ftmp,"%lf %lf\n",Xa,Xa+1);
			fgets(cbff,128,ftmp);
			fscanf(ftmp,"%lf %lf\n",Xb,Xb+1);
			fclose(ftmp);
		}


		Wd0[0]=Xb[0]-Xa[0];
		Wd0[1]=Xb[1]-Xa[1];
		Wd0[2]=0.0;
		Wd0[3]=0.0;
		Xc[0]=0.5*(Xa[0]+Xb[0]);
		Xc[1]=0.5*(Xa[1]+Xb[1]);
		for(i=0;i<4;i++) Wd[i]=Wd0[i];

		for(i=0;i<3;i++) Wd[i]=Wd0[i];
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&cntrl); // 0: strain, 1: stress driven, -1: cell contst read from a file  
	printf("cntnrl=%d\n",cntrl);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf %lf\n",exx, exx+1, txx, txx+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf %lf\n",eyy, eyy+1, tyy, tyy+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf %lf\n",exy, exy+1, txy, txy+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf %lf\n",eyx, eyx+1, tyx, tyx+1);
	for(i=0;i<2;i++){
		exy[i]=exy[i]/180.*PI;
		eyx[i]=eyx[i]/180.*PI;
	}

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf\n",&sxx,&syy,&sxy);
	printf("(sxx,syy,sxy)=(%lf %lf %lf)\n",sxx,syy,sxy);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nave);
	printf("nave=%d\n",nave);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&mu);
	printf("mu=%lf\n",mu);
	fclose(fp);
	

	Eps=Eps0*1.e-06;	// [kg (nm/ps)^2]
	K1*=1.e-24;	// [kg/ ps^2]
	K2*=1.e-24;	// [kg/ ps^2]

	Eps=Eps*24.0/(m0*sig);
	K1/=m0;
	K2/=m0;

	ma=24.*Eps/sig;	// unit force [kg nm/ps^2]

	printf("m0=%le Eps=%le\n",m0,Eps);
	printf(" K1=%le K2=%le\n",K1,K2);

//		Subcell size
	h0[0]=sig*3.5;	// [nm]
	h0[1]=sig*3.5;	// [nm]

	
};
void CNTRL :: print(){
		
	printf("rstat=%d\n",rstat);
	printf("dt=%lf,Nt=%d\n",dt,Nt);
	printf("eta=%lf (viscous constant)\n",eta);
	printf("sig=%lf (characteristic length)\n",sig);
	printf("cntrl=%d\n",cntrl);
	printf("K1=%lf K2=%lf (spring constants)\n",K1,K2);
};

void WALL :: setup(){

	int i;
	char cbff[128],fname[128];
	strcpy(fname,"wall.dat");
	FILE *fp=fopen(fname,"r");

	if(fp == NULL) show_msg(fname);
	fgets(cbff,128,fp);
	for(i=0;i<4;i++){
		fgets(cbff,128,fp);
		fscanf(fp,"%lf %lf %lf %lf\n",ux+i,uy+i,Tu+i, brst+i);
	}
	fgets(cbff,128,fp);
	nwall=0;
	for(i=0;i<4;i++){
		fscanf(fp,"%lf %lf %d\n",x0+i,y0+i,iwl+i);
		printf("%lf %lf %d\n",x0[i],y0[i],iwl[i]);
		nwall+=iwl[i];
		x[i]=x0[i];
		y[i]=y0[i];
	}
	printf("nwall=%d\n",nwall);
}

void WALL :: disp(CNTRL prms, int it){

	int jw;
	double dt=prms.dt;
	double PI=4.0*atan(1.0);

	for(jw=0;jw<4;jw++){
		if(iwl[jw] ==0) continue;
		if(it*dt <= Tu[jw]*brst[jw]){
			x[jw]=x0[jw]+ux[jw]*sin(2.*PI*(it*dt/Tu[jw]));
			y[jw]=y0[jw]+uy[jw]*sin(2.*PI*(it*dt/Tu[jw]));
		}
	}
};

void REV::set_lattice_vec(){

	double cosab;

	av.set(cos(Wd[2]),sin(Wd[2]));
	av=av.times(Wd[0]);

	bv.set(sin(Wd[3]),cos(Wd[3]));
	bv=bv.times(Wd[1]);

	cosab=cos(Wd[2]+Wd[3]);
	ah.set(cos(Wd[3]),-sin(Wd[3]));
	bh.set(-sin(Wd[2]),cos(Wd[2]));

	ah=ah.div(Wd[0]*cosab);
	bh=bh.div(Wd[1]*cosab);
};

void REV :: setup(CNTRL prms, WALL wll){

	int i;
	double Xb[2];
	FILE *fp;
	char cbff[128];

	for(i=0;i<2;i++){
		Xa[i]=prms.Xa[i];
		Wd[i]=prms.Wd[i];
		Wd0[i]=prms.Wd0[i];
		Xb[i]=Xa[i]+Wd0[i];
	}
	Wd[2]=prms.Wd[2]; // angle(rad) from x-axis
	Wd[3]=prms.Wd[3]; // angle(rad) from y-axis
	Wd0[2]=prms.Wd0[2]; // angle(rad) from x-axis
	Wd0[3]=prms.Wd0[3]; // angle(rad) from y-axis

	sxx=prms.sxx;
	syy=prms.syy;
	sxy=prms.sxy;
	cntrl=prms.cntrl;

	if(cntrl == -1){
		char fname[]="cell_const.dat";
		fp=fopen(fname,"r");
		fgets(cbff,128,fp);
		fscanf(fp,"%d\n",&nwt);
		printf("nwt=%d\n",nwt);
			w0=(double *)malloc(sizeof(double)*nwt);
			w1=(double *)malloc(sizeof(double)*nwt);
			w2=(double *)malloc(sizeof(double)*nwt);
			w3=(double *)malloc(sizeof(double)*nwt);
		fgets(cbff,128,fp);
		for(i=0;i<nwt;i++){
			fscanf(fp,"%lf,%lf,%lf,%lf\n",w0+i,w1+i,w2+i,w3+i);
		}
		printf("%lf,%lf,%lf,%lf\n",w0[i-1],w1[i-1],w2[i-1],w3[i-1]);
		fclose(fp);		

		Wd[0]=w0[0];
		Wd[1]=w1[0];
		Wd[2]=w2[0];
		Wd[3]=w3[0];
		for(i=0;i<4;i++) Wd0[i]=Wd[i];
	}

	nave=prms.nave;
	bfxx=(double *)malloc(sizeof(double)*nave);
	bfyy=(double *)malloc(sizeof(double)*nave);
	bfxy=(double *)malloc(sizeof(double)*nave);
	bfT= (double *)malloc(sizeof(double)*nave);
	bfU= (double *)malloc(sizeof(double)*nave);
	for(i=0;i<nave;i++){
		bfxx[i]=0.0;
		bfyy[i]=0.0;
		bfxy[i]=0.0;
		bfT[i]=0.0;
		bfU[i]=0.0;
	};
	ibf=0;

	kxx=0.0001;	// normal strain increment/GPa
	kyy=0.0001;	// normal strain increment/GPa
	kxy=0.005;	// shear strain increment[rad]/GPa;
	//kxx=0.001;	// normal strain increment/GPa
	//kyy=0.001;	// normal strain increment/GPa
	//kxy=0.001;	// shear strain increment[rad]/GPa;

	set_lattice_vec();
	if(prms.iprd[0]==0){
		if(wll.iwl[1]==1) Xb[0]=wll.x0[1]; 
		if(wll.iwl[3]==1) Xa[0]=wll.x0[3]; 
	}
	if(prms.iprd[1]==0){
		if(wll.iwl[0]==1) Xa[1]=wll.y0[0]; 
		if(wll.iwl[2]==1) Xb[1]=wll.y0[2]; 
	}
	
	for(i=0;i<2;i++){
		Xa0[i]=Xa[i];
		txx[i]=prms.txx[i]; txy[i]=prms.txy[i];
		tyx[i]=prms.tyx[i]; tyy[i]=prms.tyy[i];

		exx[i]=prms.exx[i]; exy[i]=prms.exy[i];
		eyx[i]=prms.eyx[i]; eyy[i]=prms.eyy[i];

		iprd[i]=prms.iprd[i];
	}


	
	h0[0]=prms.h0[0];
	h0[1]=prms.h0[1];

	Nh[0]=ceil(Wd[0]/h0[0]);
	Nh[1]=ceil(Wd[1]/h0[1]);
	hd[0]=Wd[0]/Nh[0];
	hd[1]=Wd[1]/Nh[1];
	nsub=Nh[0]*Nh[1];

	sxxb_pre=0.0; syyb_pre=0.0; sxyb_pre=0.0;
	sxxb=0.0; syyb=0.0; sxyb=0.0;

};

void REV :: smooth(double S[2][2], double Tk, double UE){

	int i,nbf;

	bfxx[ibf%nave]=S[0][0];
	bfyy[ibf%nave]=S[1][1];
	bfxy[ibf%nave]=S[1][0];
	bfT[ibf%nave]=Tk;
	bfU[ibf%nave]=UE;

	ibf++;
	if((ibf%nave)==2){
		sxxb_pre=sxxb; 
		syyb_pre=syyb; 
		sxyb_pre=sxyb; // stresses at previous step
		Tb_pre=Tb;
		Ub_pre=Ub;
	}

	sxxb=0.0,syyb=0.0,sxyb=0.0;
	Tb=0.0;
	Ub=0.0;

	nbf=ibf;
	if(nbf > nave) nbf=nave;	
	for(i=0;i<nbf;i++){
		sxxb+=bfxx[i];
		syyb+=bfyy[i];
		sxyb+=bfxy[i];
		Tb+=bfT[i];
		Ub+=bfU[i];
	};

	sxxb/=nbf;
	syyb/=nbf;
	sxyb/=nbf;
	Tb/=nbf;
	Ub/=nbf;
	
};

void REV :: print(){
	printf("---- REV Parameters ----\n");
	printf("Xa=(%lf, %lf)\n",Xa[0],Xa[1]);
	printf("Wd =(%lf, %lf, %lf, %lf)\n",Wd[0],Wd[1],Wd[2],Wd[3]);
	printf("Wd0=(%lf, %lf, %lf, %lf)\n",Wd0[0],Wd0[1],Wd0[2],Wd0[3]);
	printf("iprd[0]=%d, iprd[1]=%d\n",iprd[0],iprd[1]);
	printf("exx=[%lf, %lf]  txx=[%lf, %lf]\n",exx[0],exx[1],txx[0],txx[1]);
	printf("eyy=[%lf, %lf]  tyy=[%lf, %lf]\n",eyy[0],eyy[1],tyy[0],tyy[1]);
	printf("exy=[%lf, %lf]  txy=[%lf, %lf]\n",exy[0],exy[1],txy[0],txy[1]);
	printf("eyx=[%lf, %lf]  tyx=[%lf, %lf]\n",eyx[0],eyx[1],tyx[0],tyx[1]);
	printf("av="); av.print();
	printf("bv="); bv.print();
	printf("ah="); ah.print();
	printf("bh="); bh.print();
	printf("------------------------\n");
};

void REV :: update(int it, double dt){

	double tt=it*dt;
	double epst,gmmt;
	double PI=4.0*atan(1.0);

	if(cntrl==-1){
		Wd[0]=w0[it];
		Wd[1]=w1[it];
		Wd[2]=w2[it];
		Wd[3]=w3[it];
	}else if(cntrl==0){	// strain given condition
		if(iprd[0] ==1){
			epst=exx[0];
			if(tt >= txx[0]) epst+=(tt-txx[0])/(txx[1]-txx[0])*(exx[1]-exx[0]);
			if(tt >= txx[1]) epst=exx[1];

			gmmt=exy[0];
			if(tt >= txy[0]) gmmt+=(tt-txy[0])/(txy[1]-txy[0])*(exy[1]-exy[0]);
			if(tt >= txy[1]) gmmt=exy[1];

			Wd[0]=Wd0[0]*(1.0+epst);
			Wd[2]=Wd0[2]+gmmt;
		};
		if(iprd[1]==1){
			epst=eyy[0];
			if(tt >= tyy[0]) epst+=(tt-tyy[0])/(tyy[1]-tyy[0])*(eyy[1]-eyy[0]);
			if(tt >= tyy[1]) epst=eyy[1];
			gmmt=eyx[0];
			if(tt >= tyx[0]) gmmt+=(tt-tyx[0])/(tyx[1]-tyx[0])*(eyx[1]-eyx[0]);
			if(tt >= tyx[1]) gmmt=eyx[1];

			Wd[1]=Wd0[1]*(1.0+epst);
			Wd[3]=Wd0[3]+gmmt;
		};
	}else if(cntrl==1){	// stress given condition
		//if( (ibf%nave)==2 ){
		if( ibf>=2 ){
		double dWxx,dWyy,dWxy;
		double sgn;

		sgn=1.0;
		if((sxxb-sxx)<0.0) sgn=-1.0;
		dWxx=kxx*(sxxb-sxx)*Wd[0];
		//dWxx=kxx*sgn*Wd[0];

		sgn=1.0;
		if((syyb-syy)<0.0) sgn=-1.0;
		dWyy=kyy*(syyb-syy)*Wd[1];
		//dWyy=kyy*sgn*Wd[1];

		sgn=0.0;
		dWxy=kxy*(sxyb-sxy);
		//dWxy=kxy*sgn;

		Wd[0]+=dWxx;
		Wd[1]+=dWyy;
		Wd[2]+=dWxy;
		Wd[3]+=dWxy;
		}
	
	};

	Nh[0]=ceil(Wd[0]/h0[0]);
	Nh[1]=ceil(Wd[1]/h0[1]);
	hd[0]=Wd[0]/Nh[0];
	hd[1]=Wd[1]/Nh[1];
	nsub=Nh[0]*Nh[1];

	Vol=Wd[0]*Wd[1]*sin(PI*0.5-Wd[2]-Wd[3]);

	set_lattice_vec();
};

void REV :: update(int it, double dt, WALL wll){	// <-- Currently in use (01/17/2017)

	int i;
	double tt=it*dt;	// current time
	double epst,gmmt;	// ramp function values
	double Xb[2];
	double PI=4.0*atan(1.0);

	if(cntrl==-1){
		Wd[0]=w0[it];
		Wd[1]=w1[it];
		Wd[2]=w2[it];
		Wd[3]=w3[it];
	}else if(cntrl==0){
		if(iprd[0] ==1){	// Periodic in X(a)-direction
			epst=exx[0];	// for normal strain
			if(tt >= txx[0]) epst+=(tt-txx[0])/(txx[1]-txx[0])*(exx[1]-exx[0]);
			if(tt >= txx[1]) epst=exx[1];

			gmmt=exy[0];	// for shear strain
			if(tt >= txy[0]) gmmt+=(tt-txy[0])/(txy[1]-txy[0])*(exy[1]-exy[0]);
			if(tt >= txy[1]) gmmt=exy[1];

			Wd[0]=Wd0[0]*(1.0+epst);
			Wd[2]=Wd0[2]+gmmt;

		}else{	// Non-Periodic in X-direction
			if(wll.iwl[1]==1) Xb[0]=wll.x[1]; 
			if(wll.iwl[3]==1) Xa[0]=wll.x[3]; 
		}

		if(iprd[1]==1){	// Periodic in Y(b)-direction
			epst=eyy[0];	// for normal strain
			if(tt >= tyy[0]) epst+=(tt-tyy[0])/(tyy[1]-tyy[0])*(eyy[1]-eyy[0]);
			if(tt >= tyy[1]) epst=eyy[1];

			gmmt=eyx[0];	// for shear strain
			if(tt >= tyx[0]) gmmt+=(tt-tyx[0])/(tyx[1]-tyx[0])*(eyx[1]-eyx[0]);
			if(tt >= tyx[1]) gmmt=eyx[1];

			Wd[1]=Wd0[1]*(1.0+epst);
			Wd[3]=Wd0[3]+gmmt;

		}else{
			if(wll.iwl[0]==1) Xa[1]=wll.y[0]; 
			if(wll.iwl[2]==1) Xb[1]=wll.y[2]; 
		};
	}else if(cntrl==1){
		//if( (ibf%nave)==2 ){
		if( ibf>=2 ){

		dWxx=kxx*(sxxb-sxx)*Wd[0];
		dWyy=kyy*(syyb-syy)*Wd[1];
		//dWxy=kxy*(sxyb-sxy);
		dWxy=0.0;

		/*
		double sgn=1.0;
		if((sxxb-sxx)<0.0) sgn=-1.0;
		dWxx=kxx*sgn*Wd[0];
		sgn=1.0;
		if((syyb-syy)<0.0) sgn=-1.0;
		dWyy=kyy*sgn*Wd[1];
		dWxy=0.0;
		*/
		Wd[0]+=dWxx;
		Wd[1]+=dWyy;
		Wd[2]+=dWxy;
		Wd[3]+=dWxy;
		}
	};

	Nh[0]=ceil(Wd[0]/h0[0]);
	Nh[1]=ceil(Wd[1]/h0[1]);
	hd[0]=Wd[0]/Nh[0];
	hd[1]=Wd[1]/Nh[1];
	nsub=Nh[0]*Nh[1];

	set_lattice_vec();
	Xb[0]=Xa[0]+av.x[0]+bv.x[0];
	Xb[1]=Xa[1]+av.x[1]+bv.x[1];

	Vol=Wd[0]*Wd[1]*sin(PI*0.5-Wd[2]-Wd[3]);
};
Vec2 REV::unfold(double *xcod, int *ofst){
	Vec2 xf;
	double xa,xb;

	xf.set(xcod[0]-Xa[0],xcod[1]-Xa[1]);
	xa=iprod(xf,ah);
	xb=iprod(xf,bh);

	xf=vsum(av.times(xa+ofst[0]), bv.times(xb+ofst[1]));
	xf.x[0]+=Xa[0];
	xf.x[1]+=Xa[1];
	return(xf);

};

void SUBCELL:: setup(REV rev, double sig){

	double sigd=sig*0.5;
	np_max=ceil(rev.hd[0]/sigd)*ceil(rev.hd[1]/sigd);
	indx[0]=0;
	indx[1]=0;
	list=(int *)malloc(sizeof(int)*np_max);
	np=0;
};
SUBCELL::~SUBCELL(){
	if(list != NULL){
//		delete [] list;
//		free(list);
	}
};
