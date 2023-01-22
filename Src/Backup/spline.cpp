#include <stdio.h>
#include <stdlib.h>
#include "dem.h"
///////////////////////// Member Functions //////////////////////////////
Curve2D::~Curve2D(){
	/*
	free(x); free(y); 
	free(bx); free(cx); free(dx);
	free(by); free(cy); free(dy);
	*/
};
void Curve2D:: init(int N){
	np=N;
	x=(double *)malloc(sizeof(double)*np);
	y=(double *)malloc(sizeof(double)*np);

	int npd;
	
	npd=np-1;
	bx=(double *)malloc(sizeof(double)*npd);
	cx=(double *)malloc(sizeof(double)*npd);
	dx=(double *)malloc(sizeof(double)*npd);

	by=(double *)malloc(sizeof(double)*npd);
	cy=(double *)malloc(sizeof(double)*npd);
	dy=(double *)malloc(sizeof(double)*npd);

	for(int i=0;i<np;i++){
		x[i]=0.0;
		y[i]=0.0;
	};

	for(int i=0;i<npd;i++){
		bx[i]=0.0;
		cx[i]=0.0;
		dx[i]=0.0;

		by[i]=0.0;
		cy[i]=0.0;
		dy[i]=0.0;
	}
};
void Curve2D:: smooth(int nsmp2){
	int i,j,indx;
	double *xb,*yb;
	int nsmp=nsmp2*2+1;

	xb=(double *)malloc(sizeof(double)*np);
	yb=(double *)malloc(sizeof(double)*np);
	for(i=0;i<np;i++){
		xb[i]=0.0;
		yb[i]=0.0;
	for(j=-nsmp2; j<=nsmp2; j++){
		indx=j;
		if(indx<0) indx=0;
		if(indx>=np) indx=np-1;
		xb[i]+=x[indx];
		yb[i]+=y[indx];
	}
	xb[i]/=nsmp;
	yb[i]/=nsmp;
	}
	free(xb);
	free(yb);
};
void Curve2D::print(){
	puts("-----------------------------");
	for(int i=0;i<np;i++){
		printf("x,y[%d]=%lf %lf\n",i,x[i],y[i]);
	};
	puts("-----------------------------");
}
void Curve2D::print(char *fname){
	FILE *fp=fopen(fname,"a");
	for(int i=0;i<np;i++) fprintf(fp,"%lf %lf\n",x[i],y[i]);
	fprintf(fp,"\n");
	fclose(fp);
};

void Curve2D::spline(){

	int i;
	MatTriDiag B;
	B.init(np-2);
	double h0,h1;
	double x0,x1,x2;
	double y0,y1,y2;
	double u0,u1;
	double *px=(double *)malloc(sizeof(double)*(np-2));
	double *py=(double *)malloc(sizeof(double)*(np-2));

	h0=1.0;	//h0=i+1-i;
	h1=1.0; //h1=i+2-(i+1);
	for(i=0;i<np-2;i++){
		x0=x[i];
		x1=x[i+1];
		x2=x[i+2];
		y0=y[i];
		y1=y[i+1];
		y2=y[i+2];
		px[i]=(x0/h0-(1./h0+1./h1)*x1+x2/h1);
		px[i]*=6.0;
		py[i]=(y0/h0-(1./h0+1./h1)*y1+y2/h1);
		py[i]*=6.0;

		B.A0[i]=h0;
		B.A1[i]=2.*(h0+h1);
		B.A2[i]=h1;
	}

	B.LinSolve(px,py);
	for(i=0;i<np-1;i++){
		u0=0.0;	
		u1=0.0;
		if(i>0) u0=px[i-1];
		if(i<np-2) u1=px[i];
		bx[i]=x[i+1]-x[i]-(2.*u0+u1)/6.0;
		cx[i]=u0*0.5;
		dx[i]=(u1-u0)/6.0;

		u0=0.0;	
		u1=0.0;
		if(i>0) u0=py[i-1];
		if(i<np-2) u1=py[i];
		by[i]=y[i+1]-y[i]-(2.*u0+u1)/6.0;
		cy[i]=u0*0.5;
		dy[i]=(u1-u0)/6.0;
	}

	free(px);
	free(py);
	free(B.A0);
	free(B.A1);
	free(B.A2);
	free(B.y);
}
double Curve2D::intplx(double s){
	int i=int(s);
	double xi=s-i;
	if(i>=np-1){
	      	i=np-2;
		xi=1.0;
	}

	return(x[i]+xi*(bx[i]+xi*(cx[i]+xi*dx[i])));
}
double Curve2D::intply(double s){
	//int i=int(s);
	//if(i>=np-1) i=np-2;
	//double xi=s-i;
	int i=int(s);
	double xi=s-i;
	if(i>=np-1){
	      	i=np-2;
		xi=1.0;
	}
	return(y[i]+xi*(by[i]+xi*(cy[i]+xi*dy[i])));
}

double Curve2D::dxds(double s){
	int i=int(s);
	double xi=s-i;
	if(i>=np-1){
	      	i=np-2;
		xi=1.0;
	}
	//int i=int(s);
	//if(i>=np-1) i=np-2;
	//double xi=s-i;
	return(bx[i]+xi*(2.*cx[i]+3.*xi*dx[i]));
};
double Curve2D::dyds(double s){
	int i=int(s);
	double xi=s-i;
	if(i>=np-1){
	      	i=np-2;
		xi=1.0;
	}
	//int i=int(s);
	//if(i>=np-1) i=np-2;
	//double xi=s-i;
	return(by[i]+xi*(2.*cy[i]+3.*xi*dy[i]));
};

void MatTriDiag::init(int N){
	n=N;
	A0=(double*)malloc(sizeof(double)*n);
	A1=(double*)malloc(sizeof(double)*n);
	A2=(double*)malloc(sizeof(double)*n);
	y=NULL;

	int i;
	for(i=0;i<n;i++){
		A0[i]=0.0;
		A1[i]=0.0;
		A2[i]=0.0;
	}
}

void MatTriDiag::print(){

	puts("---- MatTriDiag -----");
	for(int i=0;i<n;i++){
		printf("A0,A1,A2[%d]=%lf %lf %lf\n",i,A0[i],A1[i],A2[i]);
	}
};
void MatTriDiag::LinSolve(double *b){

	int l;
	double p,q;

	if(y==NULL) y=(double *)malloc(sizeof(double)*n);

	// Forward Elimination
	for(l=0;l<n-1; l++){	// row 
		p=A1[l]; //pivot
    		b[l]/=p;
		A1[l]/=p; 
		      A2[l]/=p;
		      q=A0[l+1];
		      A0[l+1]-=A1[l]*q; A1[l+1]-=A2[l]*q;
		      b[l+1]-=b[l]*q;
	} //row

	l=n-1;
	p=A1[l]; //pivot
    	b[l]/=p;
	A1[l]/=p; 
	y[n-1]=b[n-1];
	for(l=n-2;l>=0;l--){
		y[l]=b[l];
		y[l]-=A2[l]*y[l+1];
	}

};
void MatTriDiag::LinSolve(double *b1, double *b2){

	int l;
	double p,q;

	double *sol1=(double *)malloc(sizeof(double)*n);
	double *sol2=(double *)malloc(sizeof(double)*n);

	// Forward Elimination
	for(l=0;l<n-1; l++){	// row 
		p=A1[l]; //pivot
    		b1[l]/=p;
    		b2[l]/=p;
		A1[l]/=p; 
		      A2[l]/=p;
		      q=A0[l+1];
		      A0[l+1]-=A1[l]*q; A1[l+1]-=A2[l]*q;
		      b1[l+1]-=b1[l]*q;
		      b2[l+1]-=b2[l]*q;
	} //row

	l=n-1;
	p=A1[l]; //pivot
    	b1[l]/=p; 
	b2[l]/=p;
	A1[l]/=p; 
	sol1[n-1]=b1[n-1];
	sol2[n-1]=b2[n-1];
	for(l=n-2;l>=0;l--){
		sol1[l]=b1[l];
		sol1[l]-=A2[l]*sol1[l+1];

		sol2[l]=b2[l];
		sol2[l]-=A2[l]*sol2[l+1];
	}

	for(l=0;l<n;l++){
		b1[l]=sol1[l];
		b2[l]=sol2[l];
	}
	free(sol1); 
	free(sol2);
};
