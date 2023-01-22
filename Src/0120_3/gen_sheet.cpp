#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// ----- CONSTANTS ----
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)

#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// ----- Function Prototype  -----
float ran1(long *idum);
float gasdev(long *idum);
double dist(double *x, double *y);
int is_cross(double a1, double a2, double b1, double b2, double sig, double wd);
int is_cross2(double x1[2], double x2[2],	// sheet 1 (end points)
	double y1[2], double y2[2],	// sheet 2 (end points)
	double sigx, 			// characteristic length
	double sigy 			// characteristic length
);
double iprod(double x[2], double y[2]);
int lin_cross(
	double x1[2], double x2[2],
	double ya[2], double yb[2]
);
int is_cross_prd( 
	double x1[2],double x2[2], 
	double y1[2], double y2[2], 
	double sigx, double sigy,
	int iprd[2],
	double Xa[2], double Wd[2]
);


int main(){

	long idum=-2;	// seed for random number generator
	int i_max=100000;	// 
	int i,j;
	int iprd[2],irev[2];
	int ipt,n_sheet;
	int Nst,Npt,imb,ncrs;
	int iwall[4],nwall,imb_wall[4];
	int *nps,np_sum,nsig;
	double Xa[2],Xb[2],Wx[2],Xc[2],x0[2],n1[2];
	double xf1[2],xf2[2];
	double yf1[2],yf2[2];
	double *xx, *yy,*rr,*ths;
	double PI=4.0*atan(1.0);
	double th,wd,wd0,dr,sig,mass,wdh2,dw,std;
	double sigx,sigy,*sigs;
	double xwl[4],ywl[4],Xwl[4],Ywl[4];
	char cbff[128],fname[128],Dir0[128];
	double a1,a2,b1,b2,cost,sint;

	char fn1[128]="ptc.dat";
	char fn2[128]="sheet.dat";
	char fn3[128]="gen_sheet.out";
	char fn4[128]="bbox.dat";


//	-------- OPEN FIELS ----------
	FILE *fp0=fopen("gen_sheet.inp","r");
		fgets(Dir0,128,fp0);
		fgets(Dir0,128,fp0);
		int nchar=strlen(Dir0);
		Dir0[nchar-1]=int(NULL);


		strcpy(fname,Dir0); strcat(fname,fn1);
		puts(fname);
		FILE *fp1=fopen(fname,"w");

		strcpy(fname,Dir0); strcat(fname,fn2);
		puts(fname);
		FILE *fp2=fopen(fname,"w");

		strcpy(fname,Dir0); strcat(fname,fn3);
		puts(fname);
		FILE *fp3=fopen(fname,"w");

		strcpy(fname,Dir0); strcat(fname,fn4);
		puts(fname);
		FILE *fp4=fopen(fname,"w");
/*
	FILE *fp1=fopen("ptc.dat","w");
	FILE *fp2=fopen("sheet.dat","w");
	FILE *fp3=fopen("gen_sheet.out","w");
	FILE *fp4=fopen("bbox.dat","w");
	scanf("pause %s",cbff);
*/

	if(fp0 == NULL){
		printf("Can't find gen_sheet.inp \n");
		puts("  --> process terminated");
		exit(-1);
	}
	
//	-------	READ INPUT DATA  ----------
	fgets(cbff,128,fp0);
	fscanf(fp0,"%d %d\n",iprd,iprd+1);
//		iprd[0]=1;	// periodic B.C. (0:no, 1:yes) 
//		iprd[1]=1;	// periodic B.C. (0:no, 1:yes)
	printf("iprd=%d %d\n",iprd[0],iprd[1]);
	fgets(cbff,128,fp0);
	fscanf(fp0,"%lf %lf\n",Xa,Xa+1);	// bounding box
	fscanf(fp0,"%lf %lf\n",Xb,Xb+1);	// 
	
		Wx[0]=(Xb[0]-Xa[0]);	// box width 
		Wx[1]=(Xb[1]-Xa[1]);	// box height
		Xc[0]=0.5*(Xb[0]+Xa[0]); // center
		Xc[1]=0.5*(Xb[1]+Xa[1]); // center

		xwl[0]=Xc[0]; ywl[0]=Xa[1]; imb_wall[0]=-1;
		xwl[1]=Xb[0]; ywl[1]=Xc[1]; imb_wall[1]=-4;
		xwl[2]=Xc[0]; ywl[2]=Xb[1]; imb_wall[2]=-3;
		xwl[3]=Xa[0]; ywl[3]=Xc[1]; imb_wall[3]=-2;

		fprintf(fp4,"%lf %lf\n",Xa[0],Xa[1]);
		fprintf(fp4,"%lf %lf\n",Xb[0],Xa[1]);
		fprintf(fp4,"%lf %lf\n",Xb[0],Xb[1]);
		fprintf(fp4,"%lf %lf\n",Xa[0],Xb[1]);
		fprintf(fp4,"%lf %lf\n",Xa[0],Xa[1]);
		fclose(fp4);

	fgets(cbff,128,fp0);
	fscanf(fp0,"%d\n",&Nst);	// number of sheets
		printf("Number of Sheets=%d\n",Nst);
	fgets(cbff,128,fp0);
	fscanf(fp0,"%lf\n",&dw);	// segment size (particle spacing)
		printf("Segment size =%lf\n",dw);
	fgets(cbff,128,fp0);
	fscanf(fp0,"%lf \n",&wd0); // width (mean)
	fgets(cbff,128,fp0);
	fscanf(fp0,"%lf\n",&std); // standard dev for wd0
		std/=100.0;
		std/=3.0;
	fgets(cbff,128,fp0);
	fscanf(fp0,"%lf %d\n",&mass,&imb); // mass and mibile/immobile
	fgets(cbff,128,fp0);
	nwall=0;
	for(i=0;i<4;i++){
		fscanf(fp0,"%d",iwall+i);
		printf("iwall[i]=%d\n",iwall[i]);
		nwall+=iwall[i];
	}
	printf("Number of Walls=%d\n",nwall);
	fgets(cbff,3,fp0);
	fgets(cbff,128,fp0);
	fscanf(fp0,"%d\n",&nsig);
	fgets(cbff,128,fp0);
	sigs=(double *)malloc(sizeof(double)*Nst);
	int *null_wd=(int *)malloc(sizeof(int)*Nst);

	if(nsig==0){
		int itmp;
		fscanf(fp0,"%lf %d\n",&sig,&itmp);
		for(i=0;i<Nst;i++){
		       sigs[i]=sig; 
		       null_wd[i]=itmp;
		}
	}else{
		for(i=0;i<Nst;i++) fscanf(fp0,"%lf %d\n",sigs+i,null_wd+i);
	};
	fclose(fp0);
	for(i=0;i<Nst;i++){
		printf("sig[%d]=%lf mass point? %d\n",i,sigs[i],null_wd[i]);
		sigs[i]*=pow(2.0,0.1666666667);	// minimum potential distance 
	};

//		wall particle coordinates
	for(i=0;i<4;i++){
		Xwl[i]=xwl[i]; Ywl[i]=ywl[i]; 
	}
	//sig*=pow(2.0,0.1666666667);	// minimum potential distance 
	sig=sigs[0];
	Ywl[0]-=0.5*sig;
	Ywl[2]+=0.5*sig;
	Xwl[1]+=0.5*sig;
	Xwl[3]-=0.5*sig;
//	-----------------------------------

	xx=(double *)malloc(sizeof(double)*Nst);
	yy=(double *)malloc(sizeof(double)*Nst);
	rr=(double *)malloc(sizeof(double)*Nst);
	ths=(double *)malloc(sizeof(double)*Nst);
	nps=(int *)malloc(sizeof(int)*Nst);

	FILE *tmp=fopen("temp.out","w");

	fprintf(tmp,"#REV\n");
	fprintf(tmp,"%lf %lf\n",Xa[0],Xa[1]);
	fprintf(tmp,"%lf %lf\n",Xa[0]+Wx[0],Xa[1]);
	fprintf(tmp,"%lf %lf\n",Xa[0]+Wx[0],Xa[1]+Wx[1]);
	fprintf(tmp,"%lf %lf\n",Xa[0],Xa[1]+Wx[1]);
	fprintf(tmp,"%lf %lf\n",Xa[0],Xa[1]);
	fprintf(tmp,"\n");

	fprintf(fp3,"# Position (x,y), Width wd, Angle th(deg), Particles: Number Npt, Spacing dw\n");

	n_sheet=0;
	np_sum=0;
	for(i=0;i<i_max;i++){	// START_Trial 

		x0[0]=Xa[0]+Wx[0]*ran1(&idum);
		x0[1]=Xa[1]+Wx[1]*ran1(&idum);
		th=ran1(&idum)*PI;
		wd=wd0*(1.+std*gasdev(&idum));
		if(wd < 0.0) wd=2.*dw;

		if(null_wd[i]==1) wd=0.0;

		Npt=ceil(wd/dw)+1;
		if(Npt <= 0) Npt=1;
		wd=(Npt-1)*dw;

		if(fabs(wd) < 1.e-06) wd=dw*1.e-06;

		yf1[0]=x0[0]-wd*0.5*cos(th);
		yf1[1]=x0[1]-wd*0.5*sin(th);
		yf2[0]=x0[0]+wd*0.5*cos(th);
		yf2[1]=x0[1]+wd*0.5*sin(th);

		ncrs=0;
		sigy=sigs[n_sheet];
		for(j=0;j<n_sheet;j++){	// Start_Sheet
			
			sigx=sigs[j];
			xf1[0]=xx[j]-rr[j]*cos(ths[j]);
			xf1[1]=yy[j]-rr[j]*sin(ths[j]);
			xf2[0]=xx[j]+rr[j]*cos(ths[j]);
			xf2[1]=yy[j]+rr[j]*sin(ths[j]);

			//ncrs=is_cross_prd(xf1,xf2,yf1,yf2,sig,iprd,Xa,Wx);
			ncrs=is_cross_prd(xf1,xf2,yf1,yf2,sigx,sigy,iprd,Xa,Wx);
			if(ncrs==1) break;

			//ncrs=is_cross_prd(xf1,xf2,yf1,yf2,sig,iprd,Xa,Wx);
			ncrs=is_cross_prd(yf1,yf2,xf1,xf2,sigy,sigx,iprd,Xa,Wx);
			if(ncrs==1) break;

			ncrs=lin_cross(xf1,xf2,yf1,yf2);
			if(ncrs==1) break;

		}	// End_Sheet

		if(ncrs ==0){
			fprintf(tmp,"#n_sheet=%d\n",n_sheet);
			fprintf(tmp,"%lf %lf\n",yf1[0],yf1[1]);
			fprintf(tmp,"%lf %lf\n",yf2[0],yf2[1]);
			fprintf(tmp,"\n");
		}
		
		for(j=0;j<4;j++){
			if(iwall[j]==0) continue;

			n1[0]=1.0;
			n1[1]=1.0;

			if(j%2 == 0) n1[0]=0.0;
			if(j%2 == 1) n1[1]=0.0;

			xf1[0]=x0[0]-0.5*wd*cos(th) - xwl[j];
			xf1[1]=x0[1]-0.5*wd*sin(th) - ywl[j];
			b1=xf1[0]*n1[0]+xf1[1]*n1[1];

			xf2[0]=x0[0]+0.5*wd*cos(th) - xwl[j];
			xf2[1]=x0[1]+0.5*wd*sin(th) - ywl[j];
			b2=xf2[0]*n1[0]+xf2[1]*n1[1];

			if(fabs(b1) <= sig){
				ncrs=1;
				break;
			}else if(fabs(b2) <= sig){
				ncrs=1;	
				break;
			}else if(b1*b2 > 0.0){
				continue;
			}
			ncrs=1;
		}

		if(ncrs == 1 ) continue;

		xx[n_sheet]=x0[0];
		yy[n_sheet]=x0[1];
		rr[n_sheet]=wd*0.5;
		ths[n_sheet]=th;
		nps[n_sheet]=Npt;
		np_sum+=Npt;
		n_sheet++;

		th=th/PI*180.0;
		//fprintf(fp3,"%d: (x,y)=(%5.2lf, %5.2lf),wd=%5.2lf, th=%6.2lf(deg), Npt=%d, dw=%5.2lf\n",
		fprintf(fp3,"%5.2lf, %5.2lf, %5.2lf, %6.2lf, %d, %5.2lf\n", x0[0],x0[1],wd,th,Npt,dw);
		
		printf("#i=%d, sheet %d\n",i,n_sheet);
		if(n_sheet ==Nst) break;
	} // END_Trial
	fclose(tmp);
	if(i==i_max){	
		puts("Can't pack the specified number of sheets!!");
		puts(" Increse the number of trials i_max or use larger bounding box");
		puts(" process terminated");
		exit(-1);
	}

	printf("Number of PARTICLES : %d\n",np_sum+nwall);
	printf("Number of SHEETS    : %d\n",Nst);
	fprintf(fp1,"# np (number of particles )\n");
	fprintf(fp1,"%d\n",np_sum+nwall);
	fprintf(fp1,"# (x,y), mass, type (1:mobile, -1,-2: rev(horiz, vert), sig\n");
	fprintf(fp2,"# Number of sheets\n");
	fprintf(fp2,"%d\n",Nst+nwall);

	for(i=0;i<nwall;i++){
		fprintf(fp2,"# Sheet: %d\n",i+1);
		fprintf(fp2," %d\n",1);
		fprintf(fp2," %d\n",i+1);
	}
	for(i=0;i<4;i++){
		if(iwall[i]==1){
			 fprintf(fp1,"%24.16le %24.16le %24.16le %d %d %d\n",Xwl[i],Ywl[i],1.0,imb_wall[i],0,0);
			 fprintf(fp1,"\n");
		}
	}
	
	tmp=fopen("temp2.out","w");
	ipt=0;
	for(i=0;i<Nst;i++){
		fprintf(fp2,"# Sheet: %d\n",i+nwall+1);
		fprintf(fp2," %d\n",nps[i]);
		for(j=0;j<nps[i];j++){
			xf1[0]=xx[i]+(dw*j-rr[i])*cos(ths[i]);
			xf1[1]=yy[i]+(dw*j-rr[i])*sin(ths[i]);


			irev[0]=0;
			irev[1]=0;
			if(iprd[0]==1){
				if(xf1[0] < Xa[0]) irev[0]--;
				if(xf1[0] > Xa[0]+Wx[0]) irev[0]++;
			}
			if(iprd[0]==1){
				if(xf1[1] < Xa[1]) irev[1]--;
				if(xf1[1] > Xa[1]+Wx[1]) irev[1]++;
			}
			xf1[0]-=Wx[0]*irev[0]*iprd[0];
			xf1[1]-=Wx[1]*irev[1]*iprd[1];

			fprintf(tmp,"%lf %lf %d %d\n",xf1[0],xf1[1],irev[0],irev[1]);
			ipt++;
			sig=sigs[i]/pow(2.0,0.166666666667);
			fprintf(fp1,"%24.16le %24.16le %24.16le %d %d %d %le\n",xf1[0],xf1[1],mass,imb,irev[0],irev[1],sig);
			//fprintf(fp1,"%24.16le %24.16le %24.16le %d %d %d\n",xf1[0],xf1[1],mass,imb,irev[0],irev[1]);
			fprintf(fp2," %d",ipt+nwall);

		}
		fprintf(fp1,"\n");
		fprintf(fp2,"\n");
	}
	
	fclose(tmp);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);

	return(0);
}

// ----------------RANDOM NUMBER GENERATOR---------------
float gasdev(long *idum){
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if(iset ==0){
		do{
		v1=2.0*ran1(idum)-1.0;
		v2=2.0*ran1(idum)-1.0;
		rsq=v1*v1+v2*v2;
		}while(rsq>=1.0 || rsq == 0.0);
	fac=sqrt(-2.0*log(rsq)/rsq);
	gset=v1*fac;
	iset=1;
	return v2*fac;
	}else{
		iset=0;
		return gset;
	}
}

float ran1(long *idum){
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if(*idum <= 0 || !iy){
		if(-(*idum) < 1) *idum=1;
		else *idum =-(*idum);
		for(j=NTAB+7;j>=0;j--){
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if(*idum < 0) *idum += IM;
			if(j < NTAB) iv[j]= *idum;
		}
		iy=iv[0];
 	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if(*idum <0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]= *idum;
	if((temp=AM*iy) >RNMX) return RNMX;
	else return temp;
}
// ------------------------------------------------------------
double dist(double *x, double *y){
	
	double dx,dy;
	dx=x[0]-y[0];
	dy=x[1]-y[1];
	return sqrt(dx*dx+dy*dy);
}
int is_cross2(
	double x1[2], double x2[2],	// sheet 1 (end points)
	double y1[2], double y2[2],	// sheet 2 (end points)
	double sigx, double sigy	// characteristic length
){

	double sig;
	int i,j,ncrs=0,isum;
	double ys[5][2],ya[2],yb[2];
	double lenx,leny,len,txny;
	double tx[2],ty[2],ny[2],u1[2],u;
	double xi1[2],xi2[2];
	double fac=2.0,L1,L2;
	double p1[2],p2[2]; 

	double xc[2],yc[2],wdx,wdy;

	xc[0]=0.5*(x1[0]+x2[0]);
	xc[1]=0.5*(x1[1]+x2[1]);
	yc[0]=0.5*(y1[0]+y2[0]);
	yc[1]=0.5*(y1[1]+y2[1]);
	wdx=dist(x1,x2)*0.5;
	wdy=dist(y1,y2)*0.5;

	sig=sigx;
	//if(dist(xc,yc) > wdx+wdy+sig*2*fac) return(0);
	if(dist(xc,yc) > wdx+wdy+(sigx+sigy)*fac) return(0);

	

	for(i=0;i<2;i++){
		 ty[i]=y2[i]-y1[i];
		 tx[i]=x2[i]-x1[i];
	}
	lenx=sqrt(iprod(tx,tx));
	leny=sqrt(iprod(ty,ty));
	for(i=0;i<2;i++){
		 tx[i]/=lenx;
		 ty[i]/=leny;
	}
	ny[0]=-ty[1];
	ny[1]= ty[0];

	for(i=0;i<2;i++){
		ys[0][i]=y1[i]+fac*sigy*(-ty[i]-ny[i]);
		ys[1][i]=y2[i]+fac*sigy*( ty[i]-ny[i]);
		ys[2][i]=y2[i]+fac*sigy*( ty[i]+ny[i]);
		ys[3][i]=y1[i]+fac*sigy*(-ty[i]+ny[i]);
		ys[4][i]=ys[0][i];

		p1[i]=x1[i]-ys[0][i];
		p2[i]=x2[i]-ys[0][i];
	}
	L1=dist(ys[1],ys[0]);
	L2=dist(ys[3],ys[0]);

	isum=0;
	xi1[0]=iprod(p1,ty)/L1;
	xi1[1]=iprod(p1,ny)/L2;
	if(xi1[0] >= 0.0) isum++;
	if(xi1[0] <= 1.0) isum++;
	if(xi1[1] >= 0.0) isum++;
	if(xi1[1] <= 1.0) isum++;
	if(isum==4) ncrs=1;

	isum=0;
	xi2[0]=iprod(p2,ty)/L1;
	xi2[1]=iprod(p2,ny)/L2;
	if(xi2[0] >= 0.0) isum++;
	if(xi2[0] <= 1.0) isum++;
	if(xi2[1] >= 0.0) isum++;
	if(xi2[1] <= 1.0) isum++;
	if(isum==4) ncrs=1;

	for(i=0;i<4;i++){
//		for(j=0;j<2;j++){ ya[j]=ys[i][j]; yb[j]=ys[i+1][j]; }
//		if(lin_cross(x1,x2,ys[i],ys[i+1])==1){
//			ncrs=1;
//			break;
//		}
//		ncrs+=lin_cross(x1,x2,ya,yb);
		ncrs+=lin_cross(x1,x2,ys[i],ys[i+1]);
	}
	if(ncrs>0) ncrs=1;

	return(ncrs);

};
int lin_cross(
	double x1[2], double x2[2],
	double ya[2], double yb[2]
){
	int i,ncrs=0;
	int isum;
	double u1[2],tx[2],ty[2],nx[2],ny[2];
	double lenx,leny,txny,tynx,u,v;
	double eps=1.e-08;

	for(i=0;i<2;i++){
		tx[i]=x2[i]-x1[i];
		ty[i]=yb[i]-ya[i];
		u1[i]=x1[i]-ya[i];
	}
	lenx=sqrt(iprod(tx,tx));
	leny=sqrt(iprod(ty,ty));
	for(i=0;i<2;i++){
		tx[i]/=lenx;
		ty[i]/=leny;
	}
	nx[0]=-tx[1];
	nx[1]= tx[0];
	ny[0]=-ty[1];
	ny[1]= ty[0];

	txny=iprod(tx,ny);
	tynx=iprod(ty,nx);

	if(fabs(txny) <eps){
		if(txny <0.0){
 			txny=-eps;
		}else{
 			txny=eps;
		}
	}

	if(fabs(tynx) <eps){
		if(tynx <0.0){
 			tynx=-eps;
		}else{
 			tynx=eps;
		}
	}
	u=-iprod(u1,ny)/txny;
	v= iprod(u1,nx)/tynx;
//	printf("#u,v=%lf %lf %lf %lf\n",u,v,u/lenx,v/leny);

	isum=0;
	if( (u >= 0.0) && (u<=lenx)) isum++;
	if( (v >= 0.0) && (v<=leny)) isum++;
	if(isum==2) ncrs=1;

	return(ncrs);
	
};

double iprod(
	double x[2], double y[2]
){
	return(x[0]*y[0]+x[1]*y[1]);
};

int is_cross(double a1, double a2, double b1, double b2, double sig, double wd){

	int ncrs=0;
	double wdh2=wd*0.5+sig;
/*
	if(fabs(b1) <= sig){
		if( fabs(a1) <= wdh2) ncrs=1; 
	}else if(fabs(b2) <= sig){
		if( fabs(a2) <= wdh2) ncrs=1; 
	}else if(b1*b2 > 0.0){
			ncrs=-1;
	}

	if( fabs(a1) <= wdh2) ncrs=1; 
	if( fabs(a2) <= wdh2) ncrs=1; 
*/

	if(fabs(b1) <= sig){
		if( fabs(a1) <= wdh2) ncrs+=1; 
	}
	if(fabs(b2) <= sig){
		if( fabs(a2) <= wdh2) ncrs+=1; 
	}

	if(b1*b2 <= 0.0){
		if( fabs(a1) <= wdh2) ncrs+=1; 
		if( fabs(a2) <= wdh2) ncrs+=1; 
	}


	if(ncrs >0) ncrs=1;

	return ncrs;
}

int is_cross_prd( 
	double x1[2],double x2[2], 
	double y1[2], double y2[2], 
	double sigx, double sigy,
	int iprd[2],
	double Xa[2], double Wd[2]
){

	int ncrs=0;
	int i,j;
	int i1=0,i2=1;
	int j1=0,j2=1;
	double y1d[2],y2d[2];


	if(iprd[0]==1){
		i1=-1;
		i2=2;
	} 
	if(iprd[1]==1){
		j1=-1;
		j2=2;
	} 

	for(i=i1;i<i2;i++){
		y1d[0]=y1[0]+i*Wd[0];
		y2d[0]=y2[0]+i*Wd[0];
	for(j=j1;j<j2;j++){
		y1d[1]=y1[1]+j*Wd[1];
		y2d[1]=y2[1]+j*Wd[1];
		ncrs=is_cross2(x1,x2,y1d,y2d,sigx,sigy);
		if(ncrs==1) return(ncrs);
	}
	}

	return(ncrs);
	
};


