#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dem.h"

using namespace std;

double STF1(PRTCL p1, PRTCL p2, double *Fn, double r0, double K1, REV rev, double Sab[2][2]){
	int i,j;
	double n12[2],r2,r;
	double UE=0.0,Fr=0.0;
	Vec2 x1,x2,xd;

	Fn[0]=0.0;
	Fn[1]=0.0;

	if(p1.mobile != 1) return(0.0);
	if(p2.mobile != 1) return(0.0);

	x1=rev.unfold(p1.x, p1.irev);
	x2=rev.unfold(p2.x, p2.irev);
	xd=vdiff(x1,x2);

	n12[0]=xd.x[0];
	n12[1]=xd.x[1];
	r2=n12[0]*n12[0]+n12[1]*n12[1];	

	r=sqrt(r2);
	n12[0]/=r;	
	n12[1]/=r;	

	UE=K1*(r-r0);
	Fn[0]=-n12[0]*UE;
	Fn[1]=-n12[1]*UE;

	Fr=UE*r/3.0;
	Sab[0][0]=n12[0]*n12[0]*(-Fr);
	Sab[1][0]=n12[1]*n12[0]*(-Fr);
	Sab[1][1]=n12[1]*n12[1]*(-Fr);
	//printf("r0=%lf F=%lf %lf\n",r0,Fn[0],Fn[1]);


	UE*=(r-r0);
	UE*=0.5;


	return(UE);
}

double LJ(
	PRTCL p1, PRTCL p2,// particles
	double *Fn, // VanDerWaals Force Vector
	double sig, double Eps,  // LJ parameters
	REV rev,	// Unit cell
	int *ofst1, 	// Unit cell ID (offset)
	int *ofst2,	// Unit cell ID (offset)
	double Sab[2][2], 
	int ncut	// cutoff distance (rcut=sig*ncut)
){
	double r,r2,sr,sr2,sr6,sr12,F,Fr;
	double n12[2];
	double s2=sig*sig;
	double rcut=pow(2.0,0.1666666667)*sig;
	//double rcut3=rcut*3.0; 
	double rcut3=rcut*ncut; 
	double UE=0.0;
	int mb1=p1.mobile;
	int mb2=p2.mobile;

	double Wd[2];
	Vec2 x1,x2,xd;
	
	Fn[0]=0.0; Fn[1]=0.0;

	if( (mb1!=1) && (mb2!=1)) return(UE);

	x1=rev.unfold(p1.x,ofst1);
	x2=rev.unfold(p2.x,ofst2);
	xd=vdiff(x1,x2);

	n12[0]=xd.x[0];
		if(p1.mobile==-1 || p1.mobile==-3) n12[0]=0.0;
		if(p2.mobile==-1 || p2.mobile==-3) n12[0]=0.0;
	n12[1]=xd.x[1];
		if(p1.mobile==-2 || p1.mobile==-4) n12[1]=0.0;
		if(p2.mobile==-2 || p2.mobile==-4) n12[1]=0.0;

	r2=n12[0]*n12[0]+n12[1]*n12[1];	
	r=sqrt(r2);
	n12[0]/=r;	
	n12[1]/=r;	

	sr2=s2/r2;

	sr=sqrt(sr2);
	sr6=sr2*sr2*sr2;
	sr12=sr6*sr6;

	F=Eps*(2.*sr12-sr6)*sr;
	UE=sr12-sr6;
	if( (mb1!=1)||(mb2 !=1) ){
		if(r >rcut){
			F=0.0;
			UE=0.0;
		}
	} 

	if(r > rcut3){
		F=0.0;
		UE=0.0;
	}

	Fn[0]=F*n12[0];
	Fn[1]=F*n12[1];

	Fr=F*r/3.0;

	Sab[0][0]=n12[0]*n12[0]*Fr;
	Sab[1][0]=n12[1]*n12[0]*Fr;
	Sab[1][1]=n12[1]*n12[1]*Fr;

	return(UE);
}



double VDF(
	REV rev, 
	SUBCELL *sbcll, 
	PRTCL *PTC, 
	int iprd[2],
	double sig, double Eps, 
	double Sab[2][2]
){
	int i,ix,iy,ip,ipt;
	int j,jx,jy,jp,jpt;
	int Jx,Jy;
	int npi,npj;
	int iofst[2],jofst[2];
	int nsub=rev.nsub;
	SUBCELL sbci,sbcj;
	int Nx[2];
	int idiff,ndiff,isgn;
	double dFn[2],UE=0.0,dUE;
	double dSab[2][2];
	double sigi,sigj;
	double dsig[2];

	double r0_tmp=1.0;
	dsig[0]=0.0; dsig[1]=0.0;

	Nx[0]=rev.Nh[0];
	Nx[1]=rev.Nh[1];

	iofst[0]=0; iofst[1]=0;

	Sab[0][0]=0.0; 
	Sab[1][0]=0.0; 
	Sab[1][1]=0.0; 
	for(i=0;i<nsub;i++){	// SUBCELL_i
		sbci=sbcll[i];
		npi=sbci.np;
		ix=i/Nx[1];	// 2D cell index 1
		iy=i%Nx[1];	// 2D cell index 2
		for(ip=0; ip<npi; ip++){ // PARTICLE_i 
			ipt=sbci.list[ip];
			PTC[ipt].UE[0]=0.0;
			PTC[ipt].UE[1]=0.0;
		for(jx=-1; jx<2; jx++){
			jofst[0]=0; 
			Jx=ix+jx;
			if(Jx < 0){
				Jx+=Nx[0];
				jofst[0]=-1;
			}
			if(Jx >= Nx[0]){
				Jx-=Nx[0];
				jofst[0]=1;
			}
			if(jofst[0] !=0 && iprd[0]==0) continue;
		for(jy=-1; jy<2; jy++){
			Jy=iy+jy;
			jofst[1]=0;
			if(Jy < 0){
				Jy+=Nx[1];
				jofst[1]=-1;
			}
			if(Jy >= Nx[1]){
				 Jy-=Nx[1];
				jofst[1]=1;
			}
			if(jofst[1] !=0 && iprd[1]==0) continue;

			j=Jx*Nx[1]+Jy; // get 1D subcell index
			sbcj=sbcll[j];	
			npj=sbcj.np;
		for(jp=0;jp<npj; jp++){ // PARTICLE_j
			jpt=sbcj.list[jp];
			if(ipt==jpt) continue;
			idiff=abs(PTC[ipt].ipt-PTC[jpt].ipt);

			//sigj=PTC[jpt].get_sig(PTC[ipt].x,&isgn,dsig);
			//sigi=PTC[ipt].get_sig(PTC[jpt].x,&isgn,dsig);
			sigj=PTC[jpt].get_sig_iso(PTC[ipt].x,&isgn,dsig);
			sigi=PTC[ipt].get_sig_iso(PTC[jpt].x,&isgn,dsig);
			sig=0.5*(sigi+sigj);
			ndiff=int(3.*sig/r0_tmp);
			if(PTC[ipt].ist==PTC[jpt].ist && idiff < ndiff) continue; 

			dUE=LJ(PTC[ipt],PTC[jpt],dFn,sig,Eps,rev,iofst,jofst,dSab,3);
//			printf("Fn=%le %le\n",dFn[0],dFn[1]);

			UE+=dUE;
			PTC[ipt].UE[isgn]+=dUE;
			PTC[ipt].F[0]+=dFn[0];
			PTC[ipt].F[1]+=dFn[1];
			Sab[0][0]+=dSab[0][0];
			Sab[1][0]+=dSab[1][0];
			Sab[1][1]+=dSab[1][1];
		}	// PARTILCE_j
		}
		}
		} // PARTILCE_i
	} // SUBCELL_i
	return(UE);
};
//-----------------------------------------------------------------------------
double VarUE_mu(
	REV rev, 
	SUBCELL *sbcll, 
	PRTCL *PTC, 
	int iprd[2],
	double sig, double Eps,
	int ipt,	// perturbed particles 
	int iside, 	// side (0:top, 1: bottom)
	double var_sig	// amount of water taken/supplied
){

	int i,ix,iy;
	int j,jx,jy,jp,jpt;
	int Jx,Jy;
	int npi,npj;
	int iofst[2],jofst[2];
	int nsub=rev.nsub;
	SUBCELL sbci,sbcj;
	int Nx[2];
	Nx[0]=rev.Nh[0];
	Nx[1]=rev.Nh[1];
	int idiff,ndiff,isgn;
	double UE,dUE,uvar;
	double dFn[2],dSab[2][2];
	double sigi,sigj;
	double dsig0[2],dsig[2];

	double r0_tmp=1.0;


	i=PTC[ipt].sbcll; // local particle number in
	sbci=sbcll[i];	// subcell particle belongs 
	npi=sbci.np;	// number of particles in the subcell
	ix=i/Nx[1];	// 2D cell index 1
	iy=i%Nx[1];	// 2D cell index 2
	UE=0.0;	
	uvar=0.0;

	dsig0[0]=0.0; dsig0[1]=0.0;
	dsig[0]=0.0; dsig[1]=0.0;
	dsig[iside]=var_sig;
	for(jx=-1; jx<2; jx++){	// Subcell_jx
		jofst[0]=0; 
		Jx=ix+jx;
		if(Jx < 0){
			Jx+=Nx[0];
			jofst[0]=-1;
		}
		if(Jx >= Nx[0]){
			Jx-=Nx[0];
			jofst[0]=1;
		}
		if(jofst[0] !=0 && iprd[0]==0) continue;
	for(jy=-1; jy<2; jy++){	// Subcell_jy
		Jy=iy+jy;
		jofst[1]=0;
		if(Jy < 0){
			Jy+=Nx[1];
			jofst[1]=-1;
		}
		if(Jy >= Nx[1]){
			 Jy-=Nx[1];
			jofst[1]=1;
		}
		if(jofst[1] !=0 && iprd[1]==0) continue;

		j=Jx*Nx[1]+Jy; // get 1D subcell index
		sbcj=sbcll[j];	
		npj=sbcj.np;
		for(jp=0;jp<npj; jp++){ // PARTICLE_jp 
			jpt=sbcj.list[jp];
			if(ipt==jpt) continue;
			idiff=abs(PTC[ipt].ipt-PTC[jpt].ipt);
			sigj=PTC[jpt].get_sig_iso(PTC[ipt].x,&isgn,dsig0);
			sigi=PTC[ipt].get_sig_iso(PTC[jpt].x,&isgn,dsig);
			if(isgn!=iside) continue;
			sig=0.5*(sigi+sigj);
			ndiff=int(3.*sig/r0_tmp);
			if(PTC[ipt].ist==PTC[jpt].ist && idiff < ndiff) continue; 
			dUE=LJ(PTC[ipt],PTC[jpt],dFn,sig,Eps,rev,iofst,jofst,dSab,3);
			UE+=dUE;
		}	// End_PARTICLE_jp
	}	// End_Subcell_jy
	}	// End_Subcell_jx
	uvar+=(UE-PTC[ipt].UE[iside]);
	return(uvar);

};
double VarUE(
	REV rev, 
	SUBCELL *sbcll, 
	PRTCL *PTC, 
	int iprd[2],
	double sig, double Eps,
	int ipts[2],	// perturbed particles 
	int isds[2], 	// side (0:top, 1: bottom)
	double var_sig,	// amount of water taken/supplied
	double dUE_try[2]
){
	int i,ix,iy,ip,ipt;
	int j,jx,jy,jp,jpt,jpt0;
	int Jx,Jy;
	int npi,npj;
	int iofst[2],jofst[2];
	int nsub=rev.nsub;
	SUBCELL sbci,sbcj;
	int Nx[2];
	int idiff,ndiff,isgn,iside;
	double UE[2],dUE;
	double dFn[2],dSab[2][2];
	double sigi,sigj;
	double dsig0[2],dsig1[2],dsig2[2];
	double *dsigi,*dsigj;
	double *sigs[2];

	double r0_tmp=1.0;

	for(j=0;j<2;j++){
		dsig0[j]=0.0;
		dsig1[j]=0.0;
		dsig2[j]=0.0;
	};

	dsig1[isds[0]]= var_sig;
	dsig2[isds[1]]=-var_sig;

	sigs[0]=dsig1;
	sigs[1]=dsig2;

	Nx[0]=rev.Nh[0];
	Nx[1]=rev.Nh[1];

	iofst[0]=0; iofst[1]=0;

	double uvar=0.0;
	for(int k=0;k<2;k++){	// perturbed pair of particles
		dUE_try[k]=0.0;

		ipt=ipts[k];
		jpt0=ipts[(k+1)%2];
		iside=isds[k];
		dsigi=sigs[k];

		i=PTC[ipt].sbcll;
		sbci=sbcll[i];
		npi=sbci.np;
		ix=i/Nx[1];	// 2D cell index 1
		iy=i%Nx[1];	// 2D cell index 2
		UE[0]=0.0;
		UE[1]=0.0;
		for(jx=-1; jx<2; jx++){
			jofst[0]=0; 
			Jx=ix+jx;
			if(Jx < 0){
				Jx+=Nx[0];
				jofst[0]=-1;
			}
			if(Jx >= Nx[0]){
				Jx-=Nx[0];
				jofst[0]=1;
			}
			if(jofst[0] !=0 && iprd[0]==0) continue;
		for(jy=-1; jy<2; jy++){
			Jy=iy+jy;
			jofst[1]=0;
			if(Jy < 0){
				Jy+=Nx[1];
				jofst[1]=-1;
			}
			if(Jy >= Nx[1]){
				 Jy-=Nx[1];
				jofst[1]=1;
			}
			if(jofst[1] !=0 && iprd[1]==0) continue;

			j=Jx*Nx[1]+Jy; // get 1D subcell index
			sbcj=sbcll[j];	
			npj=sbcj.np;
		for(jp=0;jp<npj; jp++){ // PARTICLE_j
			jpt=sbcj.list[jp];
			if(ipt==jpt) continue;
			idiff=abs(PTC[ipt].ipt-PTC[jpt].ipt);

			//if(PTC[ipt].ist==PTC[jpt].ist && idiff < 2) continue; 
			dsigj=dsig0;
			if(jpt==jpt0) dsigj=sigs[(k+1)%2];
			//sigj=PTC[jpt].get_sig(PTC[ipt].x,&isgn,dsigj);
			//sigi=PTC[ipt].get_sig(PTC[jpt].x,&isgn,dsigi);
			sigj=PTC[jpt].get_sig_iso(PTC[ipt].x,&isgn,dsigj);
			sigi=PTC[ipt].get_sig_iso(PTC[jpt].x,&isgn,dsigi);
			if(isgn!=iside)	continue;
			sig=0.5*(sigi+sigj);
			ndiff=int(3.*sig/r0_tmp);
			
			if(PTC[ipt].ist==PTC[jpt].ist && idiff < ndiff) continue; 
			dUE=LJ(PTC[ipt],PTC[jpt],dFn,sig,Eps,rev,iofst,jofst,dSab,3);
			UE[isgn]+=dUE;
		}	// PARTILCE_j
		}
		}
		uvar+=(UE[iside]-PTC[ipt].UE[iside]);
		//dUE_try[iside]+=(UE[iside]-PTC[ipt].UE[iside]);
		dUE_try[k]+=(UE[iside]-PTC[ipt].UE[iside]);
	}
	return(uvar);

};
//-----------------------------------------------------------------------------
double VDF_L(
	REV rev, 
	PRTCL *PTC, 
	int iprd[2],
	double Eps, 
	double Sab[2][2],
	int np, int npl, int *indx
){
	int ipt,jpt,jx,jy;
	double UE=0.0,Sig;
	double dFn[2],dSab[2][2];
	int iofst[2],jofst[2];
	int ipl,jpl;

	//for(ipt=0;ipt<npl;ipt++){
	for(ipl=0;ipl<npl;ipl++){
		ipt=indx[ipl];
		Sig=PTC[ipt].sig;
		iofst[0]=0; 
		iofst[1]=0;
		for(jpt=0;jpt<np;jpt++){
			if(PTC[jpt].sig > 10.0) continue; 
		for(jx=-1;jx<2;jx++){
			jofst[0]=jx;
		for(jy=-1;jy<2;jy++){
			jofst[1]=jy;

			UE+=LJ(PTC[ipt],PTC[jpt],dFn,Sig,Eps,rev,iofst,jofst,dSab,1);
			PTC[ipt].F[0]+=dFn[0];
			PTC[ipt].F[1]+=dFn[1];

			PTC[jpt].F[0]-=dFn[0];
			PTC[jpt].F[1]-=dFn[1];

			Sab[0][0]+=dSab[0][0];
			Sab[1][0]+=dSab[1][0];
			Sab[1][1]+=dSab[1][1];
		}
		}
		}

		for(jx=-1;jx<2;jx++){
			jofst[0]=jx;
		for(jy=-1;jy<2;jy++){
			jofst[1]=jy;
		for(jpl=0;jpl<npl;jpl++){
			jpt=indx[jpl];
			if(ipt==jpt && jx*jx+jy*jy ==0) continue;
			UE+=LJ(PTC[ipt],PTC[jpt],dFn,Sig,Eps,rev,iofst,jofst,dSab,1);
			PTC[ipt].F[0]+=dFn[0];
			PTC[ipt].F[1]+=dFn[1];

			Sab[0][0]+=dSab[0][0];
			Sab[1][0]+=dSab[1][0];
		}
		}
		}
	}

	return(UE);
};
