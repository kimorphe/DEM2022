double dist(double *x, double *y);
void show_msg(char *fname);


class Vec2{
	public:
	double x[2];
	void set(double,double);	
	void set(double*);	
	Vec2 times(double);	
	Vec2 div(double);	
	void times_me(double);	
	void div_me(double);	
	double len();
	void print();
	private:
};
double iprod(Vec2 a,Vec2 b);
double iprod(double a[2], double b[2]);
Vec2 vsum(Vec2 a, Vec2 b);
Vec2 vdiff(Vec2 a, Vec2 b);
void vsum(Vec2 a, Vec2 b, Vec2 c);
void vdiff(Vec2 a, Vec2 b, Vec2 c);

class MatTriDiag{	// triple diagonal matrix
	public:
		double *A0;	// A[i-1,i]
		double *A1;	// A[i,i]
		double *A2;	// A[i,i+1]
		double *y;	// solution vector (used when LinSolve(b) is called) 
		int n;		// matrix size (n x n)
		void init(int N); // initialize class instance (must be used firstly)
		void print(); // print matrix components
		// Linear Equation Solvers (Gaussian Elimination)
		// Note that the coeeficnet matrices are destroyed during the solution process
		void LinSolve(double *b); // take one r.h.s. vector --> solution =y
		void LinSolve(double *b1, double *b2); // take two r.h.s. vectors --> solutions b1,b2
	private:
	protected:
};

class Curve2D{	// plane (2D) curve 
	public:
		int np; // number of points
		double *x,*y; // coordinate vectors x,y
		double *bx,*cx,*dx; // cubic spline coefficients for x
		double *by,*cy,*dy; // cubic spline coefficients for y
		void init(int N); // initialize instance (must be called firstly)
		void print(); // print data 
		void print(char *fname); // print data 
		void spline(); // obtain spline coefficients
		double intplx(double s); // interpolate x(s) for s in (0, np-1)
		double intply(double s); // interpolate y(s) for s in (0, np-1)
		double dxds(double s); // interpolate dx/ds 
		double dyds(double s); // interpolate dy/ds
		void smooth(int nsmp2);
		~Curve2D();
	private:
	protected:
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
		//void trend(int type,double amp_wv);
		double eval(double xval);
		void load(char fn[128]);
		double y2x(double yval);
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
		void change_unit(double unit);
		void print();
		double nw_min,nw_max;
	private:
	protected:
};

class CNTRL{	// DEM implementation controle parameters
	public:
		char Dir[128];	// Output Data Directory
		int Nout;	// number of output steps 
		int Nt;		// total time steps
		int np;		// number of particles
		int nst;	// number of clay molecules
		int rstat;	// restart (0:No, 1:yes)
		char fnrst[128];// I.C data file for restart simulation
		double m0;	// mass / DEM particle (kg) 
		double dt;	// time increment
		double eta;	// viscous constant
		double eta2;
		double K1,K2;	// spring constants
		double sig;	// characteristic constants	
		int mvw;	// move water (1: yes, 0:No)
		double Eps;	// potential amplitude (normalized)
		double Eps0;	// potential amplitude [kg (m/s)^2]
		double Xa[2];	// cell lower left corner 
		double Xc[2];	// cell center 
		double Wd[4];	// cell size 
		double Wd0[4];	// initial cell size 
		double vmin,vmax;// initial velocity range 
		double Tset;// Temperature [K]
		double time_T_start,time_T_end; // time range [ps] for temperature control
		int T_cntrl; // Temperature control (0:No, 1:Yes)
		int cntrl;	// 0:strain 1:stress given condition
		double sxx,sxy,syy;	// stress 
		int nave;	// number of datapoints for moving average
		double exx[2], eyy[2];	// normal REV strains
		double exy[2], eyx[2];	// shear REV strains
		double txx[2],tyy[2];	// times to define ramp function
		double txy[2],tyx[2];	// times to define ramp function
		int iprd[2];	// 0: open, 1: periodic B.C. 
		int prd; 	// iprd[0]+iprd[1]
		void load(char *fname);// load data 
		void print();		// print parameters	
		double h0[2];		// subcell size (given)
		//double mu;
		double RH_ex;	// R.H. of water resarvoir
		int itime; // current time
		int Nadd; 
	private:
};
class WALL{	// WALL configuration and movement 
	public: 
		double ux[4], uy[4];	// wall displacement
		double Tu[4];		// period
		double brst[4];		// burst cycle
		void setup( );		// import wall displacement data 
		void disp(CNTRL prms, int it);
		double x[4],y[4];	// wall position
		double x0[4],y0[4];	// wall position
		int iwl[4];		// 1:enabled, 0:disabled
		int nwall;		// number of wall particles
	private:
};
class REV{	// Representative Elementary Volume 
	public: 
	double Xa[2]; 	// lower left corner (current)
	double Xa0[2]; 	// lower left corner (initial)
	//double Xc[2];	// center (fixed)
	double Wd[4];	// REV size   (current)
	double Wd0[4];	// REV size   (initial)
	double *w0,*w1,*w2,*w3; // REV size (time series)
	int nwt;	// number of time steps for w0 - w3
	int cntrl;	// 0:strain 1:stress given condition
	double sxx,sxy,syy;	// stress 
	double kxx,kyy,kxy;	// acceleration factor
	double exx[2],eyy[2];	// normal strain 
	double exy[2],eyx[2];	// normal strain 
	double txx[2]; // time defining the ramp function
	double tyy[2]; // time defining the ramp function
	double txy[2]; // time defining the ramp function
	double tyx[2]; // time defining the ramp function
	int iprd[2];	// periodic boundary (1:yes, 0:no)
	Vec2 av,bv; // basis vectors
	Vec2 ah,bh; // basis vectors (dual)
	double Vol; // volume (area)
	void setup(CNTRL prms, WALL wll);	// setup REV 
	void set_lattice_vec();
	void update(int it, double dt);
	void update(int it, double dt, WALL wll);
	Vec2 unfold(double *x, int *ofst);
	void print();
	double h0[2];		// subcell size (given)
	double hd[2];		// subcell size (adjusted)
	int Nh[2];		// numbe of subcells
	int nsub;		// total number of subcells 
	double *bfxx,*bfyy,*bfxy, *bfT, *bfU; // buffer for moving averaging  
	double sxxb,syyb,sxyb,Tb,Ub;
	double sxxb_pre,syyb_pre,sxyb_pre,Tb_pre,Ub_pre;
	double dWxx,dWyy,dWxy;
	int ibf;		// current buffer address
	int nave; 		// number of data points for moving average
	void smooth(double S[2][2], double Tk, double UE); // moving average
	private:
};
class SUBCELL{
	public:
		int indx[2]; // 2D index
		int np;		// number of particles
		int np_max; //max number of particles
		int *list; // particle number set
		void setup(REV rev, double sig);
		~SUBCELL();
	private:
};
class PRTCL{
	public:
		int mobile;	// 1 for mobile, 0 for immobile particle
		double x[2];	// position
		double nx[2];	// normal vector
		double v[2];	// velocity
		double F[2];	// force vector
		double m;	// mass
		int irev[2];	// cell index 
		double sig;	// characteristic distance
		double sigs[2]; // characteristic distances (top,bottom)
		
		double nH2O[2];	// number of hydrated H2O molecules 
		
		double UE[2];	// Pontential Energy (due associated with top/bottom side)
		int sbcll;	// subcell No.
		void init();	// initialize class instance
		void setX(double x1, double x2); // set position vector
		void setX(double x1, double x2, REV rev, CNTRL prms); // set position vector
		void incV(CNTRL prms);
		void incX(CNTRL prms, REV rev);
		void scale(int it, double dt, REV rev);
		void x2xh(REV rev);
		void xh2x(REV rev);
		int ist;	// sheet No.
		int ipt;	// particle No.
		double KE(double m0);
		double get_sig(double xf[2],int *sgn_xf,double dsig[2]);
		double get_sig_iso(double xf[2],int *sgn_xf,double dsig[2]);
	private:
 };


class SHEET{
	public:
		int Np;	// number of particles 
		double K1,K2;	// spring constants
		double r1,r2;	// natural lengths
		double sig_w,sig_s; // characteristic lengths
		void init(int N); // initialization
		int *list;	// particle numbers
		void set_prms(
			double k1, 
			double k2, 
			PRTCL*PT, 
			double sig
		);
		void set_prms(
			CNTRL prms, PRTCL *PT
		);
		int wall;	// is_wall ? --> 1:yes, 0:no		
		double xa[2],xb[2],xc[2],rd; // bounding box;
		double x0[0];	//centroid
		void bbox(PRTCL *PT, CNTRL prms);
		double v1,v2;	// constant velocity vector  
		double incN(PRTCL *PTC, REV rev, double Sab[2][2]);	// increment axial force
		double incQ(PRTCL *PTC, REV rev, double Sab[2][2]);	// increment bending force
		void set_vel(PRTCL *PTC, double vx, double vy);		// set velocity vector 
		Curve2D crv;
		void xy2crv(REV rev,PRTCL *PTC);
		void wsmooth(REV rev, PRTCL *PTC, CLAY NaMt);
	private:
}; 

void Set_Vel(PRTCL *PT, SHEET *st,CNTRL prms);


double LJ(
	PRTCL p1, 	// particle 1
	PRTCL p2, 	// particle 2
	double *Fn, 	// Interaction force vector (out) 
	double sig, 	// characteristic length
	double Eps,	// potential amplitude
	double *Wd,	// cell (REV) size 
	int ncut	// cutoff distance (rcut3=sig*ncut)
);

double LJ(
	PRTCL p1, 	// particle 1
	PRTCL p2, 	// particle 2
	double *Fn, 	// interaction force vector (out)
	double sig, 	// characteristic length
	double Eps,	// potential amplitude
	double *Wd,	// cell (REV) size
	int *ofst1,	// 2D cell index for particle 1 
	int *ofst2,	// 2D cell index for particle 2
	double Sab[2][2] // mechanical stress 
);

double STF1(	// Force due to contact spring  
	PRTCL p1, 	// particle 1 
	PRTCL p2, 	// particle 2 (* two particles belong to the same sheet)
	double *Fn,	// interaction force velctor (out)
	double r0, 	// natural length of the contact spring
	double K1, 	// spring constant 
	//double *Wd	// cell (REV) size
	REV rev,	// REV parameters
	double Sab[2][2]// mechanical stress
);

double VanDerWaals( // Van der Waals force
	SHEET stj,	// sheet 1
	SHEET stk,	// sheet 2
	PRTCL *PTC, 	// particle set 
	CNTRL prms	// DEM control parameters 
);
double VanDerWaals(
	PRTCL *PTC, 	// particle dataset
	CNTRL prms	// DEM control parameters
);

int dist_sheet(
	SHEET st1, 	// sheet 1
	SHEET st2, 	// sheet 2
	double sig	// characteristic length
);
int dist_sheet2(
	SHEET st1, 	// sheet 1
	SHEET st2, 	// sheet 2
	double sig, 	// characteristic length
	int *irev1, 	// 2D cell index for sheet 1 
	int *irev2, 	// 2D cell index for sheet 2
	double *Wd	// cell size
);
void print_ptc(
	PRTCL *PTC,	// particle class array 
	int np,		// number of particles
	int isum,	// file number
	double dt,	// time increment 
	REV rev,	// unit cell
	double rho_d,	// dry density
	double pr,	// porosity 
	char *Dir	// Output data directory
);
void regist_ptc(
	REV rev,	// REV class
	PRTCL *PTC,	// particle set
	SUBCELL *sbcll,	// subcells 
	int np		// number of particles
);
double VDF(	// Van der Waals Force (cell-based evaluation)
	REV rev, 	// REV
	SUBCELL *sbcll, // Subcell
	PRTCL *PTC, 	// Particles
	int iprd[2],	// periodic B.C. 
	double sig,	// characteristic length
	double Eps,	// potential amplitude
	double Sab[2][2] // mechanical stress
);
double VarUE(
	REV rev, 
	SUBCELL *sbcll, 
	PRTCL *PTC, 
	int iprd[2],
	double sig, double Eps,
	int ipts[2],
	int isds[2],
	double var_sig,
	double dUE[2]
);
double VarUE_mu(
	REV rev, 
	SUBCELL *sbcll, 
	PRTCL *PTC, 
	int iprd[2],
	double sig, double Eps,
	int ipt,	// perturbed particles 
	int isde, 	// side (0:top, 1: bottom)
	double var_sig	// amount of water taken/supplied
);
double VDF_L(
	REV rev, 	// REV
	PRTCL *PTC, 	// Particles
	int iprd[2],	// periodic B.C.
	double Eps, 	// pontential amplitude
	double Sab[2][2],// mechanicals stress
	int np, 	// number of particles
	int npl, 	// number of large particles
	int *indx	// index set of large particles
);
void save_ptc_data(
	int isum,	// file number
	double tt, 	// time
	double rho_d, 	// dry density  
	double pr,	// porosity
	REV rev,	// REV data
	int np,		// number of particles
	PRTCL *PTC,	// particle data
	int nst,	// number of sheets
	SHEET *st,	// sheet
	char *Dir	// output data directory
);
void restart(
	CNTRL *prms,
	int np,	// number of particles
	PRTCL *PTC,
	int nst,
	SHEET *st
);
double Uhyd_sig(double U0, double sig);

