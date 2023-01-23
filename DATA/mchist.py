
#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class PTCS:
	def __init__(self,fname):
		fp=open(fname,"r");
		fp.readline();
		tt=float(fp.readline());
		fp.readline();
		dat=fp.readline();
		dat=dat.split(" ");
		rho_d=float(dat[0]);	
		poro=float(dat[1]);
		print(str(tt)+"[ps], "+str(rho_d)+"[g/cm3]")
		self.time=tt;

		fp.readline();
		dat=fp.readline();
		dat=dat.strip().split(" ");
		Xc=list(map(float,dat));	

		fp.readline();
		dat=fp.readline();
		dat=dat.strip().split(" ");
		Wd=list(map(float,dat));	

		fp.readline();
		Np=int(fp.readline());

		fp.readline();

		x=[]; sigp=[];
		y=[]; sigm=[];
		irev=[];
		jrev=[];
		#for row in fp:
		for k in range(Np):
			dat=fp.readline();
			dat=dat.split(" ");
			x.append(float(dat[2]));	
			y.append(float(dat[3]));	
			irev.append(int(dat[0]));	
			jrev.append(int(dat[1]));	
			sigp.append(float(dat[6]));
			sigm.append(float(dat[7]));
		

		self.x=x;
		self.y=y; #print("x=",x); input("pause")
		self.irev=irev;
		self.jrev=jrev;
		self.sigp=np.array(sigp)*10 # [A]
		self.sigm=np.array(sigm)*10 # [A]
		self.Np=Np
		fp.close();
	def hist(self,ax):
            sh=np.hstack([self.sigp,self.sigm])
            #sh-=0.9;
            #sh*=0.5;
            ax.hist(sh,bins=25)
            ax.grid(True)

	def plot_w(self,ax):
		lw=0.5
		#ax.plot( (self.sigp-0.9)*0.5,"-k",linewidth=lw)
		#ax.plot( (self.sigm-0.9)*0.5,"-k",linewidth=lw)
		ax.plot( self.sigp,"-k",linewidth=lw)
		ax.plot( self.sigm,"-k",linewidth=lw)
		ax.grid(True)
	def plot_w2(self,ax,nps):
		clrs=["r","b","g","c","y","m","k"];
		nclrs=len(clrs);
		n1=0;
		st=0;
		for n in nps:
                    n2=n1+n;
                    sp=np.array(self.sigp[n1:n2])
                    sm=np.array(self.sigm[n1:n2])
                    #sp-=0.9; sm-=0.9
                    #sp*=0.5; sm*=0.5;
                    num=np.arange(n1,n2,1)
                    n1=n2
                    cl=clrs[st%nclrs]
                    ax.plot(num,sp,cl,linewidth=1.0)
                    ax.plot(num,sm,"--"+cl,linewidth=1.0)
                    st+=1
		ax.grid(True)
	def plot(self,ax,nps,Movie=False):
		if Movie == False:
			ax.cla()
	
		clrs=["r","b","g","c","y","m","k"];
		nclrs=len(clrs);
		n1=0;
		st=0;


		plts=[];
		for n in nps:
			n2=n1+n;
			irev=self.irev[n1:n2];
			jrev=self.jrev[n1:n2];
			x=self.x[n1:n2];
			y=self.y[n1:n2];
			itmp=np.abs(np.diff(irev));
			jtmp=np.abs(np.diff(jrev));
			itmp=np.append(itmp,1);
			jtmp=np.append(jtmp,1);

			tmp=itmp+jtmp;
			indx,=np.where(tmp>0)
			indx+=1;

			m1=0;
			for m2 in indx:
				#plt2,=ax.plot(x[m1:m2],y[m1:m2],"-",color="skyblue",ms=2,lw=2);
				plt,=ax.plot(x[m1:m2],y[m1:m2],"-"+clrs[st%nclrs],ms=2,lw=0.5);
				plts.append(plt);
				m1=m2;
			n1=n2;
			st+=1;

		return plts;

if __name__=="__main__":

    fig=plt.figure(figsize=[12,4])
    ax=fig.add_subplot(111)

    fp=open("ptc_nums.dat");
    nums=fp.readlines();
    nums=list(map(int,nums));

    num=4
    args=sys.argv;
    narg=len(args)
    if narg >1:
        num=int(args[1])

    nums=range(num)
    Hz=[]
    for num in nums:
        fname="x"+str(num)+".dat"
        ptc=PTCS(fname);
        #ptc.plot_w(ax)
        hz=np.concatenate([ptc.sigp,ptc.sigm])
        Hz.append(hz)
    Hz=np.array(Hz)
    print(np.shape(Hz))
    vmin=11
    vmax=16
    im=ax.imshow(Hz,aspect="auto",cmap="jet",interpolation="none",origin="lower",vmin=vmin,vmax=vmax)
    plt.colorbar(im)


    print(np.shape(hz))
    #ptc.plot_w2(ax,nums)

    #ax.set_xlim([0,ptc.Np])
    #ax.set_ylim([9,16])
    ax.tick_params(labelsize=14)

    #fig2=plt.figure()
    #ax2=fig2.add_subplot(111)
    #ptc.hist(ax2)
    plt.show()
