
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

		x=[]; sigp=[]; nwp=[];
		y=[]; sigm=[]; nwm=[];
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
			nwp.append(float(dat[8]));
			nwm.append(float(dat[9]));
		

		self.x=x;
		self.y=y; #print("x=",x); input("pause")
		self.irev=irev;
		self.jrev=jrev;
		self.sigp=np.array(sigp)*10 # [A]
		self.sigm=np.array(sigm)*10 # [A]
		self.nwp=np.array(nwp)
		self.nwm=np.array(nwm)
		self.Np=Np
		fp.close();
	def hist(self,ax,bins=30):
            nw=np.concatenate([self.nwp,self.nwm])
            ax.hist(nw,bins=bins,color="g")
            ax.grid(True)

	def plot_w(self,ax):
		lw=0.5
		ax.plot( self.nwp,"-k",linewidth=lw)
		ax.plot( self.nwm,"-k",linewidth=lw)
		ax.grid(True)
	def plot_w2(self,ax,nps):
		clrs=["r","b","g","c","y","m","k"];
		nclrs=len(clrs);
		n1=0;
		st=0;
		for n in nps:
                    n2=n1+n;
                    nwp_seg=np.array(self.nwp[n1:n2])
                    nwm_seg=np.array(self.nwm[n1:n2])
                    num=np.arange(n1,n2,1)
                    n1=n2
                    cl=clrs[st%nclrs]
                    ax.plot(num,nwp_seg,cl,linewidth=1.0)
                    ax.plot(num,nwm_seg,"--"+cl,linewidth=1.0)
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

    plt.rcParams["font.size"]=12
    fig=plt.figure(figsize=[10,3])
    ax=fig.add_subplot(111)

    fp=open("ptc_nums.dat");
    nums=fp.readlines();
    nums=list(map(int,nums));

    num=0
    args=sys.argv;
    narg=len(args)
    if narg >1:
        num=int(args[1])
    fname="x"+str(num)+".dat"

    ptc=PTCS(fname);

    ptc.plot_w(ax)
    ptc.plot_w2(ax,nums)
    ax.set_title(fname)
    ax.set_xlabel("particle No.")
    ax.set_ylabel("n(H2O)")

    ax.set_xlim([0,ptc.Np])
    ax.set_ylim([0,4])
    #ax.tick_params(labelsize=14)

    fig2=plt.figure(figsize=[6,4])
    ax2=fig2.add_subplot(111)
    ptc.hist(ax2,bins=40)
    ax2.set_title(fname)
    ax2.set_ylabel("count")
    ax2.set_xlabel("n(H2O)")
    ax2.set_xlim([0,4])
    plt.show()
