#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class SHEET:
	def __init__(self,npnt):
		self.x=[];
		self.y=[];
		self.npnt=npnt;
		self.irev=[];
		self.jrev=[];
	def plot(self,ax,clr):
		#xx=self.x+Wd[0]*i;
		#yy=self.y+Wd[0]*j;
		indx= abs(np.diff(self.irev)) >0   
		indx=np.append(indx,True);
		id=np.where(indx == True)
		id=np.append([0],id);
		print(id) 
		print(np.shape(id));
		#indx = indx.tolist();
		ax.plot(self.x,self.y,"-"+clr,lw=2);

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

		x=[];
		y=[];
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
		

		self.x=x;
		self.y=y; #print("x=",x); input("pause")
		self.irev=irev;
		self.jrev=jrev;
		fp.close();

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
				plt,=ax.plot(x[m1:m2],y[m1:m2],"-"+clrs[st%nclrs],ms=2,lw=1);
				plts.append(plt);
				m1=m2;
			n1=n2;
			st+=1;

		return plts;

#------------------------ MAIN ROUTINE -------------------------------

if __name__=="__main__":

	nfile1=80;
	nfile2=250; inc=10;

	args=sys.argv
	narg=len(args)
	if  narg > 1:
		nfile1=int(args[1]);
	if  narg > 2:
		nfile2=int(args[2]);
	if  narg > 3:
		inc=int(args[3]);

	#imv=int( raw_input("Movie(0) or Snapshots(1) ?"));
	imv=int( input("Movie(0) or Snapshots(1) ?"));
	if imv==0:
		MV=True;
	else:
		MV=False;

	fp=open("ptc_nums.dat");
	nums=fp.readlines();
	nums=list(map(int,nums));

	fig=plt.figure();
	ax=fig.add_subplot(111);
	ax.set_aspect(1.0)
	ax.grid(True);

	#fnum=np.linspace(0,100,101);
	fnum=np.arange(nfile1,nfile2+1,inc);
	fnum=fnum.astype(int)

	if MV:
		ims=[];
		Writer=animation.writers['ffmpeg']	
		writer=Writer(fps=10,metadata=dict(artist="Me"),bitrate=1800)

	for k in fnum:
		fname="x"+str(k)+".dat";
		ptc=PTCS(fname);
		plts=ptc.plot(ax,nums,Movie=MV);
		#ax.set_xlim([0,200]);
		#ax.set_ylim([0,200]);
		ax.set_xlim([0,75]);
		ax.set_ylim([0,75]);
		ax.set_xlabel("x [nm]",fontsize="12")
		ax.set_ylabel("y [nm]",fontsize="12")
		#title=ax.text(200,0,fname+", t="+str(ptc.time)+"[ps]",ha="right")
		title=ax.text(100,0,fname+", t="+str(ptc.time)+"[ps]",ha="right")
		plts.append(title);
		ax.grid(True)

		if MV:
			ims.append(plts);
		else:
			plt.savefig("x"+str(k)+".png",bbox_inches="tight");
			#plt.pause(0.1);
	plt.savefig("x"+str(k)+".png",bbox_inches="tight");

	if MV:
		fn="movie.mp4"
		ani=animation.ArtistAnimation(fig,ims,interval=1,repeat_delay=100)
		print(" Wait, creating movie file ...");
		ani.save(fn,writer=writer)
		print(" Done.")
		print(" Play "+fn+" (e.g. by typing 'totem "+fn+"')")
