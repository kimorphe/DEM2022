
#! /home/kazushi/Enthought/Canopy_64bit/User/bin/python

import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------
fp=open("energy.out","r");
fp.readline()
time=[]
K=[]; P=[];
T=[];
W=[]; H=[];
Uh=[]; Un=[];
for row in fp:
	try:
		dat=list(map(float, row.strip().split(" ")));
		T.append(dat[5]);
		H.append(dat[4]);
		W.append(dat[3]);
		P.append(dat[2]);
		K.append(dat[1]);
		Uh.append(dat[7]);
		Un.append(dat[8]);
		time.append(dat[0]);
	except:
		break;

#---------------------------------------
fp=open("stress.out","r");
fp.readline()
tau=[]
Sxx=[]; Sxy=[]; Syy=[]; 
for row in fp:
	try:
		dat=list(map(float, row.strip().split(" ")));
		Syy.append(dat[10])
		Sxy.append(dat[9])
		Sxx.append(dat[8])
		tau.append(dat[0])
	except:
		break;


tau=np.array(tau);
Sxx=np.array(Sxx)*1.e03;
Sxy=np.array(Sxy)*1.e03;
Syy=np.array(Syy)*1.e03;
#---------------------------------------



K=np.array(K);
P=np.array(P);
fig=plt.figure()
ax1=fig.add_subplot(311)
ax2=fig.add_subplot(312)
ax3=fig.add_subplot(313)
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)

ax1.plot(time,K,label="kinetic")
ax1.plot(time,P,label="potential");
ax1.plot(time,P+K,label="total");
ax2.plot(time,T,label="temperature");

ax3.plot(tau,Sxx,label="s_{xx}");
ax3.plot(tau,Sxy,label="s_{xy}");
ax3.plot(tau,Syy,label="s_{xy}");

ax3.set_xlabel("time [ps]")
ax1.set_ylabel("Energy [kJ/mol]")
ax2.set_ylabel("Temperature [K]")
ax3.set_ylabel("Stress [MPa]")

fig2=plt.figure()
ax2=fig2.add_subplot(111)
ax2.plot(time,P,label="LJ")
ax2.plot(time,Uh,label="Hyd")
ax2.plot(time,Un,label="G")
ax2.plot(time,P+Uh+Un,label="total")
ax2.legend()
ax2.grid(True)
ax2.set_xlabel("time");
ax2.set_ylabel("Energy");
plt.show()
