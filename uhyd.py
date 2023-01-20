import numpy as np
import matplotlib.pyplot as plt


class MEAS:
    def __init__(self,fname="rH_x.dat"):
        fp=open(fname,"r")
        ndat=int(fp.readline())
        print("ndat=",ndat)
        fp.readline()

        rH=[]
        hz=[]
        for row in fp:
            dat=row.strip().split(",")
            rH.append(float(dat[0]))
            hz.append(float(dat[1]))
        fp.close()

        self.rH=np.array(rH)
        self.hz=np.array(hz)
    def plot(self,ax):
        ax.plot(self.rH*100,self.hz,"ko-",markersize=8)
        ax.set_xlabel("relative humidity [%]")
        ax.set_ylabel("basal spacing [$\AA$]")
        ax.grid(True)
    def get_seg_ID(self,hz):
        if hz<self.hz[0]:
            return(-1)
        if hz>=self.hz[-1]:
            return(len(self.hz))

        indx=np.argmin(np.abs(self.hz-hz))
        h0=self.hz[indx]
        dh=hz-h0
        if dh>0:
            return(indx)
        else:
            return(indx-1)

    def rH_val(self,hz):
        i1=self.get_seg_ID(hz);
        if i1<0:
            return(self.rH[0])
        if i1>= len(self.hz):
            return(self.rH[-1])
        i2=i1+1
        xi=(hz-self.hz[i1])/(self.hz[i2]-self.hz[i1])
        eta=1.-xi
        return(self.rH[i1]*eta+self.rH[i2]*xi)

    def hz_val(self,rval):
        rH=self.rH
        hz=self.hz
        if rval<=rH[0]:
            return(hz[0])
        if rval>=rH[-1]:
            return(hz[-1])

        indx=np.argmin(np.abs(rH-rval))
        r0=self.rH[indx]
        dr=rval-r0
        if dr>0:
            i1=indx
        else:
            i1=indx-1
        i2=i1+1

        xi=(rval-rH[i1])/(rH[i2]-rH[i1])
        eta=1.-xi
        return(hz[i1]*eta+hz[i2]*xi)


class SIM:
    def __init__(self,fname):
        fp=open(fname,"r")
        fp.readline()
        hz=[]
        nH2O=[]
        dndh=[]
        for row in fp:
            data=row.strip().split(",")
            hz.append(float(data[0]))
            nH2O.append(float(data[2]))
            dndh.append(float(data[4]))

        self.hz=np.array(hz)
        self.nH2O=np.array(nH2O)
        self.dndh=np.array(dndh)

    def plot_nH2O(self,ax):
        ax.plot(self.hz,self.nH2O,"ko-",markersize=8)
        ax.set_xlabel("basal spacing [$\AA$]")
        ax.set_ylabel("n(H2O)")
        ax.grid(True)
    def plot_dndh(self,ax):
        ax.plot(self.hz,self.dndh,"ko-",markersize=8)
        ax.set_xlabel("basal spacing [$\AA$]")
        ax.set_ylabel("dn/dh[/$\AA$]")
        ax.grid(True)

    def get_seg_ID(self,hz):
        if hz<self.hz[0]:
            return(-1)
        if hz>=self.hz[-1]:
            return(len(self.hz))

        indx=np.argmin(np.abs(self.hz-hz))
        h0=self.hz[indx]
        dh=hz-h0
        if dh>0:
            return(indx)
        else:
            return(indx-1)

    def nH2O_val(self,hz):
        i1=self.get_seg_ID(hz);
        if i1<0:
            return(self.nH2O[0])
        if i1>= len(self.hz):
            return(self.nH2O[-1])
        i2=i1+1
        xi=(hz-self.hz[i1])/(self.hz[i2]-self.hz[i1])
        eta=1.-xi
        return(self.nH2O[i1]*eta+self.nH2O[i2]*xi)

    def dndh_val(self,hz):
        i1=self.get_seg_ID(hz);
        if i1<0:
            return(self.dndh[0])
        if i1>= len(self.hz):
            return(self.dndh[-1])
        i2=i1+1
        xi=(hz-self.hz[i1])/(self.hz[i2]-self.hz[i1])
        eta=1.-xi
        return(self.dndh[i1]*eta+self.dndh[i2]*xi)

    def hz_val(self,nw):
        dn=self.nH2O[1]-self.nH2O[0]
        i1=int((nw-self.nH2O[0])/dn)
        if i1<0:
            return(self.hz[0]);
        if i1>=len(self.nH2O):
            return(self.hz[-1]);

        i2=i1+1
        xi=(nw-self.nH2O[i1])/dn
        eta=1.-xi

        return(eta*self.hz[i1]+xi*self.hz[i2]);
if __name__=="__main__":
    
    plt.rcParams["font.size"]=12


    #*********** Functions of h (basal spacing) ****************
    h_xrd=MEAS(fname="rH_x.dat")    # load XRD-measured swelling curve 
    h_sim=SIM("MD_Data/mk_dndh.dat")# load MD-generated swelling curve

    hs=np.linspace(h_xrd.hz[-1],h_xrd.hz[0],201) # make h(basal spacing) axis
    dndh=[] # dn(H2O)/dh
    rH=[]   # relative humidity [0,1]
    for h in hs:    # linearly interplated data
        dndh.append(h_sim.dndh_val(h))
        rH.append(h_xrd.rH_val(h))
    dndh=np.array(dndh)
    rH=np.array(rH)
    kB=1.380649*1.e-23 #[J/K] Boltzmann constant
    NA=6.02214076*1.e23 #[/mol] Avogadro's number
    R_gas=kB*NA #[J/K/mol] Gas Constant
    T=300   # [K]
    mu_sat=-44.5*1.e03    #[J/mol]
    beta=R_gas*T/abs(mu_sat) # RT/mu_sat (R:gas constant, T:temperature, mu_sat:chemical potential at saturation)
    dUdh=(-1.0+beta*np.log(rH+1.e-08))*dndh # Gradiendt of U_hyd (hydration energy/particle)
    mu=-1+beta*np.log(rH+1.e-08)

    fig=plt.figure(figsize=[5,6])
    ax1=fig.add_subplot(311)
    ax2=fig.add_subplot(312)
    ax3=fig.add_subplot(313)
    ax1.plot(hs, dndh,"-ko",label="dn/dh")
    ax2.plot(hs,mu,"-o",label="mu/mu_sat")
    ax3.plot(hs, dUdh,"-o",label="dU/dh(h)")
    axs=[ax1,ax2,ax3]
    for ax in axs:
        ax.grid(True)
        ax.legend()
    ax2.set_ylim([-1.2,-0.9])
    ax3.set_xlabel("basal spacing $h$ [$\AA$]")


    #*********** Functions of n(H2O) (H2O molar number) ****************
    #   rH --(XRD)--> hz --(MD)--> n(H2O) 
    r1=0.0
    h=h_xrd.hz_val(r1)
    n1=h_sim.nH2O_val(h)

    r2=1.0
    h=h_xrd.hz_val(r2)
    n2=h_sim.nH2O_val(h)
    print("n(H2O) limits=",n1,n2)

    Nr=101;
    nws=np.linspace(n2,n1,Nr);   # make n(H2O) axis
    mu_n=np.zeros(Nr)
    hz_n=np.zeros(Nr)
    rh_n=np.zeros(Nr)
    k=0
    for nw in nws:
        hn=h_sim.hz_val(nw) # n(H2O) --(MD)--> hz
        rn=h_xrd.rH_val(hn) # hz --(XRD)--> rH
        mn=-1+beta*np.log(rn+1.e-08) # rH --> mu (chem. potential)
        hz_n[k]=hn
        rh_n[k]=rn
        mu_n[k]=mn
        k+=1
    mu_n+=1;        # subtract mu_sat (chem.potential at saturation)
    mu_n[-1]*=0.5   # for trapezoidal integration 
    Un=np.cumsum(mu_n)*(nws[1]-nws[0]) # integrate mu_n w.r.t. n(H2O)

    fig2=plt.figure(figsize=(5,6))
    bx1=fig2.add_subplot(311)
    bx2=fig2.add_subplot(312)
    bx3=fig2.add_subplot(313)
    bxs=[bx1,bx2,bx3]
    bx1.plot(nws,hz_n,label="h(n)")
    bx2.plot(nws,rh_n,label="R.H.(n)")
    bx3.plot(nws,Un,label="u_hyd(n)")
    for bx in bxs:
        bx.grid(True)
        bx.legend()
        bx3.set_xlabel("n(H2O)")
    bx1.set_ylabel("basal spacing [$\AA$]")
    bx2.set_ylabel("R.H.")
    bx3.set_ylabel("hydration energy")

    fig2.savefig("uhyd.png",bbox_inches="tight")
    plt.show()

    """
    nH2O=[]
    for hz in hs:
        nH2O.append(h_sim.nH2O_val(hz))
    nH2O=np.array(nH2O)

    Uh=np.cumsum(dUdh)*(hs[1]-hs[0])  #hydration energy
    fig3=plt.figure()
    bx=fig3.add_subplot(111)
    bx.plot(nH2O,Uh,"k-",linewidth=2)
    bx.grid(True)
    bx.set_xlabel("n(H2O)[mol/unit lattice]")
    bx.set_ylabel("hyndration energy/particle/$\mu_{sat}$")

    txt="# basal spacing [nm], hydration energy\n"
    ndat=len(hs)
    txt+=str(ndat)+"\n"
    for k in range(ndat):
        txt+=str(hs[ndat-1-k]*0.1)+", "+str(Uh[ndat-1-k])+"\n"

    fp=open("uhyd.dat","w")
    fp.write(txt)
    fp.close()
    plt.show()
    """

