import numpy as np
import matplotlib.pyplot as plt


class MEAS:
    def __init__(self,fname="rH_x.dat"):
        fp=open(fname,"r")
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
    """
    def rH_val(self,hz):
        if hz<self.hz[0]:
            return(self.rH[0])
        if hz>=self.hz[-1]:
            return(self.rH[-1])

        indx=np.argmin(np.abs(self.hz-hz))

        h0=self.hz[indx]
        dh=hz-h0
        if dh>0:
            i1=indx
            i2=i1+1
            xi=dh/(self.hz[i2]-self.hz[i1])
            eta=1.-xi
            return(self.rH[i1]*eta+self.rH[i2]*xi)
        else:
            i2=indx
            i1=i2-1
            eta=-dh/(self.hz[i2]-self.hz[i1])
            xi=1.-eta
            return(self.rH[i1]*eta+self.rH[i2]*xi)
    """


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


if __name__=="__main__":
    
    plt.rcParams["font.size"]=12


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

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(hs, dndh,"-ko",label="dn/dh")
    ax.plot(hs, np.log(rH+1.e-08),"-o",label="log(rH)")

    mu_sat=-1.0 # chemical potential of saturated H2O vapor 
    dUdh=(mu_sat+np.log(rH+1.e-08))*dndh # Gradiendt of U_hyd (hydration energy/particle)
    ax.plot(hs, dUdh,"-o",label="dU/dh(h)")
    ax.grid(True)
    ax.set_xlabel("basal spacing [$\AA$]")
    ax.legend()
    ax.set_ylim([-3.5,3])

    Uh=np.cumsum(dUdh)*(hs[1]-hs[0])  #hydration energy
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.plot(hs,Uh,"k-",linewidth=2)
    ax2.grid(True)
    ax2.set_xlabel("basal spacing [$\AA$]")
    ax2.set_ylabel("hyndration energy/particle/$\mu_{sat}$")

    txt="# basal spacing [A], hydration energy\n"
    for k in range(len(hs)):
        txt+=str(hs[k])+", "+str(Uh[k])+"\n"

    fp=open("uhyd.dat","w")
    fp.write(txt)
    fp.close()

    fig2.savefig("uhyd.png",bbox_inches="tight")

    plt.show()

