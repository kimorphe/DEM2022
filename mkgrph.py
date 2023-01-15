import numpy as np
import matplotlib.pyplot as plt

# Heaviside unit step function H(x-a)
def Hx(x,a=0):
    jndx=np.array(x>a)
    y=np.zeros(len(x))
    y[jndx]=1.0
    return(y)

class Hx_gss:
    def __init__(self,sig):
        N=201
        x=np.linspace(-6*sig,6*sig,N)
        arg=x/sig
        arg*=arg
        y=np.exp(-arg*0.5)/np.sqrt(np.pi*2.)/sig
        dx=x[1]-x[0];
        Y=np.cumsum(y)*dx
        self.y=Y
        self.x=x

    def show(self,ax):
        ax.plot(self.x,self.y)

    def get_val(self,xval):
        indx=np.argmin(np.abs(self.x-xval))
        return(self.y[indx])

    def eval(self,xcod):
        y=np.zeros(len(xcod));
        k=0
        for x in xcod:
            y[k]=self.get_val(x)
            k+=1
        return(y)


if __name__=="__main__":


    fname="data.txt"    # <--inkscape object geometry
    fp=open(fname,"r")
    print("Name, x, y, width, height")  # data format

    #   Set Scale
    txt=fp.readline(); print(txt)
    data=txt.strip().split(",")
    Xa=[float(data[1]),float(data[2])] # reference point 1

    txt=fp.readline(); print(txt)
    data=txt.strip().split(",")
    Xb=[float(data[1]),float(data[2])] # reference point 2
    H=Xb[0]-Xa[0] # horizontal length unit 
    V=Xb[1]-Xa[1] # vertical length unit
    Xunit=0.2 # [-] scale facgtor( relative humidity )
    Yunit=2.0 # [A] scale factor ( basal spacing )

    # ------ Set Origin --------
    txt=fp.readline()
    data=txt.strip().split(","); print(data)
    x0=float(data[1])   
    y0=float(data[2])
    # Physical quantities at the origin
    X0=0.0 # [-] (relative humidity)
    Y0=9.0 # [A] (basal spacing)

    # Read data points
    xcod=[]; ycod=[]
    for row in fp:
        data=row.strip().split(","); print(data)
        xcod.append(float(data[1]))
        ycod.append(float(data[2]))
    fp.close()
    # convertion to the physical quantities
    xcod=(np.array(xcod)-x0)/H*Xunit+X0     # relative humidity
    ycod=(np.array(ycod)-y0)/V*(Yunit)+Y0   # basal spacing  

    #------------- Plotting ----------------

    rj=np.array([0.31,0.72])
    rH=xcod;
    # basal spacings
    h0=10.0; # 0th order swelling
    h1=12.4; # 1st order swelling 
    h2=15.2; # 2nd order swelling
    y=np.zeros(len(rH))
    z=np.zeros(len(rH))
    stp1=Hx_gss(0.05)    # Error function class instance (1st jump)
    stp2=Hx_gss(0.02)    # Error function class instance (2nd jump)
    stp=[stp1,stp2]

    y+=h0;
    z+=h0;
    dh=[h1-h0,h2-h1] # jumps
    
    k=0
    for ri in rj:
        y=y+Hx(rH,ri)*dh[k] # perfect staircase function
        z=z+stp[k].eval(rH-ri)*dh[k] # smoothed staircase function
        k+=1
    #--------------------------------------

    plt.rcParams["font.size"]=12
    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    ax.grid(True)
    ax.set_xlim([-5,105.0])
    ax.set_ylim([9,16.])
    
    ax.plot(xcod*100,ycod,"ko-",markersize=8,label="measured")
    ax.set_xlabel("relative humidity [%]")
    ax.set_ylabel("basal spacing [A]")
    ax.set_title("Na-Mt at 50[deg], Morodome & Kawamura (2009)",fontsize=12)
    fig1.savefig("swelling.png",bbox_inches="tight")
    
    fig2=plt.figure()
    bx1=fig2.add_subplot(211)
    bx1.plot(xcod*100,ycod,"ko-",markersize=8,label="measured")
    bx1.plot(rH*100,y,"-o",label="fitted ( step func.)")
    bx1.plot(rH*100,z,"-s",label="fitted (error func.)")
    bx1.grid(True)
    bx1.set_xlim([-5,105.0])
    bx1.set_ylim([9,16.])
    bx1.set_ylabel("basal spacing [A]")

    bx2=fig2.add_subplot(212)
    bx2.plot(rH*100,ycod-y,"-",label="step func.")
    bx2.plot(rH*100,ycod-z,"-o",label="error func.")
    bx2.grid(True)
    bx2.set_xlim([-5,105.0])
    bx2.set_ylim([-2,2])
    bx2.set_xlabel("relative humidity [%]")
    bx2.set_ylabel("misfit [A]")

    dz0=ycod-z;

    win=np.ones(3)/3
    dz=np.convolve(dz0,win,mode="same")
    dz[0]*=1.5;
    dz[-1]*=1.5;
    bx2.plot(rH*100,dz,"-",linewidth=2,label="smoothed")
    z_fit=dz+z;
    bx1.plot(rH*100,z_fit,"-",linewidth=2,label="smoothed")

    bx1.legend()
    bx2.legend()

    plt.tight_layout()
    fig3=plt.figure(figsize=[4,4])
    cx=fig3.add_subplot(111)
    cx.grid(True)
    #cx.plot(z_fit,rH*100)
    #cx.plot(z,np.log(rH)*1.00)
    cx.plot(ycod,np.log(rH)*1.00,"-ko",markersize=6)
    cx.set_xlim([9.5,16])
    #cx.set_ylim([-5,105])
    cx.set_title("log inverse swelling curve")
    cx.set_xlabel("basal spacing $h$ [A]")
    cx.set_ylabel("log $\gamma_H(h) $")

    fig2.savefig("fitted.png",bbox_inches="tight")
    fig3.savefig("log_rH.png",bbox_inches="tight")

    txt="# relative humidity, basal spacing[A](measured, stair-case fitted)\n"
    for k in range(len(rH)):
        txt+=str(rH[k])+", "+str(ycod[k])+", "+str(z[k])+"\n"
    fp=open("rH_x.dat","w")
    fp.write(txt)
    fp.close()

    plt.show()

