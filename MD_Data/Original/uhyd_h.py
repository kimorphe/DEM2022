import numpy as np
import matplotlib.pyplot as plt


#   Smoothing by moving averaging 
def smooth(ydat,nwd,extend=True):

    ndat=len(ydat)  # numbe of data points
    if extend:
        y=odd_extension(ydat)   # odd-extend given function
    else:
        y=np.zeros(len(ydat))   
        y+=ydat; # copy data

    n=len(y)
    nsum=np.ones(n)
    I=np.ones(n)
    Y=np.zeros(n)
    Y+=y
    nwd2=int(nwd/2)
    for k in range(nwd2):
        m=k+1 
        Y[0:n-m]+=y[m:n]
        Y[m:n]+=y[0:n-m]
        nsum[0:n-m]+=I[m:n]
        nsum[m:n]+=I[0:n-m]
    Y/=nsum;

    if extend:
        return(Y[ndat:ndat*2]) # slice and return smoothed data
    else:
        return(Y)

# Make odd-extension of y
def odd_extension(y):
    n=len(y)
    indx=np.array(range(n))
    indx=indx.astype(int)
    yL=2*y[0]-y[n-indx-1]
    yR=2*y[n-1]-y[n-indx-1]
    y_ext=np.hstack([ yL,y[0:n],yR])

    return(y_ext)
#  Numerical Differentioation by central differencing
def cdiff(y):    
    n=len(y)
    Y=np.zeros(n)
    dy=np.diff(y)
    Y[0:n-1]=dy;
    Y[1:n]+=dy
    Y[1:-1]*=0.5;
    return(Y)



#-------------------------------------------------------
class DATA:
    def __init__(self,fname):
        self.fname=fname
        self.nH2O=[]
        self.Uhyd=[]
        self.b_spc=[]
        self.cod_num=[]
        self.prss=[]
    def load(self):
        fp=open(self.fname,"r")
        fp.readline()

        nH2O=self.nH2O
        Uhyd=self.Uhyd
        b_spc=self.b_spc
        cod_num=self.cod_num
        prss=self.prss
        for row in fp:
            dat=row.strip().split(",")
            nH2O.append(float(dat[0]))
            Uhyd.append(float(dat[1]))
            b_spc.append(float(dat[2]))
            cod_num.append(float(dat[3]))
            prss.append(float(dat[4]))
        fp.close()

        self.ndat=len(nH2O)     # number of data points
        self.nH2O=np.array(nH2O)# n(H2O) n in chemical formula
        self.Uhyd=np.array(Uhyd) # hydration energy 
        self.b_spc=np.array(b_spc) # basal spacing
        self.cod_num=np.array(cod_num) # coordination number 
        self.prss=np.array(prss) # pressure increment

    def who(self):
        print("Data File Name: "+self.fname)
    def print_data(self):
        print(" Trying to print data...")
        if len(self.nH2O) ==0:
            print(" --> data to show has not been loaded !\n")
        else:
            for k in range(self.ndat):
                print(self.nH2O[k], self.Uhyd[k], self.b_spc[k], self.cod_num[k], self.prss[k])
            print("Number of lines=",self.ndat)
    def plot_logU(self,ax,name=""):
        x=np.log10(np.abs(self.nH2O))
        y=np.log10(np.abs(self.Uhyd))

        n1=13
        n2=-1
        P=np.polyfit(x[n1:n2],y[n1:n2],1)
        px=np.poly1d(P)
        Y=px(x[n1:n2])
        print("Whole plot fitting--> ",P)
        ax.plot(x,y,"o",label=name)
        ax.plot(x[n1:n2],Y,"-")
        ax.set_xlabel( "log$_{10}$n")
        ax.set_ylabel( "Energy: log$_{10}$E")

        n1=0
        n2=13
        P=np.polyfit(x[n1:n2],y[n1:n2],1)
        px=np.poly1d(P)
        Y=px(x[n1:n2])
        print("Partial plot fitting--> ",P)
        ax.plot(x[n1:n2],Y,"-")
    def plot_U(self,ax):
        ax.plot(self.nH2O,self.Uhyd)
        ax.set_ylabel("$U_{hyd}$")
    def plot_b_spc(self,ax):
        ax.plot(self.nH2O,self.b_spc)
        ax.set_ylabel("basal spacing")
    def plot_cod_num(self,ax):
        ax.plot(self.nH2O,self.cod_num)
        ax.set_ylabel("coordination number")
    def plot_prss(self,ax):
        ax.plot(self.nH2O,self.prss)
        ax.set_ylabel("pressure increment")
    def plot_dev(self,ax,m=1.0,name=""):
        x=self.nH2O
        y=self.Uhyd

        P=np.polyfit(x,y,1)
        #P[1]=0
        px=np.poly1d(P)
        Y=px(x)

        ax.plot(x,(y-Y)/self.nH2O,"-",markersize=8,label=name+"(line fitting)")

        #m=1.1
        Dn=np.sum(x**(2*m));
        Nm=np.sum((x**m)*y)
        A=Nm/Dn
        print("A=",A)
        print("P[0]=",P[0])
        dev=(y-A*x**m)/self.nH2O
        #ax.plot(x,(y-A*x**m)/self.nH2O,label=name+"(curve fitting)")
        ax.plot(x,dev,label=name+"(curve fitting)")

        ax.set_title("Energy Fluctuation (monotonic trend subtracted)")
        return(dev)

#-------------------------------------------------------

if __name__=="__main__":

    Ca=DATA("CaMt.dat") 
    Ca.load() # load Ca-Montmorillonite data
    Na=DATA("NaMt.dat")
    Na.load() # load Na-Montmorillonite data


#       LogLog Energy Plot
    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    ax.grid(True)
    ax.set_aspect(1.0)
    print("Na-Mt")
    Na.plot_logU(ax,name="Na")    
    print("Ca-Mt")
    Ca.plot_logU(ax,name="Ca")
    ax.legend()


#       Deviatric components
    fig3=plt.figure()
    cx=fig3.add_subplot(111)
    #Na_dev=Na.plot_dev(cx,m=1.09,name="Na") 
    #Ca_dev=Ca.plot_dev(cx,m=1.0,name="Ca")
    #Na_dev=Na.plot_dev(cx,m=1.15,name="Na") 
    #Ca_dev=Ca.plot_dev(cx,m=1.04,name="Ca")
    Na_dev=Na.plot_dev(cx,m=1.097,name="Na") 
    Ca_dev=Ca.plot_dev(cx,m=0.963,name="Ca")
    cx.grid(True)
    cx.legend()

#   ---------------------------------------------
#      Variation in terms of h (basal spacing)
    fig4=plt.figure()
    ex1=fig4.add_subplot(141)
    ex2=fig4.add_subplot(142)
    ex3=fig4.add_subplot(143)
    ex4=fig4.add_subplot(144)
    ex=[ex1,ex2,ex3,ex4]

    ex1.plot(Ca.nH2O,Ca.b_spc,"or-",label="Ca")
    ex1.plot(Na.nH2O,Na.b_spc,"ob-",label="Na")
    ex1.legend()

    nsmp=5
    dn=Na.nH2O[1]-Na.nH2O[0];
    y1=smooth(Na.b_spc,nsmp,extend=True)
    y2=smooth(Ca.b_spc,nsmp,extend=True)
    Na_dn=dn/cdiff(y1)
    Ca_dn=dn/cdiff(y2)
    ex2.plot(Ca_dn,Ca.b_spc,"or-",label="Ca")
    ex2.plot(Na_dn,Na.b_spc,"ob-",label="Na")

    Na_mu=cdiff(smooth(Na_dev*Na.nH2O,nsmp,extend=True))/dn
    Ca_mu=cdiff(smooth(Ca_dev*Ca.nH2O,nsmp,extend=True))/dn
    ex3.plot(Ca_mu,Ca.b_spc,"ro-")
    ex3.plot(Na_mu,Na.b_spc,"bo-")

    Na_dC=cdiff(smooth(Na.cod_num,nsmp,extend=True))/dn
    Ca_dC=cdiff(smooth(Ca.cod_num,nsmp,extend=True))/dn
    ex4.plot(Ca_dC,Ca.b_spc,"ro-")
    ex4.plot(Na_dC,Na.b_spc,"bo-")


    fsz=12
    for exj in ex:
        exj.grid(True)
        exj.set_ylim([10.5,21.5])
        exj.tick_params(labelsize=fsz)

        xlim=exj.get_xlim()
        exj.hlines(12.4,xlim[0],xlim[1],linestyles="dashed")#,linewidth=2)
        exj.hlines(15.6,xlim[0],xlim[1],linestyles="dashed")#,linewidth=2)
        exj.hlines(19.0,xlim[0],xlim[1],linestyles="dashed")#,linewidth=2)
        exj.set_xlim(xlim)

    ex1.set_ylabel("Basal Spacing h",fontsize=12)
    ex1.set_xlabel("n(H$_2$O)",fontsize=12)
    ex2.set_xlabel("$dn/dh$",fontsize=12)
    ex3.set_xlabel("$dU/dn$",fontsize=12)
    ex4.set_xlabel("$dN_c/dn$",fontsize=12)

    ex2.plot(Ca_dn,Ca.b_spc,"or-",label="Ca")
    ex2.plot(Na_dn,Na.b_spc,"ob-",label="Na")


    txt="# basal spacing (Na, Ca), increment dn of H2O molecules(Na, Ca)\n"
    for k in range(len(Ca_dn)):
        txt+=str(Na.b_spc[k])+", "+str(Ca.b_spc[k])+", "+str(Na_dn[k])+", "+str(Ca_dn[k])+"\n"

    fp=open("dndh.dat","w")
    fp.write(txt)
    fp.close()
    plt.show()


