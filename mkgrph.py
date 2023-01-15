import numpy as np
import matplotlib.pyplot as plt

def Hx(x,a=0):  # Heaviside unit step function H(x-a)
    jndx=np.array(x>a)
    y=np.zeros(len(x))
    y[jndx]=1.0
    return(y)


fname="data.txt"
fp=open(fname,"r")


print("Name, x, y, width, height")


#       Scale 
txt=fp.readline()
print(txt)
data=txt.strip().split(",")
Xa=[float(data[1]),float(data[2])]

txt=fp.readline()
print(txt)
data=txt.strip().split(",")
Xb=[float(data[1]),float(data[2])]

Xunit=0.2 #[-] relative humidity 
Yunit=2 # [A] basal spacing
H=Xb[0]-Xa[0] # horizontal axis
V=Xb[1]-Xa[1] # vertical axis

# ------ Origin --------
txt=fp.readline()
data=txt.strip().split(",")
x0=float(data[1])
y0=float(data[2])
print(data)
# Physical origin 
X0=0.0 #[-] 
Y0=9.0 #[A]


xcod=[]
ycod=[]
for row in fp:
    data=row.strip().split(",")
    print(data)
    xcod.append(float(data[1]))
    ycod.append(float(data[2]))
fp.close()


xcod=(np.array(xcod)-x0)/H*Xunit+X0
ycod=(np.array(ycod)-y0)/V*(Yunit)+Y0

#--------------------------------------
#rj=np.array([0.3,0.7])
rj=np.array([0.3,0.725])
rH=xcod;
h0=10.0;
h1=12.4
h2=15.2
y=np.zeros(len(rH))
y+=h0;

dh=[h1-h0,h2-h1]
k=0
for ri in rj:
    y=y+Hx(rH,ri)*dh[k]
    k+=1
#--------------------------------------

plt.rcParams["font.size"]=14
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(xcod*100,ycod,"ko-",markersize=8)
ax.grid(True)
ax.set_xlabel("relative humidity [%]")
ax.set_ylabel("basal spacing [A]")
ax.set_title("Na-Mt at 50[deg], (Morodome & Kawamura (2009)",fontsize=12)
ax.plot(rH*100,y,"-o")
ax.plot(rH*100,ycod-y,"-s")
fig.savefig("swelling.png",bbox_inches="tight")
plt.show()

