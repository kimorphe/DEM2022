import numpy as np
import matplotlib.pyplot as plt


def Hx(x,a=0):  # Heaviside unit step function H(x-a)
    jndx=np.array(x>a)
    y=np.zeros(len(x))
    y[jndx]=1.0
    return(y)

if __name__=="__main__":

    rj=np.array([0.3,0.7])

    rH=np.linspace(0,1,101)
    y=np.zeros(len(rH))

    fig=plt.figure()
    ax=fig.add_subplot(111)

    h0=10.0;
    h1=12.4
    h2=15.2

    for ri in rj:
        y=y+Hx(rH,ri)

    ax.plot(rH,y,"-o")
    ax.grid(True)
    plt.show()

