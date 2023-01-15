import numpy as np
import matplotlib.pyplot as plt


x=np.linspace(0,1,101)

sig=0.025
xb=0.25
beta=0.5

arg=(x-xb)/sig
arg*=arg
y=beta*np.exp(-0.5*arg)/np.sqrt(2.*np.pi)/sig

xb=0.75
arg=(x-xb)/sig
arg*=arg
y+=beta*np.exp(-0.5*arg)/np.sqrt(2.*np.pi)/sig


fig=plt.figure()
ax=fig.add_subplot(111)
dx=x[1]-x[0]
Y=np.cumsum(y)*dx
#ax.plot(x,y)
ax.plot(x,Y)
ax.grid(True)
plt.show()
