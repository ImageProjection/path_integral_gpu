import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace
from scipy.fftpack import fft,fftshift,fftfreq
from scipy.fftpack.basic import ifft
from scipy.stats import norm

fig=plt.figure()
ax1=fig.add_subplot(2,1,1)
f=open("out_dens_plot.txt",'r')
x_data=[]
y_data=[]
NFFT=0
for line in f:
    NFFT+=1#number of points for IFFT to process
    tmp_list=list(map(float,line.split(",")))
    x_data.append(tmp_list[0])
    y_data.append(tmp_list[1])


#for testing
'''it works for this gauss
y_data=np.zeros(NFFT)
for i in range(250):
    y_data[i]=norm.pdf((i-0)/20)
'''
ax1.set_ylabel("|P(x)|^2")
ax1.set_xlabel("coordinate x")
ax1.plot(x_data,y_data)
ax1.grid(color = 'black', linestyle = '--', linewidth = 0.5)


ax2=fig.add_subplot(2,1,2)
ax2.set_xlabel("momentum p (units not to scale to anything)")
ift_ydata=ifft(y_data)
ift_xdata=np.linspace(-1,1,NFFT,endpoint=False)
#normalising
ax2.plot(ift_xdata[512-60:512+60],fftshift(ift_ydata)[512-60:512+60])
ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
plt.locator_params(nbins=20)
plt.tight_layout()
plt.savefig("histogram.png",bbox_inches='tight',dpi=350)
plt.show()