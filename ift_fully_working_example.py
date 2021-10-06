import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace
from scipy.fftpack import fft,fftshift,fftfreq
from scipy.fftpack.basic import ifft
from scipy.stats import norm





NFFT=1024 #NFFT-point DFT      
fVals_half=np.zeros(512)#amplitudes for each frequency
fVals=np.zeros(NFFT)
for i in range(512):
    fVals_half[i]=np.exp(-(i/20)**2)*(5*(i/20)**3+5*(i/20)**2+0.9)
#fVals=np.concatenate((fVals_half,np.zeros(512)))
fVals=np.concatenate((fVals_half,np.flip(fVals_half)))

#t=np.linspace(0,NFFT*1/fs,NFFT,endpoint=False)
t=np.linspace(0,1,NFFT,endpoint=False)

x=ifft(fVals,NFFT)

fig=plt.figure()
ax1=fig.add_subplot(311)
ax2=fig.add_subplot(312)
ax3=fig.add_subplot(313)

ax1.plot(fVals)
ax1.set_title('spectrum without fftshift')

ax2.plot(fftshift(fVals))
ax2.set_title('spectrum with fftshift')

ax3.plot(t,fftshift(x))
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('signal value')
plt.show()
