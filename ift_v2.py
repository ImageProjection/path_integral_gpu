import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace
from scipy.fftpack import fft,fftshift,fftfreq
from scipy.fftpack.basic import ifft
from scipy.stats import norm





NFFT=1024 #NFFT-point DFT      
fVals_half=np.zeros(512)
#fVals=np.zeros(NFFT)
for i in range(250):
    fVals_half[i]=norm.pdf((i-125)/25)
fVals=np.concatenate((fVals_half,np.flip(fVals_half)))
'''
fVals[1]=1
fVals[5]=2
fVals[9]=3
fVals[13]=4
fVals[17]=5
fVals[21]=10
fVals[25]=10
fVals[29]=5
fVals[33]=4
fVals[37]=3
fVals[41]=2
fVals[45]=1
'''
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

ax3.plot(t,x)
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('signal value')
plt.show()
