import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace
from scipy.fftpack import fft,fftshift,fftfreq
from scipy.fftpack.basic import ifft
from scipy.stats import norm




f = 5 #frequency = 10 Hz
overSampRate = 30 #oversammpling rate
fs = f*overSampRate #sampling frequency


NFFT=1024 #NFFT-point DFT      

fVals=np.zeros(NFFT)
for i in range(NFFT):
    fVals[i]=norm.pdf((i-NFFT/2)/NFFT*10)
#t=np.linspace(0,NFFT*1/fs,NFFT,endpoint=False)
t=np.linspace(0,1,NFFT,endpoint=False)

x=ifft(fftshift(fVals),NFFT)

plt.plot(t,x) # plot using pyplot library from matplotlib package
plt.title('IFFT') # plot title
plt.xlabel('Time (s)') # x-axis label
plt.ylabel('sig value') # y-axis label
plt.show() # display the figure