import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace # library for plotting
from scipy.fftpack import fft,fftshift,fftfreq


def sine_wave(f,overSampRate,phase,nCyl):
	"""
	Generate sine wave signal with the following parameters
	Parameters:
		f : frequency of sine wave in Hertz
		overSampRate : oversampling rate (integer)
		phase : desired phase shift in radians
		nCyl : number of cycles of sine wave to generate
	Returns:
		(t,g) : time base (t) and the signal g(t) as tuple
	Example:
		f=10; overSampRate=30;
		phase = 1/3*np.pi;nCyl = 5;
		(t,g) = sine_wave(f,overSampRate,phase,nCyl)
	"""
	fs = overSampRate*f # sampling frequency
	t = np.arange(0,nCyl*1/f-1/fs,1/fs) # time base
	g = np.sin(2*np.pi*f*t+phase) + np.cos(2*np.pi*f*5*t) # replace with cos if a cosine wave is desired
	return (t,g) # return time base and signal g(t) as tuple



f = 10 #frequency = 10 Hz
overSampRate = 30 #oversammpling rate
fs = f*overSampRate #sampling frequency
phase = 1/3*np.pi #phase shift in radians
nCyl = 5 # desired number of cycles of the sine wave

(t,x) = sine_wave(f,overSampRate,phase,nCyl) #function call

plt.plot(t,x) # plot using pyplot library from matplotlib package
plt.title('Sine wave f='+str(f)+' Hz') # plot title
plt.xlabel('Time (s)') # x-axis label
plt.ylabel('Amplitude') # y-axis label
#plt.show() # display the figure

NFFT=1024 #NFFT-point DFT      
X=fft(x,NFFT) #compute DFT using FFT    

fig1, ax = plt.subplots(nrows=1, ncols=1) #create figure handle
#fVals_ar=1/NFFT*fs*np.arange(start = -NFFT/2,stop = NFFT/2) # raw index for FFT plot
#print(fVals_ar)
fVals_ft=fftshift(fftfreq(NFFT,1/fs))
print(fVals_ft)
ax.plot((fVals_ft),fftshift(np.abs(X)))      
ax.set_title('Double Sided FFT - withFFTShift')
ax.set_xlabel('freq in Hz')        
ax.set_ylabel('DFT Values')
plt.locator_params(nbins=40)
plt.show()