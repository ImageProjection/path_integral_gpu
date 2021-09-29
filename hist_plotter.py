import matplotlib
import numpy as np
import matplotlib.pyplot as plt

fig=plt.figure()
ax1=fig.add_subplot(2,1,1)
f=open("out_dens_plot.txt",'r')
x_data=[]
y_data=[]
NFFT=0
for line in f:
    NFFT+=1#number of points for FFT to process
    tmp_list=list(map(float,line.split(",")))
    x_data.append(tmp_list[0])
    y_data.append(tmp_list[1])

ax1.set_ylabel("|P(x)|^2")
ax1.set_xlabel("coordinate x")
ax1.plot(x_data,y_data)
ax1.grid(color = 'black', linestyle = '--', linewidth = 0.5)

ax2=fig.add_subplot(2,1,2)
sample_spacing=x_data[1]-x_data[0]
ft_ydata=
ft_xdata=
ax2.plot()
ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
plt.locator_params(nbins=40)
plt.show()
plt.savefig("histogram.png")