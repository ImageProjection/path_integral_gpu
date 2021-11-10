import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace


fig=plt.figure()
#plot momentum density plot
ax1=fig.add_subplot(2,1,1)
fp=open("out_p_dens_plot.txt",'r')
x_data=[]
y_data=[]
for line in fp:
    tmp_list=list(map(float,line.split(",")))
    x_data.append(tmp_list[0])
    y_data.append(tmp_list[1])

ax1.set_ylabel("|P(p)|^2")
ax1.set_xlabel("coordinate p")
ax1.plot(x_data,y_data, color='red',linewidth = 1.1)
ax1.grid(color = 'black', linestyle = '--', linewidth = 0.5)


#plot coordinate density plot
ax2=fig.add_subplot(2,1,2)
fx=open("out_x_dens_plot.txt",'r')
x_data=[]
y_data=[]
for line in fx:
    tmp_list=list(map(float,line.split(",")))
    x_data.append(tmp_list[0])
    y_data.append(tmp_list[1])

ax2.set_ylabel("|P(x)|^2")
ax2.set_xlabel("coordinate x")
ax2.plot(x_data,y_data, color='green',linewidth = 1.1)
ax2.grid(color = 'black', linestyle = '--', linewidth = 0.5)


#display results
plt.locator_params(nbins=20)
plt.tight_layout()
plt.savefig("p_and_x_dens_plots.png",bbox_inches='tight',dpi=350)
#plt.show()