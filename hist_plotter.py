import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace
from pathlib import Path

#simulaton parameters
values=[]
spf=open("out_gen_des.txt",'r')
for line in spf:
    tmp_list=list(line.split(","))
    values.append(float(tmp_list[1]))

N_spots=round(values[0])
N_waiting_trajectories=round(values[1])
N_sample_trajectories=round(values[2])
Traj_sample_period=round(values[3])
a=values[4]
beta=values[5]
v_fermi=values[6]
m=values[7]
omega=values[8]
p_bottom=values[9]
p_range=values[10]
x_range=values[11]
traj_p_range=values[12]
traj_x_range=values[12]



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

Path("traj_hist").mkdir(exist_ok=True)
plt.savefig("traj_hist/m="+str(m)+"_p_and_x_dens_plots.png",bbox_inches='tight',dpi=350)
#plt.show()