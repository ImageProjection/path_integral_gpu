'''
input: p and x trajectories files and gen_params_file
output: -p trajectories density plot
        -x trajectories density plot
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import time

#grab simulaton parameters
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
traj_x_range=values[13]
sigma=values[14]
print_termo_traj_flag=values[15]

N_bins=1024
discarded_p_points=0
discarded_x_points=0

#h_hist holds counts for bins
def add_hist(p_traj, h_hist, p_range):
        bin_width=2.0*p_range/N_bins
        discarded_points=0
        for i in range(N_spots):
                if abs(p_traj[i]) < p_range:
                        bin_i=math.floor((p_traj[i]-(-p_range))/bin_width)
                        h_hist[bin_i]+=1
                else:
                        discarded_points+=1
        return discarded_points

    
fig=plt.figure()
h_p_hist=np.zeros(N_bins)
#main for p
ax1=fig.add_subplot(1,1,1)
fp=open("out_p_traj.txt",'r')
#skip to sampling trajectories part in both files
if print_termo_traj_flag:
    for i in range(N_waiting_trajectories):
        fp.readline()
#read each line add to cumulative histogram
for line in fp:
        full_line=list(map(float,line.split(",")))
        p_traj=full_line[0:N_spots]
        discarded_p_points+=add_hist(p_traj,h_p_hist, p_range)

h_p_dense_plot=h_p_hist/(np.sum(h_p_hist)*2.0*p_range/N_bins)


x_data=np.linspace(-p_range,p_range,N_bins)
y_data=h_p_dense_plot
ax1.grid(color = 'black', linestyle = '--', linewidth = 0.5)
ax1.set_ylabel("|P(p)|^2")
ax1.set_xlabel("координата p")
ax1.set_xlim(-p_range,p_range)
ax1.set_xticks(np.arange(-p_range,p_range,1.0))
ax1.step(x_data,y_data,color='red',linewidth = 1.1)
