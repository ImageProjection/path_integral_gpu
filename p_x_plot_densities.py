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
#main for p
h_hist=np.zeros(N_bins)
ax1=fig.add_subplot(2,1,1)
fp=open("out_p_traj.txt",'r')
#read each line add to cumulative histogram
n_lines=0
for line in fp:
        full_line=list(map(float,line.split(",")))
        p_traj=full_line[0:N_spots]
        discarded_p_points+=add_hist(p_traj,h_hist, p_range)

#h_dense_plot=h_hist/(np.sum(h_hist)*2.0*p_range/N_bins)
h_dense_plot=h_hist

x_data=np.linspace(-p_range,p_range,N_bins)
y_data=h_dense_plot
ax1.set_xlim(-p_range,p_range)
ax1.step(x_data,y_data)
print("n_lines")
plt.show()

