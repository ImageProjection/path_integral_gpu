'''
input: trajectories file and gen_params_file
output: -into separate file, each line for 1 sample traj
            -aver_T
            -aver_V (over 1 sample traj)
            -aver_p_dot term
            -aver_E
            -rel
            -number of sign changes ("kink metric")
        -evaluates and puts results into yet another file:
            -global average T
            -global average V
            -global average p_dot term
            -global average E=V-T
            -global average rel
            -global average kink metric
            with statystical errors for each of these quantities
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


#list of values for each observable, later convert to np.arrays
aver_T_vals=[]#p
aver_V_vals=[]#x
aver_p_dot_vals=[]#p
aver_E_vals=[]#from both
aver_rel_vals=[]#from both
kink_metr_vals=[]#p
fp=open("out_p_traj.txt",'r')

#getting p-related local data
for line in fp:
    full_line=list(map(float,line.split(",")))
    p_traj=full_line[0:N_spots]
    aver_T=aver_T(p_traj)
    aver_p_dot=aver_p_dot(p_traj)
    kink_metr=kink_metr(p_traj)

#getting x related local data

#evaluating "both" global values

#at this point all lists of values a filled
#and ready for further use

aver_T_vals=np.array(aver_T_vals)
aver_V_vals=np.array(aver_V_vals)
aver_p_dot_vals=np.array(aver_p_dot_vals)
aver_E_vals=np.array(aver_E_vals)
aver_rel_vals=np.array(aver_rel_vals)
kink_metr_vals=np.array(kink_metr_vals)