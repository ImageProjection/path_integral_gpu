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
            -global average E
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



#main
#list of values for each observable, later convert to np.arrays
aver_T_vals=[]#p
aver_V_vals=[]#x
aver_p_dot_vals=[]#p
aver_E_vals=[]#from both
aver_rel_vals=[]#from both
kink_metr_vals=[]#p
fp=open("out_p_traj.txt",'r')
fx=open("out_x_traj.txt",'r')
#getting p-related local data
for line in fp:
    full_line=list(map(float,line.split(",")))
    p_traj=full_line[0:N_spots]

    aver_T=aver_T(p_traj)
    aver_p_dot=aver_p_dot(p_traj)
    kink_metr=kink_metr(p_traj)

    aver_T_vals.append(aver_T)
    aver_p_dot_vals.append(aver_p_dot)
    kink_metr_vals.append(kink_metr)

#getting x related local data
for line in fx:
    full_line=list(map(float,line.split(",")))
    x_traj=full_line[0:N_spots]

    aver_V=aver_V(x_traj)
    
    aver_V_vals.append(aver_V)

#evaluating "both" global values
for i in range(N_sample_trajectories):
    aver_E_vals.append(aver_T_vals[i]-aver_p_dot_vals[i])
    aver_rel_vals.append(aver_T_vals[i]/aver_V_vals[i])


#at this point all lists of values a filled
#and ready for further print and use

aver_T_vals=np.array(aver_T_vals)
aver_V_vals=np.array(aver_V_vals)
aver_p_dot_vals=np.array(aver_p_dot_vals)
aver_E_vals=np.array(aver_E_vals)
aver_rel_vals=np.array(aver_rel_vals)
kink_metr_vals=np.array(kink_metr_vals)

global_aver_T=np.average(aver_T_vals)
global_aver_V=np.average(aver_V_vals)
global_aver_p_dot=np.average(aver_p_dot_vals)
global_aver_E=np.average(aver_E_vals)
global_aver_rel=np.average(aver_rel_vals)
global_aver_kink_metr=np.average(kink_metr_vals)

#evaluating errors
global_aver_T_error=np.std(aver_T_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_V_error=np.std(aver_V_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_p_dot_error=np.std(aver_p_dot_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_E_error=np.std(aver_E_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_rel_error=np.std(aver_rel_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_kink_metr_error=np.std(kink_metr_vals)/math.sqrt(N_sample_trajectories-1)

#format is
#N_sample_traj_index, E, T, V, p_dot, rel, kink_metr
f_locals=open("local_averages.txt","w")

#format is
#E, T, V, p_dot, rel, kink_metr, beta
#same for errors, except error for beta is 0
f_summary=open("global_averages")

for i in range(N_sample_trajectories):
    f_locals.write(str(i)+", "
        +str(aver_E_vals[i])+", "
        +str(aver_T_vals[i])+", "
        +str(aver_V_vals[i])+", "
        +str(aver_p_dot_vals[i])+", "
        +str(aver_rel_vals[i])+", "
        +str(kink_metr_vals[i])+", "
        +str(beta)+"\n")


f_summary.write(str(global_aver_E)+", "
    +str(global_avet_T)+", "
    +str(global_avet_V)+", "
    +str(global_avet_p_dot)+", "
    +str(global_avet_rel)+", "
    +str(global_avet_kink_metr)+"\n")

f_summary.write(str(global_aver_E_error)+", "
    +str(global_avet_T_error)+", "
    +str(global_avet_V_error)+", "
    +str(global_avet_p_dot_error)+", "
    +str(global_avet_rel_error)+", "
    +str(global_avet_kink_metr_error)+", "
    +str(0.0)+"\n")    


f_locals.close()
f_summary.close()
