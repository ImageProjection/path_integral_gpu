'''
input: log.txt file
output: -evaluates and puts results into yet another file:
            -global average E
            -global average T
            -global average p_dot term
            with statystical errors for each of these quantities
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import time

import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble']=r'\usepackage[utf8]{inputenc}'
plt.rcParams['text.latex.preamble']=r'\usepackage[russian]{babel}'

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
traj_p_range=values[11]




#main
#list of values for each observable, later convert to np.arrays
aver_T_vals=[]#p
aver_p_dot_vals=[]#p
aver_E_vals=[]#from both

flog=open("out_log.txt",'r')
for line in flog:
    full_line=list(map(float,line.split(',')))
    aver_E_vals.append(full_line[0])
    aver_T_vals.append(full_line[1])
    aver_p_dot_vals.append(full_line[2])

#at this point all lists of values a filled
#and ready for further print and use

aver_T_vals=np.array(aver_T_vals)
aver_p_dot_vals=np.array(aver_p_dot_vals)
aver_E_vals=np.array(aver_E_vals)

global_aver_T=np.average(aver_T_vals)
global_aver_p_dot=np.average(aver_p_dot_vals)
global_aver_E=np.average(aver_E_vals)

#evaluating errors
global_aver_T_error=np.std(aver_T_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_p_dot_error=np.std(aver_p_dot_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_E_error=np.std(aver_E_vals)/math.sqrt(N_sample_trajectories-1)

#format is
#N_sample_traj_index, E, T, V, p_dot, rel, kink_metr
f_locals=open("local_averages_log.txt","w")

#format is
#E, T, p_dot, kink_metr, beta, 1/beta
#same for errors, except error for beta is 0
f_summary=open("global_averages_log.txt","a")

for i in range(N_sample_trajectories):
    f_locals.write(str(i)+", "
        +str(aver_E_vals[i])+", "
        +str(aver_T_vals[i])+", "
        +str(aver_p_dot_vals[i])+", "
        +str(beta)+"\n")


f_summary.write(str(global_aver_E)+", "
    +str(global_aver_T)+", "
    +str(global_aver_p_dot)+", "
    +str(beta)+", "
    +str(1/beta)+"\n")

f_summary.write(str(global_aver_E_error)+", "
    +str(global_aver_T_error)+", "
    +str(global_aver_p_dot_error)+", "
    +str(0.0)+", "
    +str(0.0)+"\n")    


f_locals.close()
f_summary.close()
