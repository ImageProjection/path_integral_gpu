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

#functions
def aver_T_func(p_traj):
    pb=p_bottom
    result=0
    for i in range(N_spots):
        result+= v_fermi*math.sqrt( (p_traj[i]**2-pb**2)**2/(4*pb**2)
             + m*m*v_fermi*v_fermi )
    return result/N_spots

def aver_V_func(x_traj):
    result=0
    for i in range(N_spots):
        result+=0.5*m*m*omega*omega*x_traj[i]*x_traj[i]
    return result/N_spots

def aver_p_dot_func(p_traj):
    result=0
    for i in range(N_spots):
        result+= ( p_traj[i] - p_traj[(i-1+N_spots)%N_spots] )**2
    result= result/(2*a*a*m*omega*omega)
    return result/N_spots

def kink_metr_func(p_traj):
    result=0
    for i in range(N_spots):
        if np.sign(p_traj[i]) != np.sign(p_traj[(i-1+N_spots)%N_spots]):
            result+=1
    return result

#main
#list of values for each observable, later convert to np.arrays
aver_T_vals=[]#p
aver_p_dot_vals=[]#p
aver_E_vals=[]#from both
kink_metr_vals=[]#p

fp=open("out_p_traj.txt",'r')
x=open("out_x_traj.txt",'r')
#skip to sampling trajectories part in both files
#if print_termo_traj_flag:
#    for i in range(N_waiting_trajectories):
#        fx.readline()
#        fp.readline()

#getting p-related local data
for line in fp:
    full_line=list(map(float,line.split()))
    p_traj=full_line[0:N_spots]

    aver_T=aver_T_func(p_traj)
    aver_p_dot=aver_p_dot_func(p_traj)
    kink_metr=kink_metr_func(p_traj)

    aver_T_vals.append(aver_T)
    aver_p_dot_vals.append(aver_p_dot)
    kink_metr_vals.append(kink_metr)

#evaluating "both" global values
for i in range(N_sample_trajectories):
    aver_E_vals.append(aver_T_vals[i]-aver_p_dot_vals[i])


#at this point all lists of values a filled
#and ready for further print and use

aver_T_vals=np.array(aver_T_vals)
aver_p_dot_vals=np.array(aver_p_dot_vals)
aver_E_vals=np.array(aver_E_vals)
kink_metr_vals=np.array(kink_metr_vals)

global_aver_T=np.average(aver_T_vals)
global_aver_p_dot=np.average(aver_p_dot_vals)
global_aver_E=np.average(aver_E_vals)
global_aver_kink_metr=np.average(kink_metr_vals)

#evaluating errors
global_aver_T_error=np.std(aver_T_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_p_dot_error=np.std(aver_p_dot_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_E_error=np.std(aver_E_vals)/math.sqrt(N_sample_trajectories-1)
global_aver_kink_metr_error=np.std(kink_metr_vals)/math.sqrt(N_sample_trajectories-1)

#format is
#N_sample_traj_index, E, T, V, p_dot, rel, kink_metr
f_locals=open("local_averages.txt","w")

#format is
#E, T, V, p_dot, rel, kink_metr, beta, 1/beta
#same for errors, except error for beta is 0
f_summary=open("global_averages.txt","a")

for i in range(N_sample_trajectories):
    f_locals.write(str(i)+", "
        +str(aver_E_vals[i])+", "
        +str(aver_T_vals[i])+", "
        +str(aver_p_dot_vals[i])+", "
        +str(kink_metr_vals[i])+", "
        +str(beta)+"\n")


f_summary.write(str(global_aver_E)+", "
    +str(global_aver_T)+", "
    +str(global_aver_p_dot)+", "
    +str(global_aver_kink_metr)+", "
    +str(beta)+", "
    +str(1/beta)+"\n")

f_summary.write(str(global_aver_E_error)+", "
    +str(global_aver_T_error)+", "
    +str(global_aver_p_dot_error)+", "
    +str(global_aver_kink_metr_error)+", "
    +str(0.0)+", "
    +str(0.0)+"\n")    


f_locals.close()
f_summary.close()
