'''
list of all output files to be copied to separate storage folder:
global_averages.txt
local_averages.txt
p_and_x_dens_plots.png
out_gen_des.txt
out_p_traj.txt
out_x_traj.txt
p_traj_evolution.mp4
x_traj_evolution.mp4
'''


files_list=("global_averages.txt "+
"local_averages.txt "+
"p_and_x_dens_plots.png "+
"out_gen_des.txt "+
"out_p_traj.txt "+
"out_x_traj.txt "+
"p_traj_evolution.mp4 "+
"x_traj_evolution.mp4")


import os
from time import localtime, strftime
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble']=r'\usepackage[utf8]{inputenc}'
plt.rcParams['text.latex.preamble']=r'\usepackage[russian]{babel}'
#eg ax1.set_xlabel(r'значение параметра $\beta$')

uniq_id=1
#clean folder before launch
os.system("git clean -fx")

beta_start=9
beta_stop=16
n_beta_points=8
beta_list=np.linspace(beta_start,beta_stop,n_beta_points,endpoint=True)
n_periods_list=[25,25,25,25,25,25,25,25]
#launch for first point
os.system("make nb_long_run beta_val="+str(beta_start)+" n_periods="+str(n_periods_list[0]))

#grab simulation parameters from out_gen_des.txt
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

#create overhead launch folder
date_time=strftime("%d.%m_%H:%M", localtime())
multi_beta_folder_name=(date_time
    +"_beta_"+str(round(beta_start,2))+"to"+str(round(beta_stop,2))
    +"_vf="+str(round(v_fermi,2))
    +"_m="+str(round(m,1))
    +"_w="+str(round(omega,2))
    +"_pb="+str(round(p_bottom,1))
    +"_N_wait="+str(N_waiting_trajectories)
    +"_N_sample"+str(N_sample_trajectories)
    +"/")
os.system("mkdir ../path_integral_gpu_results/"+multi_beta_folder_name)

#create local results folder  for first beta value, put results into it
single_beta_folder_name=("beta="+str(round(beta_start,2))+"a="
    +str(round(beta/N_spots,2))+"_uid"+str(uniq_id)+"/")
uniq_id+=1
os.system(("mkdir ../path_integral_gpu_results/"+multi_beta_folder_name
    +single_beta_folder_name))
os.system(("cp "+files_list+" ../path_integral_gpu_results/"
    +multi_beta_folder_name+single_beta_folder_name))

#for each of the rest beta values
#launch
#put results into folder
for i in range(1,len(beta_list)):
    #launch
    os.system("make nb_long_run beta_val="+str(beta_list[i])+" n_periods="+str(n_periods_list[i]))
    #create folder
    single_beta_folder_name=("beta="+str(round(beta_list[i],2))+"a="
    +str(round(beta/N_spots,2))+"_uid"+str(uniq_id)+"/")
    uniq_id+=1
    os.system(("mkdir ../path_integral_gpu_results/"+multi_beta_folder_name
    +single_beta_folder_name))
    #copy results
    os.system(("cp "+files_list+" ../path_integral_gpu_results/"
    +multi_beta_folder_name+single_beta_folder_name))

#copy last termod summary to overhead folder
os.system(("cp "+"global_averages.txt " + "../path_integral_gpu_results/"
    +multi_beta_folder_name))

