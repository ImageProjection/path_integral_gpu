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

def repl(a,b):#a->b
    # Read in the file
    with open('main1.cpp', 'r') as file :
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace(a, b)

    # Write the file out again
    with open('main1.cpp', 'w') as file:
        file.write(filedata)

#main
beta_start=200
beta_stop=7
n_beta_points=5
beta_list=np.linspace(beta_start,beta_stop,n_beta_points,endpoint=True)
n_periods_list=[1.2, 1.4, 1.6, 1.8, 2.0]#p_bottom for now

#create overhead launch folder
date_time=strftime("%d.%m_%H:%M", localtime())
multi_beta_folder_name=(date_time
    +"_beta_"+str(round(beta_start,2))+"to"+str(round(beta_stop,2))
    +"/")
os.system("mkdir ../path_integral_gpu_results/"+multi_beta_folder_name)

for i in range(0,len(beta_list)):
    #launch(ie reconfigure cpp file, then do nb_long_run)
    repl("const int N=000;","const int N="+str(int(beta_lis[i]))+";")
    os.system("make nb_long_run")
    #create folder
    single_beta_folder_name=("N="+str(round(beta_list[i],0))+"_uid"+str(uniq_id)+"/")
    uniq_id+=1
    os.system(("mkdir ../path_integral_gpu_results/"+multi_beta_folder_name
    +single_beta_folder_name))
    #copy results
    os.system(("cp "+files_list+" ../path_integral_gpu_results/"
    +multi_beta_folder_name+single_beta_folder_name))

#copy last termod summary to overhead folder
os.system(("cp "+"global_averages.txt " + "../path_integral_gpu_results/"
    +multi_beta_folder_name))

