'''
list of all output files to be copied to separate storage folder:
global_averages.txt
local_averages.txt
p_and_x_dens_plots.png
out_gen_des.txt
out_p_traj.txt
p_traj_evolution.mp4
'''


files_list=("global_averages.txt "+
"local_averages.txt "+
"p_and_x_dens_plots.png "+
"out_gen_des.txt "+
"p_traj_evolution.mp4 ")

import sys
from subprocess import Popen
import os
import time
from glob import glob
from time import localtime, strftime
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble']=r'\usepackage[utf8]{inputenc}'
plt.rcParams['text.latex.preamble']=r'\usepackage[russian]{babel}'
#eg ax1.set_xlabel(r'значение параметра $\beta$')

num_cores=4
uniq_id=1
#clean folder before launch
os.system("git clean -fx")
os.system("rm -rf work")
#sys.exit()
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
beta_start=190
beta_stop=520
n_beta_points=5
beta_list=np.linspace(beta_start,beta_stop,n_beta_points,endpoint=True)
n_periods_list=[1.2, 1.4, 1.6, 1.8, 2.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]#p_bottom for now

#create overhead launch folder
date_time=strftime("%d.%m_%H:%M", localtime())
multi_beta_folder_name=(date_time
    +"_beta_"+str(round(beta_start,2))+"to"+str(round(beta_stop,2))
    +"/")
os.system("mkdir ../path_integral_gpu_results/"+multi_beta_folder_name)

core_folders_list=[]
for i in range(num_cores):
    core_folders_list.append("t"+str(i)+"/")
for i in range(0,len(beta_list)):
    os.system("rm -rf work")
    #launch(ie reconfigure cpp file, then do nb_long_run)
    repl("const int N=490;","const int N="+str(int(beta_list[i]))+";")
    #os.system("make nb_long_run") this line is now a long procedure, aimed to create out p traj and then
    #go on with long run
    os.system("make nb_compile")#other part after folder merge
    os.system("mkdir work")
    for j in range(len(core_folders_list)):
        os.system("mkdir work/"+core_folders_list[j])
        os.system("cp a.out work/"+core_folders_list[j])
    os.system("rm a.out")
    
    for j in range(len(core_folders_list)):
        aout_path="/home/artem/Documents/path_integral_gpu/work/"+core_folders_list[j]+"a.out"
        fold_path="/home/artem/Documents/path_integral_gpu/work/"+core_folders_list[j]
        Popen(aout_path,cwd=fold_path)
    
    #sys.exit()
    while(1):
        time.sleep(1)
        fileCounter = 0
        for root, dirs, files in os.walk("work/"):
            for file in files:    
                if file.endswith('dummy.txt'):
                    fileCounter += 1
        if fileCounter>=num_cores:
            break
    os.system("touch out_p_traj.txt")
    for j in range(len(core_folders_list)):
        os.system("work/"+core_folders_list[j]+"out_p_traj.txt "+">> "+"../out_p_traj.txt")            
    os.system("cd ..")    
    #can now go on
    os.system("make nb_long_run_no_c")
    repl("const int N="+str(int(beta_list[i]))+";","const int N=490;")
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
os.system(("cp "+"energy_plots.py " + "../path_integral_gpu_results/"
    +multi_beta_folder_name))

