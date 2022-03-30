import os
from time import localtime, strftime

date_time=strftime("%d.%m_%H:%M", localtime())
os.system("git clean -fx")
multi_beta_folder_name="many_betas_"+date_time
os.system("mkdir ../path_integral_gpu_results/"+multi_beta_folder_name)