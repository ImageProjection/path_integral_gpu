'''
For input: all files produced during worksession (including raw output) and gen_params files
out: puts results of launch into ../pi_results/<descriptive-folder-name> 
'''
import os
from datetime import datetime

now = datetime.now()
folder_name=now.strftime("%d/%m/%Y %H:%M:%S")
os.mkdir(folder_name)
os.system()
