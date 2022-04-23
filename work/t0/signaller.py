import os
import pathlib
file_path = os.path.split(os.getcwd())[1]
os.system("touch ../"+file_path+".txt")