import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble']=r'\usepackage[utf8]{inputenc}'
plt.rcParams['text.latex.preamble']=r'\usepackage[russian]{babel}'
plt.rcParams['font.size']=r'20'

energies_list=[]
energy_error_list=[]
beta_list=[]
temperatures_list=[]
#parse global averages file for E and \beta
f=open("../path_integral_gpu/global_averages.txt",'r')
lines=f.readlines()
for i in range(0,len(lines),2):
    val_line=list(map(float,lines[i].split(", ")))
    error_line=list(map(float,lines[i+1].split(", ")))
    energies_list.append(val_line[0])
    energy_error_list.append(error_line[0])
    beta_list.append(val_line[4])
    temperatures_list.append(1/val_line[4])


#plot energies
fig=plt.figure()
fig.subplots_adjust(hspace=.5)
ax1=fig.add_subplot(2,1,1)
ax1.set_xlabel(r'значение параметра $\beta$',fontsize='20')
ax1.set_ylabel(r'энергия $E$ (условные единицы)',fontsize='20')
ax1.grid()
ax1.plot(beta_list,energies_list,marker='o')

ax2=fig.add_subplot(2,1,2)
ax2.set_xlabel(r'значение температуры (условные единицы)',fontsize='20')
ax2.set_ylabel(r'энергия $E$ (условные единицы)',fontsize='20')
ax2.grid()
ax2.plot(temperatures_list,energies_list,marker='o')
plt.show()