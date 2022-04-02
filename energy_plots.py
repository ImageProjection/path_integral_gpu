import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble']=r'\usepackage[utf8]{inputenc}'
plt.rcParams['text.latex.preamble']=r'\usepackage[russian]{babel}'


#parse global averages file for E and \beta
f=open("../path_integral_gpu_results/global_averages.txt")

#plot energies
fig=plt.figure()
ax1=fig.add_subplot(1,1,1)
ax1.plot(beta_list,,marker='o')

ax1.set_xlabel(r'значение параметра $\beta$',fontsize='22')
ax1.set_ylabel(r'энергия $E$',fontsize='22')
ax1.set_title(r'зависимость энергии $E$ от обратной температуры $\beta$')