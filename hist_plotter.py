import matplotlib
import numpy as np
import matplotlib.pyplot as plt

fig=plt.figure()
ax1=fig.add_subplot(1,1,1)
f=open("out_dens_plot.txt",'r')
x_data=[]
y_data=[]
for line in f:
    tmp_list=list(map(float,line.split(",")))
    x_data.append(tmp_list[0])
    y_data.append(tmp_list[1])

plt.ylabel("|\psi(p)|^2")
plt.xlabel("coordinate")
plt.plot(x_data,y_data)
plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
plt.show()