import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N_spots=1024

fig=plt.figure()
ax1=fig.add_subplot(2,1,1)

def init():
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-4,4])
    ax1.set_xticks(ticks=list(range(0,N_spots,round(N_spots/8)))+[N_spots])
    line,=ax1.plot([],[],color="blue")
    return line,

def upd(frame_i):
    line,=ax1.plot(range(1,N_spots+1),dig_ar[frame_i],color="blue")
    return line,

dig_ar=[]
f=open("out_traj.txt",'r')
n_lines=0
for line in f:
    dig_ar.append(list(map(float,line.split())))
    n_lines+=1

ani=animation.FuncAnimation(fig, upd, init_func=init, interval=450,frames=n_lines, repeat=False, blit=1, save_count=50)
plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
plt.show()
f.close()

