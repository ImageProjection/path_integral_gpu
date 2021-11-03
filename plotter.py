import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

N_spots=1024
p_range=5

fig=plt.figure()
ax1=fig.add_subplot(1,1,1)

def init():
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-p_range,p_range])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_xlabel("")
    line,=ax1.plot([],[],color="blue")
    return line,

def upd(frame_i):
    ax1.clear()
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-p_range,p_range])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_xlabel("")
    ax1.grid(color = 'black', linestyle = '--', linewidth = 0.5)
    sigma=dig_ar[frame_i][N_spots]
    my_xlabel=f"sample traj No={frame_i}\n"
    my_xlabel+=f"sigma={sigma}\n"
    line,=ax1.plot(range(1,N_spots+1),dig_ar[frame_i][0:N_spots],color="blue",lw=0.8)
    ax1.set_xlabel(my_xlabel,loc='left')
    return line,

#main
start_time=time.time()
dig_ar=[]
f=open("out_p_traj.txt",'r')
n_lines=0
for line in f:
    if (n_lines<50):
        dig_ar.append(list(map(float,line.split())))
    n_lines+=1

ani=animation.FuncAnimation(fig, upd, init_func=init, interval=200,frames=50, repeat=False, blit=0)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
#plt.show()
writervideo = animation.FFMpegWriter(fps=4)
ani.save('p_traj_evolution.mp4', writer=writervideo)
f.close()
end_time=time.time()
print("elapsed time plotting (seconds):",round(end_time-start_time,1))

