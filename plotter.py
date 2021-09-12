import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

N_spots=1024

fig=plt.figure()
ax1=fig.add_subplot(2,1,1)

def init():
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-4,4])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_xlabel("")
    line,=ax1.plot([],[],color="blue")
    return line,

def upd(frame_i,sigma,acc_rate):
    ax1.clear()
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-4,4])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_xlabel("")
    ax1.grid(color = 'black', linestyle = '--', linewidth = 0.5)
    my_xlabel=f"sample traj No={frame_i}\n"
    my_xlabel+=f"sigma={sigma*frame_i}\n"
    my_xlabel+=f"acc_rate={acc_rate}\n"
    line,=ax1.plot(range(1,N_spots+1),dig_ar[frame_i][0:N_spots],color="blue")
    ax1.set_xlabel(my_xlabel,loc='left')
    return line,

#main
start_time=time.time()
dig_ar=[]
f=open("out_traj.txt",'r')
n_lines=0
for line in f:
    dig_ar.append(list(map(float,line.split())))
    n_lines+=1

sigma=5.5
acc_rate=10.222
ani=animation.FuncAnimation(fig, upd, init_func=init,fargs=(sigma,acc_rate), interval=200,frames=n_lines, repeat=False, blit=0)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
plt.show()
#writervideo = animation.FFMpegWriter(fps=30)
#ani.save('traj_evolution.mp4', writer=writervideo)
f.close()
end_time=time.time()
print("elapsed time plotting (seconds):",round(end_time-start_time,1))

