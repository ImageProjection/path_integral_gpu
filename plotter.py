import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from pathlib import Path

#simulaton parameters
values=[]
spf=open("out_gen_des.txt",'r')
for line in spf:
    tmp_list=list(line.split(","))
    values.append(float(tmp_list[1]))

N_spots=round(values[0])
N_sweeps_waiting=round(values[1])
N_sample_trajectories=round(values[2])
Traj_sample_period=round(values[3])
a=values[4]
beta=values[5]
v_fermi=values[6]
m=values[7]
omega=values[8]
p_bottom=values[9]
p_range=values[10]
x_range=values[11]



fig=plt.figure()
ax1=fig.add_subplot(1,1,1)

def init_p():
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-p_range,p_range])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_xlabel("")
    line,=ax1.plot([],[],color="blue")
    return line,

def init_x():
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-x_range,x_range])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_xlabel("")
    line,=ax1.plot([],[],color="blue")
    return line,

def upd_p(frame_i):
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

def upd_x(frame_i):
    ax1.clear()
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-x_range,x_range])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_xlabel("")
    ax1.grid(color = 'black', linestyle = '--', linewidth = 0.5)
    sigma=dig_ar[frame_i][N_spots]
    my_xlabel=f"sample traj No={frame_i}\n"
    my_xlabel+=f"sigma={sigma}\n"
    line,=ax1.plot(range(1,N_spots+1),dig_ar[frame_i][0:N_spots],color="blue",lw=0.8)
    ax1.set_xlabel(my_xlabel,loc='left')
    return line,

#main for p
start_time=time.time()
dig_ar=[]
f=open("out_p_traj.txt",'r')
n_lines=0
for line in f:
    dig_ar.append(list(map(float,line.split())))
    n_lines+=1

dig_ar=dig_ar[n_lines-50:]
ani=animation.FuncAnimation(fig, upd_p, init_func=init_p, interval=200,frames=50, repeat=False, blit=0)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
#plt.show()
writervideo = animation.FFMpegWriter(fps=4)
Path("traj_mp4").mkdir(exist_ok=True)
ani.save('traj_mp4/p_traj_evolution.mp4', writer=writervideo)
f.close()
end_time=time.time()
print("elapsed time plotting p (seconds):",round(end_time-start_time,1))


#main for x
start_time=time.time()
dig_ar=[]
f=open("out_x_traj.txt",'r')
n_lines=0
for line in f:
    dig_ar.append(list(map(float,line.split())))
    n_lines+=1

dig_ar=dig_ar[n_lines-50:]
ani=animation.FuncAnimation(fig, upd_x, init_func=init_x, interval=200,frames=50, repeat=False, blit=0)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
#plt.show()
writervideo = animation.FFMpegWriter(fps=4)
Path("traj_mp4").mkdir(exist_ok=True)
ani.save('traj_mp4/x_traj_evolution.mp4', writer=writervideo)
f.close()
end_time=time.time()
print("elapsed time plotting x (seconds):",round(end_time-start_time,1))

