'''
input: p and x trajectories files and gen_params_file
output: -p trajectories mp4 plot
        -x trajectories mp4 plot
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

#grab simulaton parameters
values=[]
spf=open("out_gen_des.txt",'r')
for line in spf:
    tmp_list=list(line.split(","))
    values.append(float(tmp_list[1]))

N_spots=round(values[0])
N_waiting_trajectories=round(values[1])
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
traj_p_range=values[12]
traj_x_range=values[13]

#amount of frames used for video
N_vid_fr=N_waiting_trajectories+N_sample_trajectories#min(300,N_sample_trajectories);
#N_vid_fr=min(200,N_sample_trajectories);

#figure, first used for p, then for x
fig=plt.figure()
ax1=fig.add_subplot(1,1,1)

#init and update function for p and x

def init_p():
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-traj_p_range,traj_p_range])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_xlabel("p_traj")
    line,=ax1.plot([],[],color="blue")
    return line,

def upd_p(frame_i):
    ax1.clear()
    ax1.set_xlim([1,N_spots+1])
    ax1.set_ylim([-traj_p_range,traj_p_range])
    ax1.set_xticks(ticks=list(range(0,N_spots,N_spots//8))+[N_spots])
    ax1.set_yticks(list(plt.yticks()[0]) + [p_bottom])
    ax1.set_xlabel("")
    ax1.grid(color = 'black', linestyle = '--', linewidth = 0.5)
    acc_rate=dig_ar[frame_i][N_spots]
    my_xlabel="номер траектории:"+str(frame_i)+"\n"
    line,=ax1.plot(range(1,N_spots+1),dig_ar[frame_i][0:N_spots],color="blue",lw=0.8)
    ax1.set_xlabel(my_xlabel)
    return line,



#main for p: read data, put animation into mp4 file
start_time=time.time()
dig_ar=[] #2d array, first index is trajectory number, second is node in trajectory
f=open("out_p_traj.txt",'r')
n_lines=0
for line in f:
    dig_ar.append(list(map(float,line.split(","))))
    n_lines+=1

dig_ar=dig_ar[ (n_lines-N_vid_fr if n_lines-N_vid_fr>0 else 0) :]#so that only last are used, if not all
ani=animation.FuncAnimation(fig, upd_p, init_func=init_p, interval=200,frames=(N_vid_fr if n_lines-N_vid_fr>0 else n_lines), repeat=False, blit=0)
plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
writervideo = animation.FFMpegWriter(fps=4)
#ani.save("traj_mp4/a="+str(a)+"_m="+str(m)+"_w="+str(omega)+"_vf="+str(v_fermi)+"_p_traj_evolution.mp4", writer=writervideo)
ani.save("beta="+str(round(beta,2))+"a="
        +str(round(beta/N_spots,2))+"pb="+str(round(p_bottom,2))+"w="+str(round(omega,2))+"p_traj_evolution.mp4", writer=writervideo)

f.close()
end_time=time.time()
print("elapsed time plotting p (seconds):",round(end_time-start_time,1))