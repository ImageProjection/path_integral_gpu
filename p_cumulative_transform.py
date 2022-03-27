'''
input: p-trajectories file
output: x-trajectories file of same format
'''

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

def cumulative_transform(p_traj):
    sum=0
    x_traj=[]
    x_traj.append(p_traj[0]*a/m)
    for i in range(N_spots):
        x_traj.append(x_traj[i-1]+p_traj[i]*a/m)
    x_traj.append(p_traj[N_spots])#carry over acc_rate
    return x_traj

        
#main
p_file=open("out_p_traj.txt",'r')
x_file=open("out_x_traj.txt",'w')

p_traj=[]
for line in p_file:
    #obtain list of points
    p_traj=list(map(float,line.split(",")))
    #transform
    x_traj=cumulative_transform(p_traj)
    #put into x_traj file    
    x_file.write(", ".join(str(round(t,8)) for t in x_traj) +"\n")

p_file.close()
x_file.close()