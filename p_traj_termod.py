'''
input: trajectories file and gen_params_file
output: -into separate file, each line for 1 sample traj
            -aver_T
            -aver_V (over 1 sample traj)
            -aver_p_dot term
            -aver_E
            -rel
            -number of sign changes ("kink metric")
        -evaluates and puts results into yet another file:
            -global average T
            -global average V
            -global average p_dot term
            -global average E=V-T
            -global average rel
            -global average kink metric
            with statystical errors for each of these quantities
'''