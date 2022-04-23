compile_and_run_cpp:
	g++ main.cpp && ./a.out

compile:
	nvcc -Wno-deprecated-gpu-targets -o main -arch=sm_35 main.cu

compile_L1_on:
	nvcc -o main -Xptxas -dlcm=ca main.cu

cr:
	nvcc -Wno-deprecated-gpu-targets -o main -arch=sm_35 main.cu && ./main

run:
	./main

#for NB
nb_multi_beta_run:
	/usr/bin/python governor.py

#used for single beta, pase like: make recipe beta_val=11.06
nb_long_run:
	g++ -O2 main1.cpp && ./a.out $(beta_val) $(n_periods)\
	&& /usr/bin/python p_x_mp4_generator.py\
	&& /usr/bin/python p_x_plot_densities.py\
	&& /usr/bin/python p_traj_termod.py

nb_long_run_no_c:
	/usr/bin/python p_x_mp4_generator.py\
	&& /usr/bin/python p_x_plot_densities.py\
	&& /usr/bin/python p_traj_termod.py

nb_long_run_python_only:
	/usr/bin/python p_x_plot_densities.py\
	&& /usr/bin/python p_traj_termod.py

nb_full_run:
	g++ -O2 main.cpp && ./a.out && /usr/bin/python plotter.py && /usr/bin/python hist_plotter.py

nb_compile:
	g++ -O2 main.cpp

nb_run:
	./a.out && /usr/bin/python plotter.py && /usr/bin/python hist_plotter.py
	
#for PC
pc_full_run:
	g++ main.cpp && ./a.out && python3 plotter.py && python3 hist_plotter.py

git_log:
	git log --all --graph --decorate

#force, remove folders with contents, remove ignored
clean:
	git clean -fx
