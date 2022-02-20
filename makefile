compile:
	nvcc -Wno-deprecated-gpu-targets -o main -arch=sm_35 main.cu

compile_L1_on:
	nvcc -o main -Xptxas -dlcm=ca main.cu

cr:
	nvcc -Wno-deprecated-gpu-targets -o main -arch=sm_35 main.cu && ./main

run:
	./main


#for NB
nb_full_run:
	g++ main.cpp && ./a.out && /usr/bin/python plotter.py && /usr/bin/python hist_plotter.py 
	
#for PC
pc_full_run:
	g++ main.cpp && ./a.out && python3 plotter.py && python3 hist_plotter.py

git_log:
	git log --all --graph --decorate

#force, remove folders with contents, remove ignored
clean:
	git clean -fdx
