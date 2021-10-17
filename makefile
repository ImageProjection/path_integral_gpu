compile:
	nvcc -Wno-deprecated-gpu-targets -o main -arch=sm_35 main.cu

compile_L1_on:
	nvcc -o main -Xptxas -dlcm=ca main.cu

cr:
	nvcc -Wno-deprecated-gpu-targets -o main -arch=sm_35 main.cu && ./main

run:
	./main

plot:
	python3 plotter.py

histogram:
	python3 hist_plotter.py

full_run:
	nvcc -Wno-deprecated-gpu-targets -o main -arch=sm_35 main.cu && ./main && python3 plotter.py && python3 hist_plotter.py

clean:
	rm main && find . -name 'out*' -delete
#to find regular expressions, rm can't

git_log:
	git log --all --graph --decorate
