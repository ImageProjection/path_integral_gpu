compile_L1_on:
	nvcc -o main -Xptxas -dlcm=ca main.cu

compile:
	nvcc -o main main.cu

run:
	./main

clean:
	rm main && find . -name 'out*' -delete
#to find regular expressions, rm can't

git_log:
	git log --all --graph --decorate
