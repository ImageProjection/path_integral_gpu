#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <climits>
#include <curand.h>
#include <curand_kernel.h>
//#include "cuda_functions.h"
using namespace std;

#define N_spots 1024
#define Traj_sample_period 100000 //it takes this time to evolve into new trajectory
#define N_bins 1024 //number of bins on x axis for histogram //not used yet
#define hist_batch 512//how many points are classified simultaniously

void print_traj(FILE* out_traj,double* traj)
{
    for (int i = 0; i < N_spots; i++)
    {
        fprintf(out_traj,"%lf\n",traj[i]);
    }    
}

void print_hist(FILE* out_hist, unsigned int* hist, double range_start, double range_end)
{
	double bin_width=(double)(range_end-range_start)/N_bins;
	for (int i = 0; i < N_bins; i++)
	{
		fprintf(out_hist,"%.2lf %d\n",range_start+i*bin_width,hist[i]);
	}
}

__global__ void histogram(double* d_traj, unsigned int* d_hist, double range_start, double range_end)//N_bins and N_spots also
{
	int id=threadIdx.x;
	int bin_i;
	__shared__ unsigned int hist[N_bins];//array of counters
	double bin_width=(double)(range_end-range_start)/N_bins;
	//init shmem
	for(int i=0; i<N_bins/hist_batch; i++)
	{
		hist[id+i*hist_batch]=0;
	}
	__syncthreads();
	//fill counters
	for(int i=0; i<N_spots/hist_batch; i++)
	{
		bin_i=int( (d_traj[id+i*hist_batch]-range_start)/bin_width );
		atomicInc(&hist[bin_i],INT_MAX-1);
	}
	__syncthreads();
	//move to dram
	for(int i=0; i<N_bins/hist_batch; i++)
	{
		d_hist[id+i*hist_batch]+=hist[id+i*hist_batch];
	}
}

__global__ void init_kernel(double* d_traj,)

__global__ void perform_sweeps(double* d_traj, double a, double omega, double e,
	double bot,double p0, double sigma_coef, int sigma_sweeps_period,
	double acc_rate_up_border, double acc_rate_low_border, int N_sweeps_waiting, curandState *rng_states)
{
	int id=threadIdx.x;// all threads must be in 1 block
    __shared__ double traj[N_spots];
    __shared__ double traj_new[N_spots];
    __shared__ int accepted_tmp_st[N_spots];
    __shared__ double sigma;//will be dram, but cached immediately
    __shared__ double acc_rate;
    __shared__ int accepted;//try register type or rely on L1 cache
    double B;
	double p_old,p_new,S_old,S_new,prob_acc,gamma;
    
	curand_init(id, id, 0, &rng_states[id]);
	//init trajectory
	////instead load data from dram
    traj[id]=p0;
    accepted_tmp_st[id]=0;
    if (id==0)
    {
        sigma=0.01*sqrt(0.5/omega);
		accepted=(sigma_sweeps_period*N_spots)*0.5*(acc_rate_low_border+acc_rate_up_border);//to not update on first sweep, when no data to eval acc_rate yet		
    }
    __syncthreads();
    for (int sweeps_counter=0; sweeps_counter < N_sweeps_waiting; sweeps_counter++)
    {
        //update sigma
        if ( (sweeps_counter % sigma_sweeps_period) == 0)
		{
			for(int ps=N_spots/2; ps>=1; ps/=2)
			{
				if(id<ps)
					accepted_tmp_st[id]+=accepted_tmp_st[id+ps];
				__syncthreads();
			}
			accepted=accepted_tmp_st[0];
			if (id==0)
			{
				acc_rate=(double)accepted/(sigma_sweeps_period*N_spots);
				if (acc_rate < acc_rate_low_border)
				{
					sigma=sigma/sigma_coef;
				}
				if (acc_rate > acc_rate_up_border)
				{
					sigma=sigma*sigma_coef;
				}
			}
			accepted_tmp_st[id]=0;
		}
		__syncthreads();
        //local update for each
        p_old=traj[id];	
        p_new=p_old+sigma*curand_normal_double(&rng_states[id]);
        B=(traj[(id-1+N_spots)%N_spots]+traj[(id+1+N_spots)%N_spots]);
        S_old=(p_old*p_old-p_old*B)/a + a*omega*sqrt((p_old*p_old-1)*(p_old*p_old-1)+e*e);
        S_new=(p_new*p_new-p_new*B)/a + a*omega*sqrt((p_new*p_new-1)*(p_new*p_new-1)+e*e);
        if (S_new < S_old)
		{
			traj_new[id]=p_new;
			accepted_tmp_st[id]++;
		}
		else
		{
			prob_acc=1.0/exp(S_new-S_old);
			gamma=curand_uniform_double(&rng_states[id]);
			if (gamma < prob_acc)
			{
				traj_new[id]=p_new;
				accepted_tmp_st[id]++;
			}
		}
		//new -> old
		__syncthreads();
		traj[id]=traj_new[id];
	}
	//load to dram from shared
	d_traj[id]=traj[id];    
}

int main()
{
    clock_t start,end;
	start=clock();

	const int N_sweeps_waiting=800000;
	const double a=0.035;
	//const int N_spots=1024;
	//double beta=a*N_spots;
	const double omega=7.0;
	const double e=0.0;
	double bot=1.0;
	double p0=bot;
	const double range_start=-4.0;
	const double range_end=4.0;

	const int sigma_local_updates_period=2000;
	const int sigma_sweeps_period=ceil((double)sigma_local_updates_period/N_spots);
	const double sigma_coef=1.2;
	const double acc_rate_up_border=0.3;
	const double acc_rate_low_border=0.2;

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	printf("kernel timeout enabled: %d\n",
	prop.kernelExecTimeoutEnabled);

	FILE *out_traj;
	out_traj=fopen("out_traj.txt","w");
	FILE *out_hist;
	out_hist=fopen("out_hist.txt","w");

	double* h_traj;
	h_traj=(double*)malloc(N_spots*sizeof(double));
	double* d_traj;
	cudaMalloc((void**)&d_traj, N_spots*sizeof(double));

	unsigned int* h_hist;
	h_hist=(unsigned int*)malloc(N_bins*sizeof(int));
	unsigned int* d_hist;
	cudaMalloc((void**)&d_hist, N_bins*sizeof(unsigned int));
	cudaMemset(d_hist,0,N_bins*sizeof(unsigned int));

	dim3 grid(1,1,1);
	dim3 block(N_spots,1,1);
	
	curandState *devStates;
    cudaMalloc((void**)&devStates, N_spots*sizeof(curandState));

	//perform sweeps
	perform_sweeps<<<grid,block>>>(d_traj, a, omega, e, bot, p0, sigma_coef, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, N_sweeps_waiting, devStates);
	cudaMemcpy(h_traj,d_traj,N_spots*sizeof(double),cudaMemcpyDeviceToHost);	
	print_traj(out_traj,h_traj);

	//build histogram out of single trajectory //later it will run in cycle
	block.x=hist_batch;
	histogram<<<grid,block>>>(d_traj, d_hist, range_start,range_end);
	cudaMemcpy(h_hist,d_hist,N_bins*sizeof(unsigned int),cudaMemcpyDeviceToHost);
	print_hist(out_hist,h_hist,range_start,range_end);
		
	free(h_traj);
	cudaFree(d_traj);
	fclose(out_traj);
	fclose(out_hist);
		
	//check for errors
	cudaError_t err=cudaGetLastError();
	if (err != cudaSuccess)
	{
		printf("CUDA ERROR!!!\n");
		printf("err code: %d\n",err);
		if (err == 702)
		{
			printf("702 is similar to WDDM TDR false trigger; suggest running from tty3\n");
		}
	}
	else
	{
		printf("CUDA OK!!!\n");
	}

    end=clock();
	printf("TIME: %.2lf ms\n",(double)(end-start)/CLOCKS_PER_SEC*1000);
}
