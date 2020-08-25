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
#define N_bins 512 //number of bins on x axis for histogram //not used yet
#define hist_batch 512//how many points are classified simultaniously

void print_traj(FILE* out_traj,double* traj)
{
    for (int i = 0; i < N_spots; i++)
    {
        fprintf(out_traj,"%lf\n",traj[i]);
    }    
}

void print_hist(FILE* out_hist,int* hist, double range_start, double range_end)
{
	double bin_width=(double)(range_end-range_start)/N_bins;
	for (int i = 0; i < N_bins; i++)
	{
		fprintf(out_hist,"%.2lf %d\n",range_start+i*bin_width,hist[i]);
	}
}

__global__ void histogram(double* d_traj, int* d_hist, double range_start, double range_end)//N_bins and N_spots also
{
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	int bin_i;
	__shared__ unsigned int hist[N_bins];
	double bin_width=(double)(range_end-range_start)/N_bins;
	hist[id]=0;
	for(int i=0; i<N_spots/hist_batch; i++)
	{
		bin_i=int( (d_traj[id+i*hist_batch]-range_start)/bin_width );
		//bin_i=hist[int( (d_traj[id+i*hist_batch]-range_start)/bin_width )];
		atomicInc(&hist[bin_i],INT_MAX-1);
	}
	d_hist[id]+=hist[id];
}

__global__ void perform_sweeps(double* d_traj, double a, double omega,
	double bot,double x0, double sigma_coef, int sigma_sweeps_period,
	double acc_rate_up_border, double acc_rate_low_border, int N_sweeps_waiting, curandState *rng_states)
{
	int id=threadIdx.x;// all threads must be in 1 block
    __shared__ double traj[N_spots];
    __shared__ double traj_new[N_spots];
    __shared__ int accepted_tmp_st[N_spots];
    __shared__ double sigma;
    __shared__ double acc_rate;
    __shared__ int accepted;//try register type or rely on L1 cache
    double A=1.0-a*a*omega*omega*0.25;//try defferent memory types
    double B;
    double C=a*a*omega*omega/(bot*bot)*0.125;
	double x_old,x_new,S_old,S_new,prob_acc,gamma;
    
	curand_init(id, id, 0, &rng_states[id]);
    //init trajectory
    traj[id]=x0;
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
        x_old=traj[id];
        x_new=x_old+sigma*curand_normal_double(&rng_states[id]);
        B=(traj[(id-1+N_spots)%N_spots]+traj[(id+1+N_spots)%N_spots]);
        S_old=(A*x_old*x_old-B*x_old+C*x_old*x_old*x_old*x_old)/a;
        S_new=(A*x_new*x_new-B*x_new+C*x_new*x_new*x_new*x_new)/a;
        if (S_new < S_old)
		{
			traj_new[id]=x_new;
			accepted_tmp_st[id]++;
		}
		else
		{
			prob_acc=1.0/exp(S_new-S_old);
			gamma=curand_uniform_double(&rng_states[id]);
			if (gamma < prob_acc)
			{
				traj_new[id]=x_new;
				accepted_tmp_st[id]++;
			}
		}
		//new -> old
		__syncthreads();
		traj[id]=traj_new[id];
	}
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
	double bot=1.0;
	double x0=bot;
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

	int* h_hist;
	h_hist=(int*)malloc(N_bins*sizeof(int));
	int* d_hist;
	cudaMalloc((void**)&d_hist, N_bins*sizeof(int));
	cudaMemset(d_hist,0,N_bins*sizeof(int));

	dim3 grid(1,1,1);
	dim3 block(N_spots,1,1);
	
	curandState *devStates;
    cudaMalloc((void**)&devStates, N_spots*sizeof(curandState));

	//perform sweeps
	perform_sweeps<<<grid,block>>>(d_traj, a, omega, bot, x0, sigma_coef, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, N_sweeps_waiting, devStates);
	cudaMemcpy(h_traj,d_traj,N_spots*sizeof(double),cudaMemcpyDeviceToHost);	
	print_traj(out_traj,h_traj);

	//build histogram out of single trajectory //later it will run in cycle
	block.x=N_bins;
	histogram<<<grid,block>>>(d_traj, d_hist, range_start,range_end);
	cudaMemcpy(h_hist,d_hist,N_bins*sizeof(int),cudaMemcpyDeviceToHost);
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
	}
	else
	{
		printf("CUDA OK!!!\n");
	}

    end=clock();
	printf("TIME: %.2lf ms\n",(double)(end-start)/CLOCKS_PER_SEC*1000);
}
