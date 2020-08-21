#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <curand.h>
#include <curand_kernel.h>
//#include "cuda_functions.h"
using namespace std;

#define N_spots 1024

void print_traj(FILE* out_traj,double* traj)
{
    for (int i = 0; i < N_spots; i++)
    {
        fprintf(out_traj,"%lf\n",traj[i]);
    }    
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
        S_old=1;//(A*x_old*x_old-B*x_old+C*x_old*x_old*x_old*x_old)/a;
        S_new=1;//(A*x_new*x_new-B*x_new+C*x_new*x_new*x_new*x_new)/a;
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

	double* h_traj;
	h_traj=(double*)malloc(N_spots*sizeof(double));
	double* d_traj;
	cudaMalloc((void**)&d_traj, N_spots*sizeof(double));

	dim3 grid(1,1,1);
	dim3 block(N_spots,1,1);
	
	curandState *devStates;
    cudaMalloc((void**)&devStates, N_spots*sizeof(curandState));


	perform_sweeps<<<grid,block>>>(d_traj, a, omega, bot, x0, sigma_coef, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, N_sweeps_waiting, devStates);
		cudaMemcpy(h_traj,d_traj,N_spots*sizeof(double),cudaMemcpyDeviceToHost);
		
	print_traj(out_traj,h_traj);
		
	free(h_traj);
	cudaFree(d_traj);
	fclose(out_traj);
		
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
