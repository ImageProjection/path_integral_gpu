/*Program models the behaviour of a particle in Twin Peaks
potential using Monte-Carlo methods. It produces 1 trajectory p(t)
and a |\psi(p)|^2 graph. Computationally intensive code
runs on GPU (programmed with CUDA).*/

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <climits>
#include <curand.h>
#include <curand_kernel.h>
using namespace std;

#define N_spots 1024
#define N_bins 1024 //number of bins on x axis for histogram //not used yet
#define hist_batch 512//how many points are classified simultaniously

void print_traj(FILE* out_traj,double* traj)
{
    for (int i = 0; i < N_spots; i++)
    {
        fprintf(out_traj,"%d %lf\n", i, traj[i]);
    }    
}

void print_hist(FILE* out_dens_plot, double* h_dens_plot, double range_start, double range_end)
{
	double bin_width=(double)(range_end-range_start)/N_bins;
	for (int i = 0; i < N_bins; i++)
	{
		fprintf(out_dens_plot,"%.8lf %.8lf\n",range_start+i*bin_width,h_dens_plot[i]);
	}
}

int normalize_hist(unsigned int* h_hist, double* h_dens_plot, double range_start, double range_end)//also returns number of trajectory points used
{
	int n_points=0;
	double integral;//==n_points*delta_x
	double bin_width=(double)(range_end-range_start)/N_bins;//delta_x
	for(int i=0;i<N_bins;i++)
	{
		n_points+=h_hist[i];
	}
	integral=n_points*bin_width;
	for(int i=0;i<N_bins;i++)
	{
		h_dens_plot[i]=(double)h_hist[i]/integral;
	}
	return n_points;
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

__global__ void init_kernel(double* d_traj,double omega, double p0, int sigma_sweeps_period,
	double acc_rate_up_border, double acc_rate_low_border, double* d_sigma, int* d_accepted, curandState *d_rng_states)
{
	int id=threadIdx.x;// all threads must be in 1 block
	if (id==0)
    {
        *d_sigma=0.01*sqrt(0.5/omega);
		*d_accepted=(sigma_sweeps_period*N_spots)*0.5*(acc_rate_low_border+acc_rate_up_border);//to not update on first sweep, when no data to eval acc_rate yet		
    }
	curand_init(id, id, 0, &d_rng_states[id]);
	d_traj[id]=p0;
}

__global__ void perform_sweeps(double* d_traj, double a, double omega, double e,
	double sigma_coef, int sigma_sweeps_period,
	double acc_rate_up_border, double acc_rate_low_border, int N_sweeps_waiting,
	double* d_sigma, int* d_accepted, curandState *d_rng_states)
{
	int id=threadIdx.x;// all threads must be in 1 block
    __shared__ double traj[N_spots];
    __shared__ double traj_new[N_spots];
    __shared__ int accepted_tmp_st[N_spots];
    __shared__ double sigma;
    __shared__ double acc_rate;
    __shared__ int accepted;
    double B;
	double p_old,p_new,S_old,S_new,prob_acc,gamma;
    
    accepted_tmp_st[id]=0;
    if (id==0)
    {
        sigma=*d_sigma;
		accepted=*d_accepted;
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
        p_new=p_old+sigma*curand_normal_double(&d_rng_states[id]);
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
			gamma=curand_uniform_double(&d_rng_states[id]);
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
	if (id==0)
    {
        *d_sigma=sigma;
		*d_accepted=accepted;
    }
}

int main()
{
    clock_t start,end;
	start=clock();

	const int N_sweeps_waiting=800000;//initial termolisation length
	const int N_sample_trajectories=1;//this many traj-s are used to build histogram
	const int Traj_sample_period=100000;//it takes this time to evolve into new trajectory
	const double a=0.035;
	//const int N_spots=1024;//it's a define
	double beta=a*N_spots;
	const double omega=7.0;
	const double e=0.0;
	double bot=1.0;//corresponds to 'bottom' of potential
	double p0=bot;
	const double range_start=-4.0;//for histogram
	const double range_end=4.0;

	const int sigma_local_updates_period=2000;
	const int sigma_sweeps_period=ceil((double)sigma_local_updates_period/N_spots);
	const double sigma_coef=1.2;
	const double acc_rate_up_border=0.3;
	const double acc_rate_low_border=0.2;

	printf("===Particle in Twin Peaks potential===\n");
	printf("beta=%.2lf with a=%.4lf and N_spots=%d\n",beta,a,N_spots);
	printf("regularisation parameter e=%.2lf\n",e);
	printf("density plot resolution delta_x=%.5lf\n",(double)(range_end-range_start)/N_bins);
	printf("number of sample trajectories=%d\n",N_sample_trajectories);

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	printf("kernel timeout enabled: %d\n",prop.kernelExecTimeoutEnabled);

	//files
	FILE *out_traj;
	out_traj=fopen("out_traj.txt","w");
	FILE *out_dens_plot;
	out_dens_plot=fopen("out_dens_plot.txt","w");
	//trajectory
	double* h_traj;
	h_traj=(double*)malloc(N_spots*sizeof(double));
	double* d_traj;
	cudaMalloc((void**)&d_traj, N_spots*sizeof(double));
	//histogram
	unsigned int* h_hist;
	h_hist=(unsigned int*)malloc(N_bins*sizeof(int));
	double* h_dens_plot;
	h_dens_plot=(double*)malloc(N_bins*sizeof(double));
	unsigned int* d_hist;
	cudaMalloc((void**)&d_hist, N_bins*sizeof(unsigned int));
	cudaMemset(d_hist,0,N_bins*sizeof(unsigned int));
	//"globals" kept between evolve calls
	double* d_sigma;
	cudaMalloc((void**)&d_sigma, sizeof(double));
	int* d_accepted;
	cudaMalloc((void**)&d_accepted, sizeof(int));	
	curandState *d_rng_states;
    cudaMalloc((void**)&d_rng_states, N_spots*sizeof(curandState));
	
	//kernel launch config
	dim3 grid_sweeps(1,1,1);
	dim3 block_sweeps(N_spots,1,1);
	dim3 grid_hist(1,1,1);
	dim3 block_hist(hist_batch,1,1);
	

	//init kernel
	init_kernel<<<grid_sweeps,block_sweeps>>>(d_traj, omega, p0, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, d_sigma, d_accepted, d_rng_states);
	//termolise
	perform_sweeps<<<grid_sweeps,block_sweeps>>>(d_traj, a, omega, e, sigma_coef, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, N_sweeps_waiting, d_sigma, d_accepted, d_rng_states);
	//perform sweeps to build histogram
	for (int i=0; i<N_sample_trajectories; i++)
	{
		//evolve
		perform_sweeps<<<grid_sweeps,block_sweeps>>>(d_traj, a, omega, e, sigma_coef, sigma_sweeps_period,
			acc_rate_up_border, acc_rate_low_border, Traj_sample_period, d_sigma, d_accepted, d_rng_states);
		//add to cumulative histogram
		histogram<<<grid_hist,block_hist>>>(d_traj, d_hist, range_start,range_end);
	}


	//build last trajectory
	cudaMemcpy(h_traj,d_traj,N_spots*sizeof(double),cudaMemcpyDeviceToHost);	
	print_traj(out_traj,h_traj);
	//copy histogram, normalize, build
	cudaMemcpy(h_hist,d_hist,N_bins*sizeof(unsigned int),cudaMemcpyDeviceToHost);
	normalize_hist(h_hist, h_dens_plot, range_start, range_end);
	print_hist(out_dens_plot,h_dens_plot,range_start,range_end);
		
	free(h_traj);
	free(h_dens_plot);
	free(h_hist);
	cudaFree(d_traj);
	cudaFree(d_hist);
	cudaFree(d_sigma);
	cudaFree(d_accepted);
	cudaFree(d_rng_states);
	fclose(out_traj);
	fclose(out_dens_plot);

	printf("===launch status report===\n");
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
		if (err == 700)
		{
			printf("700 is out of range call\n");
		}
	}
	else
	{
		printf("No CUDA errors!!!\n");
	}

	end=clock();
	double total_time=(double)(end-start)/CLOCKS_PER_SEC;//in seconds
	printf("TOTAL TIME: %.1lf seconds (%.1lf minutes)\n",total_time,total_time/60);
}
