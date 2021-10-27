/*Program models the behaviour of a particle in smooth phi4-like 
coordinate potential.
It produces trajectories
and a |\psi(x)|^2 graph. Computationally intensive code
runs on GPU (programmed with CUDA).*/

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <climits>
#include <curand.h>
#include <curand_kernel.h>
using namespace std;

#define print_traj_flag 0
#define N_spots 1024
#define N_bins 1024 //number of bins on x axis for histogram //not used yet
//#define hist_batch 512//how many points are classified simultaniously

void print_traj(FILE* out_traj,double* traj,double h_sigma)
{
	for (int i = 0; i < N_spots; i++)
	{
		fprintf(out_traj,"%.3lf ",traj[i]);
	}
	fprintf(out_traj,"%.6lf\n",h_sigma);
}

void print_hist(FILE* out_dens_plot, double* h_dens_plot, double range_start, double range_end)
{
	double bin_width=(double)(range_end-range_start)/N_bins;
	for (int i = 0; i < N_bins; i++)
	{
		fprintf(out_dens_plot,"%.8lf,%.8lf\n",range_start+i*bin_width,h_dens_plot[i]);
	}
}

void normalize_hist(unsigned int* h_hist, double* h_dens_plot, double range_start, double range_end)//also returns number of trajectory points used
{
	int val_sum=0;
	double integral;//==val_sum*delta_x
	double bin_width=(double)(range_end-range_start)/N_bins;//delta_x
	for(int i=0;i<N_bins;i++)
	{
		val_sum+=h_hist[i];
	}
	integral=val_sum*bin_width;
	for(int i=0;i<N_bins;i++)
	{
		h_dens_plot[i]=(double)h_hist[i]/integral;
	}
}

void h_histogram(double* h_traj, unsigned int* h_hist, double range_start, double range_end)
{
	int bin_i;
	double bin_width=(double)(range_end-range_start)/N_bins;
	for (int i = 0; i < N_spots; i++)
	{
		bin_i=int( (h_traj[i]-range_start)/bin_width );
		h_hist[bin_i]++;
	}	
}
/*
//hopefully this is not performance sensitive, because this is single threaded
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

	////new revision
	int id=threadIdx.x;
	int bin_i;
	__shared__ unsigned int hist[N_bins];
	__shared__ int desired_inc[N_spots];//i-th element contains which bin i-th thread wants to increment
	double bin_width=(double)(range_end-range_start)/N_bins;
	//init shmem
	for(int i=0; i<N_bins/N_spots; i++)
	{
		hist[id+i*N_spots]=0;
	}
	__syncthreads();
	//fill counters
	bin_i=int( (d_traj[id]-range_start)/bin_width );
	desired_inc[id]=bin_i;
	//atomicInc(&hist[bin_i],INT_MAX-1);maybe this is causing 700?
	__syncthreads();
	if (id==0)
	{
		for (int i = 0; i < N_spots; i++)
		{
			hist[desired_inc[i]]+=1;
		}		
	}
	__syncthreads();
	//move to dram
	for(int i=0; i<N_bins/N_spots; i++)
	{
		d_hist[id+i*N_spots]+=hist[id+i*N_spots];
	}
	__syncthreads();

	//TODO remove unnesessary synchs + maybe remove shared hist array
}
*/
__global__ void init_kernel(double* d_p_traj,double omega, double p_initial, int sigma_sweeps_period,
	double acc_rate_up_border, double acc_rate_low_border, double* d_sigma, int* d_accepted, curandState *d_rng_states)
{
	int id=threadIdx.x;// all threads must be in 1 block
	if (id==0)
    {
        *d_sigma=p_initial/3;
		*d_accepted=(sigma_sweeps_period*N_spots)*0.5*(acc_rate_low_border+acc_rate_up_border);//to not update on first sweep, when no data to eval acc_rate yet		
    }
	curand_init(id, id, 0, &d_rng_states[id]);
	d_p_traj[id]=p_initial;
}

__global__ void perform_sweeps(double* d_p_traj, double a, double v_fermi, double m, double omega, double p_bottom,
	double sigma_coef, int sigma_sweeps_period,
	double acc_rate_up_border, double acc_rate_low_border, int N_sweeps,
	double* d_sigma, int* d_accepted, curandState *d_rng_states)
{
	int id=threadIdx.x;// all threads must be in 1 block
    __shared__ double traj[N_spots];
    __shared__ double traj_new[N_spots];
    __shared__ int accepted_tmp_st[N_spots];
    __shared__ double sigma;//why shared? ans: all threads must have access to smae instant
    __shared__ double acc_rate;
    __shared__ int accepted;
	double p_left_node,p_old,p_new,S_old,S_new,prob_acc,gamma;
    
    accepted_tmp_st[id]=0;
	//load variables kept between calls
    if (id==0)
    {
        sigma=*d_sigma;
		accepted=*d_accepted;
    }
   __syncthreads();
   traj[id]=d_p_traj[id];
   __syncthreads();//redundant remove later
    for (int sweeps_counter=0; sweeps_counter < N_sweeps; sweeps_counter++)
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
		p_left_node=traj[(id-1+N_spots)%N_spots];
        p_old=traj[id];	
        p_new=p_old+sigma*curand_normal_double(&d_rng_states[id]);
		S_old=-(p_old-p_left_node)*(p_old-p_left_node)/(2*a*a*m*omega*omega) +v_fermi*sqrt(  (p_old*p_old-p_bottom*p_bottom)*(p_old*p_old-p_bottom*p_bottom)/(4*p_bottom*p_bottom) + m*m*v_fermi*v_fermi  );
        S_new=-(p_new-p_left_node)*(p_new-p_left_node)/(2*a*a*m*omega*omega) +v_fermi*sqrt(  (p_new*p_new-p_bottom*p_bottom)*(p_new*p_new-p_bottom*p_bottom)/(4*p_bottom*p_bottom) + m*m*v_fermi*v_fermi  );
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
	d_p_traj[id]=traj[id];    
	if (id==0)
    {
        *d_sigma=sigma;
		*d_accepted=accepted;
    }
}

__global__ void cumulative_transform(double* d_p_traj, double* d_x_traj)
{
	d_x_traj[0]=0;
	for (int j = 0; j < N_spots; j++)
	{
		for (int i = 0; i < j; i++)
		{
			d_x_traj[j]+=d_p_traj[i];
		}		
	}	
}

int main()
{
    clock_t start,end;
	start=clock();
	
	//metropolis parameters
	const int N_sweeps_waiting=300000;//initial termolisation length (in sweeps)
	const int N_sample_trajectories=500;//this many traj-s are used to build histogram
	const int Traj_sample_period=200;//it takes this time to evolve into new trajectory //do not choose 1
	const double a=0.035*2;
	double beta=a*N_spots;

	//sigma generation parameters for metropolis
	const int sigma_sweeps_period=1;
	const double sigma_coef=1.2;
	const double acc_rate_up_border=0.3;
	const double acc_rate_low_border=0.2;

	//hamiltonian parameters
	const double v_fermi=5.0;
	const double m=1.0;
	const double omega=20.0;
	const double p_bottom=1.0;//corresponds to 'bottom' of potential
	const double p_initial=p_bottom;//starting momentum value

	//histogram parameters
	const double x_range_start=-4.0;
	const double x_range_end=4.0;
	const double p_range_start=-4.0;
	const double p_range_end=4.0;

	//display parameters to terminal
	printf("===Particle in (actual) Twin Peaks potential===\n");
	printf("beta=%.2lf with a=%.4lf and N_spots=%d\n",beta,a,N_spots);
	printf("v_fermi=%.2lf\n",v_fermi);
	printf("p_bottom=%.2lf\n",p_bottom);
	printf("mass m=%.2lf\n",m);
	printf("omega=%.2lf\n",omega);
	printf("number of sample trajectories=%d\n",N_sample_trajectories);
	printf("Traj_sample_period=%d\n",Traj_sample_period);
	printf("cuda traj build ETA estimate (seconds): %.1f\n",(N_sweeps_waiting+N_sample_trajectories*Traj_sample_period)*8.6/410000);
	printf("subsequent python plotting ETA: %.1f\n",N_sample_trajectories*14.6/200);
	printf("total estimated ETA: %.1f\n",(N_sweeps_waiting+N_sample_trajectories*Traj_sample_period)*8.6/410000+N_sample_trajectories*14.6/200);
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	printf("kernel timeout enabled: %d\n",prop.kernelExecTimeoutEnabled);

	//open files for output
	FILE *out_p_traj;
	out_p_traj=fopen("out_p_traj.txt","w");
	FILE *out_p_dens_plot;
	out_p_dens_plot=fopen("out_p_dens_plot.txt","w");
	FILE *out_x_traj;
	out_x_traj=fopen("out_x_traj.txt","w");
	FILE *out_x_dens_plot;
	out_x_dens_plot=fopen("out_x_dens_plot.txt","w");
	
	//allocate memory for p and x trajs on cpu and gpu
	double* h_p_traj;
	h_p_traj=(double*)malloc(N_spots*sizeof(double));
	double* h_x_traj;
	h_x_traj=(double*)malloc(N_spots*sizeof(double));
	double* d_p_traj;
	cudaMalloc((void**)&d_p_traj, N_spots*sizeof(double));
	double* d_x_traj;
	cudaMalloc((void**)&d_x_traj, N_spots*sizeof(double));

	//allocate memory for p and x histograms and density plots (normalised histograms) on cpu and gpu
	unsigned int* h_p_hist;
	h_p_hist=(unsigned int*)malloc(N_bins*sizeof(int));
	for(int i; i<N_bins; h_p_hist[i++]=0);
	unsigned int* h_x_hist;
	h_x_hist=(unsigned int*)malloc(N_bins*sizeof(int));
	for(int i; i<N_bins; h_x_hist[i++]=0);
	double* h_p_dens_plot;
	h_p_dens_plot=(double*)malloc(N_bins*sizeof(double));
	double* h_x_dens_plot;
	h_x_dens_plot=(double*)malloc(N_bins*sizeof(double));
	unsigned int* d_p_hist;
	cudaMalloc((void**)&d_p_hist, N_bins*sizeof(unsigned int));
	cudaMemset(d_p_hist,0,N_bins*sizeof(unsigned int));
	unsigned int* d_x_hist;
	cudaMalloc((void**)&d_x_hist, N_bins*sizeof(unsigned int));
	cudaMemset(d_x_hist,0,N_bins*sizeof(unsigned int));

	//variables preserved between perf_sweeps calls (only for p)
	double h_sigma;
	double* d_sigma;
	cudaMalloc((void**)&d_sigma, sizeof(double));
	int* d_accepted;
	cudaMalloc((void**)&d_accepted, sizeof(int));	
	curandState *d_rng_states;
    cudaMalloc((void**)&d_rng_states, N_spots*sizeof(curandState));
	
	//sweeps and histograms kernel launch config
	dim3 grid_sweeps(1,1,1);
	dim3 block_sweeps(N_spots,1,1);
	dim3 grid_hist(1,1,1);//same for x and p
	dim3 block_hist(N_spots,1,1);
	dim3 grid_trans(1,1,1);
	dim3 block_trans(1,1,1);//only 1 thread, it's a small kernel

	//initialise p-trajectory, sigma, accepted and rng
	init_kernel<<<grid_sweeps,block_sweeps>>>(d_p_traj, omega, p_initial, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, d_sigma, d_accepted, d_rng_states);

	//run termolisation sweeps
	perform_sweeps<<<grid_sweeps,block_sweeps>>>(d_p_traj, a, v_fermi, m, omega, p_bottom, sigma_coef, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, N_sweeps_waiting, d_sigma, d_accepted, d_rng_states);
	//perform sweeps to build histogram and optionaly output trajectories
	for (int i=0; i<N_sample_trajectories; i++)
	{
		//plan
		//evolve p-trajectory
		//evaluate x-trajectory from it on gpu with small 1 thread kernel
		//add both trajectories data on cpu, add it to cumulative histograms
		//if flag is set, print both trajectories to files

		//evolve p-trajectory
		perform_sweeps<<<grid_sweeps,block_sweeps>>>(d_p_traj, a, v_fermi, m, omega, p_bottom, sigma_coef, sigma_sweeps_period,
			acc_rate_up_border, acc_rate_low_border, Traj_sample_period, d_sigma, d_accepted, d_rng_states);

		//compute x-trajectory
		////cumulative_transform<<<grid_trans,block_trans>>>(d_p_traj,d_x_traj);

		//add to cumulative histograms
		cudaMemcpy(h_p_traj,d_p_traj,N_spots*sizeof(double),cudaMemcpyDeviceToHost);
		h_histogram(h_p_traj, h_p_hist, p_range_start,p_range_end);
		////histogram<<<grid_hist,block_hist>>>(d_x_traj, d_x_hist, x_range_start, x_range_end);

		//print trajectories with appended sigma
		
		if (print_traj_flag)
		{
			cudaMemcpy(&h_sigma,d_sigma,sizeof(double),cudaMemcpyDeviceToHost);
			//////cudaMemcpy(h_p_traj,d_p_traj,N_spots*sizeof(double),cudaMemcpyDeviceToHost);
			////cudaMemcpy(h_x_traj,d_x_traj,N_spots*sizeof(double),cudaMemcpyDeviceToHost);
			print_traj(out_p_traj,h_p_traj,h_sigma);
			////print_traj(out_x_traj,h_x_traj,h_sigma);
		}
	}
	
	//copy histogram, normalize, build
	//////cudaMemcpy(h_p_hist,d_p_hist,N_bins*sizeof(unsigned int),cudaMemcpyDeviceToHost);
	////cudaMemcpy(h_x_hist,d_x_hist,N_bins*sizeof(unsigned int),cudaMemcpyDeviceToHost);
	normalize_hist(h_p_hist, h_p_dens_plot, p_range_start, p_range_end);
	////normalize_hist(h_x_hist, h_x_dens_plot, x_range_start, x_range_end);
	print_hist(out_p_dens_plot,h_p_dens_plot,p_range_start,p_range_end);
	////print_hist(out_x_dens_plot,h_x_dens_plot,x_range_start,x_range_end);
	
	//free memory
	free(h_p_traj);
	free(h_x_traj);
	free(h_p_dens_plot);
	free(h_x_dens_plot);
	free(h_p_hist);
	free(h_x_hist);
	cudaFree(d_p_traj);
	cudaFree(d_x_traj);
	cudaFree(d_p_hist);
	cudaFree(d_x_hist);
	cudaFree(d_sigma);
	cudaFree(d_accepted);
	cudaFree(d_rng_states);	

	//close files
	fclose(out_p_traj);
	fclose(out_x_traj);
	fclose(out_p_dens_plot);
	fclose(out_x_dens_plot);

	//check for errors and print report
	printf("===launch status report===\n");
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

	//test printf
	//double p_old=p_bottom;
	//printf("S is of order=%.6lf\n",-(p_old-p_bottom*0.95)*(p_old-p_bottom*0.95)/(2*a*a*m*omega*omega) +v_fermi*sqrt(  (p_old*p_old-p_bottom*p_bottom)*(p_old*p_old-p_bottom*p_bottom)/(4*p_bottom*p_bottom) + m*m*v_fermi*v_fermi  ));
}
