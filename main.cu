/*Program models the behaviour of a particle with Twin Peaks hamiltonian
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

#define print_traj_flag 1
#define N_spots 1024
#define N_bins 1024
int discarded_x_points=0;//number of x-traj points which did not fit into histogram range

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
	double bin_width=double(range_end-range_start)/N_bins;
	double abs;
	for (int i = 0; i < N_spots; i++)
	{
		abs= ( (h_traj[i] >= 0) ? h_traj[i] : -h_traj[i] );
		if (abs < range_end)
		{
			bin_i=int( (h_traj[i]-range_start)/bin_width );
			h_hist[bin_i]+=1;
		}
		else
		{
			discarded_x_points++;
		}		
	}	
}

void h_cumulative_transform(double* h_p_traj, double* h_x_traj)
{
	//memset to 0
	for (int k = 0; k < N_spots; k++)
	{
		h_x_traj[k]=0;
	}
	//transform
	for (int j = 1; j < N_spots; j++)
	{
		for (int i = 1; i <= j; i++)
		{
			h_x_traj[j]+=h_p_traj[i];
		}		
	}	
}

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
		S_old=+(p_old-p_left_node)*(p_old-p_left_node)/(2*a*a*m*omega*omega) +v_fermi*sqrt(  (p_old*p_old-p_bottom*p_bottom)*(p_old*p_old-p_bottom*p_bottom)/(4*p_bottom*p_bottom) + m*m*v_fermi*v_fermi  );
        S_new=+(p_new-p_left_node)*(p_new-p_left_node)/(2*a*a*m*omega*omega) +v_fermi*sqrt(  (p_new*p_new-p_bottom*p_bottom)*(p_new*p_new-p_bottom*p_bottom)/(4*p_bottom*p_bottom) + m*m*v_fermi*v_fermi  );
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
	const int N_sweeps_waiting=100000;//initial termolisation length (in sweeps)
	const int N_sample_trajectories=10000;//this many traj-s are used to build histogram
	const int Traj_sample_period=500;//it takes this time to evolve into new trajectory //do not choose 1
	const double a=0.035*2;
	double beta=a*N_spots;

	//sigma generation parameters for metropolis
	const int sigma_sweeps_period=1;
	const double sigma_coef=1.2;
	const double acc_rate_up_border=0.3;
	const double acc_rate_low_border=0.2;

	//hamiltonian parameters
	const double v_fermi=500;
	const double m=0.05;
	const double omega=200;//200 is dense kinks
	const double p_bottom=2;//corresponds to 'bottom' of potential
	const double p_initial=p_bottom;//starting momentum value

	//histogram parameters, will be updated
	const double p_range=5;
	const double x_range=175;//tweaked manually, values outside are discarded

	//display parameters to terminal
	printf("===Particle with (actual) Twin Peaks hamiltonian===\n");
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
	double* d_x_traj;//TODO remove
	cudaMalloc((void**)&d_x_traj, N_spots*sizeof(double));

	//allocate memory for p and x histograms and density plots (normalised histograms) on cpu and gpu
	unsigned int* h_p_hist;
	h_p_hist=(unsigned int*)malloc(N_bins*sizeof(int));
	for(int i=0; i<N_bins; i++)
		h_p_hist[i]=0;
	unsigned int* h_x_hist;
	h_x_hist=(unsigned int*)malloc(N_bins*sizeof(int));
	for(int i=0; i<N_bins; i++)
		h_x_hist[i]=0;
	double* h_p_dens_plot;
	h_p_dens_plot=(double*)malloc(N_bins*sizeof(double));
	double* h_x_dens_plot;
	h_x_dens_plot=(double*)malloc(N_bins*sizeof(double));

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

	//initialise p-trajectory, sigma, accepted and rng
	init_kernel<<<grid_sweeps,block_sweeps>>>(d_p_traj, omega, p_initial, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, d_sigma, d_accepted, d_rng_states);

	//run termolisation sweeps
	perform_sweeps<<<grid_sweeps,block_sweeps>>>(d_p_traj, a, v_fermi, m, omega, p_bottom, sigma_coef, sigma_sweeps_period,
		acc_rate_up_border, acc_rate_low_border, N_sweeps_waiting, d_sigma, d_accepted, d_rng_states);
	
	//perform sweeps to build histogram and optionaly output trajectories
	for (int i=0; i<N_sample_trajectories; i++)
	{
		//evolve p-trajectory
		perform_sweeps<<<grid_sweeps,block_sweeps>>>(d_p_traj, a, v_fermi, m, omega, p_bottom, sigma_coef, sigma_sweeps_period,
			acc_rate_up_border, acc_rate_low_border, Traj_sample_period, d_sigma, d_accepted, d_rng_states);

		//copy it to host
		cudaMemcpy(h_p_traj,d_p_traj,N_spots*sizeof(double),cudaMemcpyDeviceToHost);

		//evaluate x-trajectory
		h_cumulative_transform(h_p_traj,h_x_traj);

		//add both trajectories points to cumulative histograms
		h_histogram(h_p_traj, h_p_hist, -p_range, p_range);
		h_histogram(h_x_traj, h_x_hist, -x_range, x_range);

		//print trajectories with appended sigma		
		if (print_traj_flag)
		{
			cudaMemcpy(&h_sigma,d_sigma,sizeof(double),cudaMemcpyDeviceToHost);
			print_traj(out_p_traj,h_p_traj,h_sigma);
			print_traj(out_x_traj,h_x_traj,h_sigma);
		}
	}
	
	//copy, normalize and plot histograms to file
	normalize_hist(h_p_hist, h_p_dens_plot, -p_range, p_range);
	normalize_hist(h_x_hist, h_x_dens_plot, -x_range, x_range);
	print_hist(out_p_dens_plot,h_p_dens_plot,-p_range,p_range);
	print_hist(out_x_dens_plot,h_x_dens_plot,-x_range,x_range);
	
	//free memory
	free(h_p_traj);
	free(h_x_traj);
	free(h_p_dens_plot);
	free(h_x_dens_plot);
	free(h_p_hist);
	free(h_x_hist);
	cudaFree(d_p_traj);
	cudaFree(d_x_traj);
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

	printf("total number of histogram-discarded x points: %d (%.2lf%)\n",discarded_x_points,(double)discarded_x_points/(N_sample_trajectories*N_spots));
	end=clock();
	double total_time=(double)(end-start)/CLOCKS_PER_SEC;//in seconds
	printf("TOTAL TIME: %.1lf seconds (%.1lf minutes)\n",total_time,total_time/60);
}
