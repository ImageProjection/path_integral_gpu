/*Program models the behaviour of a particle with Twin Peaks hamiltonian
It produces trajectories
and a |\psi(x)|^2 graph.*/

//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

//using namespace std;

#define print_traj_flag 1
#define N_spots 3000
#define N_bins 2000
int discarded_x_points=0;//number of x-traj points which did not fit into histogram range

struct hamiltonian_params_container
{
	double v_fermi;
	double m;
	double omega;
	double p_bottom;
	double a;
};
struct metrop_params_container
{
	int sigma_sweeps_period;
	double sigma_coef;
	double acc_rate_up_border;
	double acc_rate_low_border;
	double p_initial;
};

double my_normal_double()
{
	double g1,g2;
	g1=(double)rand()/RAND_MAX;
	g2=(double)rand()/RAND_MAX;
	return sqrt(-2*log(g2))*sin(2*M_PI*g1);
}

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
			bin_i=floor( (h_traj[i]-range_start)/bin_width );
			h_hist[bin_i]+=1;
		}
		else
		{
			discarded_x_points++;
		}		
	}	
}

void h_cumulative_transform(double* h_p_traj, double* h_x_traj,double a,double m)
{
	//memset to 0
	for (int k = 0; k < N_spots; k++)
	{
		h_x_traj[k]=0;
	}
	//transform
	h_x_traj[0]=h_p_traj[0]*a/m;
	for (int j = 1; j < N_spots; j++)
	{
		h_x_traj[j]=h_x_traj[j-1]+h_p_traj[j]*a/m;
	}
	/*
	for (int j = 1; j < N_spots; j++)
	{
		for (int i = 1; i <= j; i++)
		{
			h_x_traj[j]+=h_p_traj[i];
		}		
	}
	*/	
}

double perform_sweeps(double* h_p_traj, double* h_p_traj_new, int N_sweeps,
	struct hamiltonian_params_container ham_params,
	struct metrop_params_container met_params)//h_p_traj_new is purely for internal usage,
{										 //but is allocated outside since it's 1 time	
    static double acc_rate;
    static double sigma=met_params.p_initial/3;
    static int accepted=(met_params.sigma_sweeps_period*N_spots)*
		0.5*(met_params.acc_rate_low_border+met_params.acc_rate_up_border);
	
	double p_left_node,p_right_node,p_old,p_new,S_old,S_new,prob_acc,gamma;
    
    for (int sweeps_counter=0; sweeps_counter < N_sweeps; sweeps_counter++)
    {
        //update sigma
        if ( (sweeps_counter % met_params.sigma_sweeps_period) == 0)
		{   
			acc_rate=(double)accepted/(met_params.sigma_sweeps_period*N_spots);
			accepted=0;
			if (acc_rate < met_params.acc_rate_low_border)
			{
				sigma/=met_params.sigma_coef;
			}
			if (acc_rate > met_params.acc_rate_up_border)
			{
				sigma*=met_params.sigma_coef;
			}
		}
        //local update for each, write into new traj
		for (int i = 0; i < N_spots; i++)
		{
			p_left_node=h_p_traj[(i-1+N_spots)%N_spots];
			p_right_node=h_p_traj[(i+1+N_spots)%N_spots];
        	p_old=h_p_traj[i];	
        	p_new=p_old+sigma*my_normal_double();
			S_old=(p_old*p_old-p_old*(p_left_node+p_right_node))/(ham_params.a*ham_params.m*ham_params.omega*ham_params.omega) + p_old*p_old/2/ham_params.m;//sqrt(p_old*p_old+m*m);
			S_new=(p_new*p_new-p_new*(p_left_node+p_right_node))/(ham_params.a*ham_params.m*ham_params.omega*ham_params.omega) + p_new*p_new/2/ham_params.m;//sqrt(p_new*p_new+m*m);

			if (S_new < S_old)
			{
				h_p_traj_new[i]=p_new;
				accepted++;
			}
			else
			{
				prob_acc=1.0/exp(S_new-S_old);
				gamma=(double)rand()/RAND_MAX;
				if (gamma < prob_acc)
				{
					h_p_traj_new[i]=p_new;
					accepted++;
				}
			}
		}
		//new -> old
        for (int k = 0; k < N_spots; k++)
        {
            h_p_traj[k]=h_p_traj_new[k];
        }
	}
	return sigma;
}

double average_square(double* h_traj)
{
	double sum=0;
	for (int i = 0; i < N_spots; i++)
	{
		sum+=h_traj[i]*h_traj[i];
	}
	return sum/N_spots;
}

int main()
{
	srand(time(NULL));
    clock_t start,end;
	start=clock();
	//termo parameters
	const int N_sweeps_waiting=10000;//initial termolisation length (in sweeps)
	const int N_sample_trajectories=100;//this many traj-s are used to build histogram
	const int Traj_sample_period=200;//it takes this time to evolve into new trajectory //do not choose 1
	const double a=1;//0.035*2;
	double beta=a*N_spots;

	//hamiltonian parameters
	struct hamiltonian_params_container ham_params;
	ham_params.v_fermi=50;
	ham_params.m=1;
	ham_params.omega=2000;//200 is dense kinks
	ham_params.p_bottom=2;//corresponds to 'bottom' of potential
	ham_params.a=a;

	//sigma generation parameters for metropolis
	struct metrop_params_container met_params;
	met_params.sigma_sweeps_period=10;//try bigger
	met_params.sigma_coef=1.2;
	met_params.acc_rate_up_border=0.3;
	met_params.acc_rate_low_border=0.2;
	met_params.p_initial=ham_params.p_bottom/3;

	//histogram parameters, will be updated
	const double p_range=10;
	const double x_range=300;//tweaked manually, values outside are discarded

	//display parameters to terminal
	printf("===Particle with (actual) Twin Peaks hamiltonian===\n");
	printf("beta=%.2lf with a=%.4lf and N_spots=%d\n",beta,a,N_spots);
	printf("v_fermi=%.2lf\n",ham_params.v_fermi);
	printf("p_bottom=%.2lf\n",ham_params.p_bottom);
	printf("mass m=%.2lf\n",ham_params.m);
	printf("omega=%.2lf\n",ham_params.omega);
	printf("number of sample trajectories=%d\n",N_sample_trajectories);
	printf("Traj_sample_period=%d\n",Traj_sample_period);
	printf("python plotting ETA: %.1f\n",N_sample_trajectories*14.6/200);

	//open files for output
	FILE *out_gen_des;//lists simulation parameters
	out_gen_des=fopen("out_gen_des.txt","w");
	FILE *out_energies;
	out_energies=fopen("out_energies.txt","w");
	double aver_T,aver_V;
	FILE *out_p_traj;
	out_p_traj=fopen("out_p_traj.txt","w");
	FILE *out_p_dens_plot;
	out_p_dens_plot=fopen("out_p_dens_plot.txt","w");
	FILE *out_x_traj;
	out_x_traj=fopen("out_x_traj.txt","w");
	FILE *out_x_dens_plot;
	out_x_dens_plot=fopen("out_x_dens_plot.txt","w");

	//print general simulation description to file
	fprintf(out_gen_des,"N_spots,%d\n",N_spots);
	fprintf(out_gen_des,"N_sweeps_waiting,%d\n",N_sweeps_waiting);
	fprintf(out_gen_des,"N_sample_trajectories,%d\n",N_sample_trajectories);
	fprintf(out_gen_des,"Traj_sample_period,%d\n",Traj_sample_period);
	fprintf(out_gen_des,"a,%.4lf\n",a);
	fprintf(out_gen_des,"beta,%.4lf\n",beta);
	fprintf(out_gen_des,"v_fermi,%.4lf\n",ham_params.v_fermi);
	fprintf(out_gen_des,"m,%.4lf\n",ham_params.m);
	fprintf(out_gen_des,"omega,%.4lf\n",ham_params.omega);
	fprintf(out_gen_des,"p_bottom,%.4lf\n",ham_params.p_bottom);
	fprintf(out_gen_des,"p_range,%.4lf\n",p_range);
	fprintf(out_gen_des,"x_range,%.4lf\n",x_range);

	
	
	//allocate memory for p and x trajs on cpu (using heap for size)
	double* h_p_traj;
	h_p_traj=(double*)malloc(N_spots*sizeof(double));
    double* h_p_traj_new;
	h_p_traj_new=(double*)malloc(N_spots*sizeof(double));
	double* h_x_traj;
	h_x_traj=(double*)malloc(N_spots*sizeof(double));

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
/*
	//variables preserved between perf_sweeps calls (only for p)
	double h_sigma;
    double h_accepted;

      maybe its worth generating rng on gpu
	cudaMalloc((void**)&d_accepted, sizeof(int));	
	curandState *d_rng_states;
    cudaMalloc((void**)&d_rng_states, N_spots*sizeof(curandState));
    
	
	//sweeps and histograms kernel launch config

	//initialise p-trajectory, sigma, accepted and rng
    h_sigma=p_initial/3;
    h_accepted=(sigma_sweeps_period*N_spots)*0.5*(acc_rate_low_border+acc_rate_up_border);//to not update on first sweep, when no data to eval acc_rate yet		
    for (int i = 0; i < N_spots; i++)
    {
        h_p_traj[i]=p_initial;
    }
 */   

	double h_sigma;//extracts a value of sigma for later printing
	//run termolisation sweeps
	h_sigma=perform_sweeps(h_p_traj, h_p_traj_new, N_sweeps_waiting, ham_params, met_params);
	
	//perform sweeps to build histogram and optionaly output trajectories
	for (int i=0; i<N_sample_trajectories; i++)
	{
		//evolve p-trajectory
        h_sigma=perform_sweeps(h_p_traj, h_p_traj_new, Traj_sample_period, ham_params, met_params);

		//evaluate x-trajectory
		h_cumulative_transform(h_p_traj,h_x_traj,ham_params.a,ham_params.m);

		//add both trajectories points to cumulative histograms
		h_histogram(h_p_traj, h_p_hist, -p_range, p_range);
		h_histogram(h_x_traj, h_x_hist, -x_range, x_range);

		//print trajectories with appended sigma		
		if (print_traj_flag)
		{
			print_traj(out_p_traj,h_p_traj,h_sigma);
			print_traj(out_x_traj,h_x_traj,h_sigma);
		}

		//evaluate energies corresponding to each trajectory
		aver_T=average_square(h_p_traj)/(2*ham_params.m);
		aver_V=average_square(h_x_traj)*(ham_params.m*ham_params.omega*ham_params.omega/2);
		fprintf(out_energies,"%.3lf, %.3lf, %.3lf\n",aver_T,aver_V,ham_params.omega/4);
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
	
	//close files
	fclose(out_p_traj);
	fclose(out_x_traj);
	fclose(out_p_dens_plot);
	fclose(out_x_dens_plot);

	//check for errors and print report
	printf("===launch status report===\n");
	
	printf("total number of histogram-discarded x points: %d (%.2lf%)\n",discarded_x_points,(double)discarded_x_points/(N_sample_trajectories*N_spots));
	end=clock();
	double total_time=(double)(end-start)/CLOCKS_PER_SEC;//in seconds
	printf("TOTAL TIME: %.1lf seconds (%.1lf minutes)\n",total_time,total_time/60);

	printf("my double asdnormal: %.3lf\n",(double)rand()/RAND_MAX);
	printf("my double asdnormal: %.3lf\n",(double)rand()/RAND_MAX);
	printf("my double asdnormal: %.3lf\n",(double)rand()/RAND_MAX);
	printf("my double asdnormal: %.3lf\n",(double)rand()/RAND_MAX);
	printf("my double asdnormal: %.3lf\n",(double)rand()/RAND_MAX);
	printf("%d\n",RAND_MAX+1);

}
