#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define print_traj_flag 1
#define N_spots 1024
#define N_bins 1024
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
	int N_cycles_per_step;
	int T_molec;
	int T_lang;
	double e_molec;
	double e_lang;
};

double my_normal_double()//TODO try cuda for rng
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

void copy_traj(double* destination, double* source)
{
	for(int i=0; i<N_spots; i++)
	{
		destination[i]=source[i];
	}
}

double S(double* h_traj, struct hamiltonian_params_container ham_params)//action, PBC trajectory
{
	int S=0;
	double a=ham_params.a;
	double m=ham_params.m;
	double omega=ham_params.omega;
	double prev_node;

	for(int k=0; k<N_spots; k++)
	{
		prev_node=h_traj[(k-1+N_spots)%N_spots];
		S+= (h_traj[k]-prev_node)*(h_traj[k]-prev_node) / (2*a*a*m*omega*omega) + h_traj[k]*h_traj[k]/(2*m);
	}
	S *= a; 
	return S;
}

int perform_sweeps(double* h_p_traj, double* h_p_traj_new, 
	double* h_pi_vect, double* h_pi_vect_new, int N_steps,
	struct hamiltonian_params_container ham_params,
	struct metrop_params_container met_params)//h_p_traj_new (and both pi vectors) is purely for internal usage, but is allocated outside since it's 1 time	
{										 
	int accepted=0;	
	double a=ham_params.a;
	double m=ham_params.m;
	double omega=ham_params.omega;
	double p_prev_node,p_next_node,p_old,p_new,S_old,S_new,prob_acc,gamma;
    //TODO init pi somehow
    for (int steps_counter=0; steps_counter < N_steps; steps_counter++)
    {
		//evaluate trajectory for metropolis proposition
		for (int cycles_counter=0; cycles_counter < met_params.N_cycles_per_step; cycles_counter++)
		{
			//perform iterations using molecular dynamics algo
			for (int iteration_counter=0; iteration_counter < met_params.T_molec; iteration_counter++)
			{
				//TODO make propoper initialisation for pi and decide what to do with pi when using langevin,
				//especially when doing more then 1 langevin step (is it just translation, so we can keep pi?)
				//phi(1)=phi(0)+eps*pi(1/2)
				for(int i=0; i<N_spots; i++)
				{
					h_p_traj_new[i]=h_p_traj[i] + met_params.e_molec*h_pi_vect[i];
				}
				//pi(3/2)=pi(1/2)-eps*{ds}/{dphi(1)}
				for(int i=0; i<N_spots; i++)
				{
					p_prev_node=h_p_traj_new[(i-1+N_spots)%N_spots];
					p_next_node=h_p_traj_new[(i+1+N_spots)%N_spots];
					h_pi_vect_new[i]=h_pi_vect[i] - met_params.e_molec*( a*h_p_traj_new[i]/m + (2*h_p_traj_new[i]-(p_prev_node+p_next_node))/(a*m*omega*omega)  );
				}
				copy_traj(h_p_traj, h_p_traj_new);
				copy_traj(h_pi_vect, h_pi_vect_new);
			}
			//perform iterations using langevin algo
			for (int iteration_counter=0; iteration_counter < met_params.T_lang; iteration_counter++)
			{
				//phi(1)=phi(0)-eps_lang*{ds}/{dphi(0)} + sqrt(2eps_lang)*etta
				for(int i=0; i<N_spots; i++)
				{
					p_prev_node=h_p_traj_new[(i-1+N_spots)%N_spots];
					p_next_node=h_p_traj_new[(i+1+N_spots)%N_spots];
					h_p_traj_new[i]=h_p_traj[i] - met_params.e_lang*( a*h_p_traj[i]/m + (2*h_p_traj[i]-(p_prev_node+p_next_node))/(a*m*omega*omega)  ) + sqrt(2*met_params.e_lang)*my_normal_double();
				}
				copy_traj(h_p_traj, h_p_traj_new);
			}
		}
		//accept or discard this trajectory using standard metropolis fork
		S_old=S(h_p_traj, ham_params);
		S_new=S(h_p_traj_new, ham_params);
		if (S_new < S_old)
		{
			copy_traj(h_p_traj,h_p_traj_new);
		}
		else
		{
			prob_acc=1.0/exp(S_new-S_old);
			gamma=(double)rand()/RAND_MAX;
			if (gamma < prob_acc)
				{
					copy_traj(h_p_traj,h_p_traj_new);
					accepted++;
				}
		}

	}
	return accepted;//how many trajs of N_steps_per_traj were accepted
}

int main()
{
	srand(time(NULL));
    clock_t start,end;
	start=clock();
	//termo parameters
	const int N_steps_waiting=2000; //number of Metropolis steps to termolise the system
	const int N_sample_trajectories=300;//this many traj-s are used to build histogram
	const int N_steps_per_traj=100;//this many metropolis propositions are made for each of this traj-s
	const double a=0.015;//0.035*2;
	double beta=a*N_spots;

	//hamiltonian parameters
	struct hamiltonian_params_container ham_params;
	ham_params.v_fermi=50;
	ham_params.m=1;
	ham_params.omega=1;
	ham_params.p_bottom=2;//corresponds to 'bottom' of potential
	ham_params.a=a;

	//generation parameters for metropolis
	struct metrop_params_container met_params;
	met_params.p_initial=ham_params.p_bottom/3;
	met_params.N_cycles_per_step=10;
	met_params.T_molec=9;
	met_params.T_lang=1;//do not touch, unless it is pure Langevin
	met_params.e_lang=0.0005;
	met_params.e_molec=met_params.e_lang;//for correspondence

	//histogram parameters
	const double p_range=10;
	const double x_range=30;//tweaked manually, values outside are discarded
	
	//traj range for plotter
	const double traj_p_range=10;
	const double traj_x_range=30;

	//display parameters to terminal
	printf("===CPP CODE LAUNCH===\n");
	printf("beta=%.2lf with a=%.4lf and N_spots=%d\n",beta,a,N_spots);
	printf("v_fermi=%.2lf\n",ham_params.v_fermi);
	printf("p_bottom=%.2lf\n",ham_params.p_bottom);
	printf("mass m=%.2lf\n",ham_params.m);
	printf("omega=%.2lf\n",ham_params.omega);
	printf("number of sample trajectories=%d\n",N_sample_trajectories);
	printf("N_steps_waiting=%d\n",N_steps_waiting);
	printf("N_sample_trajectories=%d\n",N_sample_trajectories);
	printf("N_steps_per_traj=%d\n",N_steps_per_traj);
	printf("N_cycles_per_step=%d\n",met_params.N_cycles_per_step);
	printf("T_molec=%d\n",met_params.T_molec);
	printf("T_lang=%d\n",met_params.T_lang);
	//printf("cpp code ETA (seconds): %.1f\n",0.15*1e-6*N_sweeps_waiting*N_spots+3.78/3.2*1e-7*Traj_sample_period*N_sample_trajectories*N_spots);

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
	fprintf(out_gen_des,"N_steps_waiting,%d\n",N_steps_waiting);
	fprintf(out_gen_des,"N_sample_trajectories,%d\n",N_sample_trajectories);
	fprintf(out_gen_des,"N_steps_per_traj,%d\n",N_steps_per_traj);
	fprintf(out_gen_des,"a,%.4lf\n",a);
	fprintf(out_gen_des,"beta,%.4lf\n",beta);
	fprintf(out_gen_des,"v_fermi,%.4lf\n",ham_params.v_fermi);
	fprintf(out_gen_des,"m,%.4lf\n",ham_params.m);
	fprintf(out_gen_des,"omega,%.4lf\n",ham_params.omega);
	fprintf(out_gen_des,"p_bottom,%.4lf\n",ham_params.p_bottom);
	fprintf(out_gen_des,"p_range,%.4lf\n",p_range);
	fprintf(out_gen_des,"x_range,%.4lf\n",x_range);
	fprintf(out_gen_des,"p_range,%.4lf\n",traj_p_range);
	fprintf(out_gen_des,"p_range,%.4lf\n",traj_x_range);

	
	
	//allocate memory for p and x trajs on cpu (using heap for size)
	double* h_p_traj;
	h_p_traj=(double*)malloc(N_spots*sizeof(double));
    double* h_p_traj_new;
	h_p_traj_new=(double*)malloc(N_spots*sizeof(double));
	double* h_pi_vect;
	h_pi_vect=(double*)malloc(N_spots*sizeof(double));
	double* h_pi_vect_new;
	h_pi_vect_new=(double*)malloc(N_spots*sizeof(double));
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

	//init h_p_traj and h_pi_vect
	for(int i=0; i<N_spots; i++)
	{
		h_p_traj[i]=0;
	}
	double p_prev_node,p_next_node;
	for(int i=0; i<N_spots; i++)
	{
		p_prev_node=h_p_traj_new[(i-1+N_spots)%N_spots];
		p_next_node=h_p_traj_new[(i+1+N_spots)%N_spots];
		h_pi_vect[i]=-met_params.e_molec*0.5*( a*h_p_traj[i]/ham_params.m + (2*h_p_traj[i]-(p_prev_node+p_next_node))/(ham_params.a*ham_params.m*ham_params.omega*ham_params.omega)  );
	}

	//run termolisation sweeps
	//TODO configure to increase langevin part?
	double accepted=perform_sweeps(h_p_traj, h_p_traj_new, h_pi_vect, h_pi_vect_new, N_steps_waiting, ham_params, met_params);
	printf("Acceptance rate after running termolisation steps: %.4lf",accepted/N_steps_waiting);
	
	//perform sweeps to build histogram and optionaly output trajectories
	for (int i=0; i<N_sample_trajectories; i++)
	{
		//evolve p-trajectory
        accepted=perform_sweeps(h_p_traj, h_p_traj_new, h_pi_vect, h_pi_vect_new, N_steps_per_traj, ham_params, met_params);
		if (i%10==9)
		{
			printf("Acceptance rate after reaching p-traj No (%d) %.4lf\n",i,accepted/N_steps_per_traj);
		}

		//evaluate x-trajectory
		h_cumulative_transform(h_p_traj,h_x_traj,ham_params.a,ham_params.m);

		//add both trajectories points to cumulative histograms
		h_histogram(h_p_traj, h_p_hist, -p_range, p_range);
		h_histogram(h_x_traj, h_x_hist, -x_range, x_range);

		//print trajectories with appended acc rate (evaluated over steps made for this traj)		
		if (print_traj_flag)
		{
			print_traj(out_p_traj,h_p_traj,accepted/N_steps_per_traj);
			print_traj(out_x_traj,h_x_traj,accepted/N_steps_per_traj);
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
	printf("===CPP CODE FINISHED WORKING===\n");
}
