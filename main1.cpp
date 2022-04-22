#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <cmath>
#include <random>
using namespace std;
random_device rd; 
mt19937_64 gen(rd()); 
normal_distribution<double> my_normal_double(0, 1); 

#define print_traj_flag 1//sample traj
#define print_termo_traj_flag 1//termo traj

const int N=490;
const double m=0.1;
const double omega=3.0;
const double p0=2.0;
const double a=0.01;
const double dtau=0.01;
const double beta= a*N_spots;
const int N_waiting_trajectories=200;
const int N_sample_trajectories=200;
const int T_md=20;

double p[N],oldp[N];
double P[N];
double action()
{
    double res=0.0;
    for(int i=0; i<N; i++)
    {
        int ir= (i+1)%N;
        res+=(p[ir]-p[i])*(p[ir]-p[i])/(2.0*a*m*omega*omega);
        res+=a*sqrt((p[i]*p[i] - p0*p0)*(p[i]*p[i] - p0*p0)/4.0/p0/p0 + mass*mass);
    }
    return res;
}

void update_moment(double dtau)
{
    for(int i=0; i<N; i++)
    {
        int ir=(i+1)%N;
        int il=(i-1+N)%N;
        double denomin = sqrt((p[i]*p[i] - p0*p0)*(p[i]*p[i] - p0*p0)/p0/p0 + 4.*mass*mass);
		P[i] -= dtau*((-p[ir]+2*p[i]-p[il])/a/mass/omega/omega + a*(p[i]*p[i] - p0*p0)*p[i]/p0/p0/denomin);

    }
}

void update_phi(double dtau)//md molecular update for traj
{
	for(int i = 0; i<N; i++)
		p[i] += dtau*P[i];
}

void run_md()
{
    update_moment(dtau/2);
	for(int t = 0; t<Ntau-1; t++) {
		update_phi(dtau);
		update_moment(dtau);
	}
	update_phi(dtau);
	update_moment(dtau/2);
}

double perform_sweeps()
{
    double accepted=0;
    for(int i=0;i<N;i++)
        oldp[i]=p[i];
    for(int i=0;i<N;i++)
        P[i]=my_normal_double(gen);
    double H_old=action()
    for(int i=0;i<N;i++)
        H_old+=P[i]*P[i]/2.0;

    run_md();

    double H_new=action();
    for(int i=0;i<N;i++)
        H_new+=P[i]*P[i]/2.0;
    if (H_new < H_old)
        accepted++;
        else
        {
            prob_acc=exp(H_old-H_new);
            gamma=(double)rand()/RAND_MAX;
            if (gamma < prob_acc)//then accept
                    accepted++;
                else//do not accept, thus revert
                    for(int i = 0; i<N; i++)
                        p[i] = oldp[i];
        }
}

void print_traj(FILE* out_traj,double* traj)
{
	for (int i = 0; i < N_spots; i++)
	{
		fprintf(out_traj,"%.8lf, ",traj[i]);
	}
	fprintf(out_traj,"%.6lf\n");
}



int main(int argc, char *argv[])
{
    struct timeval start, end;
	gettimeofday(&start, NULL);
	srand(start.tv_usec);

	//histogram parameters
	const double p_range=3;
	const double x_range=1000;//tweaked manually, values outside are discarded
	
	//traj range for plotter
	const double traj_p_range=3;
	const double traj_x_range=1000;
    /*
	//display parameters to terminal
	printf("===CPP CODE LAUNCH===\n");
	printf("beta=%.2lf with a=%.4lf and N_spots=%d\n",beta,a,N);
	printf("v_fermi=%.2lf\n",ham_params.v_fermi);
	printf("p_b=%.2lf\n",ham_params.p_b);
	printf("mass m=%.2lf\n",ham_params.m);
	printf("omega=%.2lf\n",ham_params.omega);
	printf("number of sample trajectories=%d\n",N_sample_trajectories);
	printf("N_waiting_trajectories=%d\n",N_waiting_trajectories);
	printf("N_sample_trajectories=%d\n",N_sample_trajectories);
	printf("N_steps_per_traj=%d\n",N_steps_per_traj);
	printf("N_cycles_per_step=%d\n",met_params.N_cycles_per_step);
	printf("T_molec=%d\n",met_params.T_molec);
	printf("T_lang=%d\n",met_params.T_lang);
    */
	//open files for output
	FILE *out_gen_des;//lists simulation parameters
	out_gen_des=fopen("out_gen_des.txt","w");
	FILE *out_p_traj;
	out_p_traj=fopen("out_p_traj.txt","w");
    /*
	//print general simulation description to file
	fprintf(out_gen_des,"N_spots,%d\n",N_spots);
	fprintf(out_gen_des,"N_waiting_trajectories,%d\n",N_waiting_trajectories);
	fprintf(out_gen_des,"N_sample_trajectories,%d\n",N_sample_trajectories);
	fprintf(out_gen_des,"N_steps_per_traj,%d\n",N_steps_per_traj);
	fprintf(out_gen_des,"a,%.8lf\n",a);
	fprintf(out_gen_des,"beta,%.8lf\n",beta);
	fprintf(out_gen_des,"v_fermi,%.8lf\n",ham_params.v_fermi);
	fprintf(out_gen_des,"m,%.8lf\n",ham_params.m);
	fprintf(out_gen_des,"omega,%.8lf\n",ham_params.omega);
	fprintf(out_gen_des,"p_b,%.8lf\n",ham_params.p_b);
	fprintf(out_gen_des,"p_range,%.8lf\n",p_range);
	fprintf(out_gen_des,"x_range,%.8lf\n",x_range);
	fprintf(out_gen_des,"traj_p_range,%.8lf\n",traj_p_range);
	fprintf(out_gen_des,"traj_x_range,%.8lf\n",traj_x_range);
	fprintf(out_gen_des,"sigma,%.8lf\n",sigma);
	fprintf(out_gen_des,"print_termo_traj_flag,%d\n",print_termo_traj_flag);



    */

    //init
    for(int i=0;i<N;i++)
        p[i]=oldp[i]=0.0;
    
    double accepted;
    //termo
    for(int i=0;i<N_waiting_trajectories;i++)
    {
        accepted=perform_sweeps();
		//print_traj(out_p_traj,h_p_traj);        
    }
    //sampling
    for(int i=0;i<N_waiting_trajectories;i++)
    {
        accepted=perform_sweeps();
		print_traj(out_p_traj,h_p_traj);        
    }



	fclose(out_p_traj);

	//check for errors and print report
	printf("===launch status report===\n");
	
	gettimeofday(&end, NULL);
	double total_time=((end.tv_sec  - start.tv_sec) * 1000000u + 
        end.tv_usec - start.tv_usec) / 1.e6;//in seconds
	printf("TOTAL TIME: %.1lf seconds (%.1lf minutes)\n",total_time,total_time/60);
	printf("===CPP CODE FINISHED WORKING===\n");
}
