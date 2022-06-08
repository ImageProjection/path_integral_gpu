#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <cmath>
#include <random>

using namespace std;

random_device rd; 
mt19937_64 gen(rd()); 
normal_distribution<double> my_normal_double(0, 1); 

#define rare_printing_flag 1//used for very large number of trajectories
#define print_traj_flag 1//sample traj
#define print_termo_traj_flag 1//termo traj

int vid_traj=700; //it will be num_cores times that

const int N=340;//must be -333 for governor script to work
const double m=0.1;//not used
const double omega=3.0;
const double v_fermi=0.01;//not used
const double m_v_fermi=0.01;
//actually only these matter
const double p0=8.7;
const double A=0.66;//0.66; //eV times nm //h times v_fermi
const double B=0.1;//B //m over h
const double d=11.9;//d //h*h*mw^2


const double a=0.01;//a=0.01;
const double dtau=0.04;//epsilon for md
const double beta= a*N;
const int N_waiting_trajectories=200;//100
const int N_sample_trajectories=220;//40000 for 15 min
const int N_steps_per_traj=100;//200

const int T_md=20;

double p[N],oldp[N];
double P[N];

double action()
{
    double res=0.0;
    for(int i=0; i<N; i++)
    {
        int ir= (i+1)%N;
        res+=(p[ir]-p[i])*(p[ir]-p[i])/(2.0*a*d);
        res+=a*A*sqrt((p[i]*p[i] - p0*p0)*(p[i]*p[i] - p0*p0)/4.0/p0/p0 + B*B*A*A);
    }
    return res;
}
/* old vers
void update_moment(double dtau)
{
    for(int i=0; i<N; i++)
    {
        int ir=(i+1)%N;
        int il=(i-1+N)%N;
        double denomin = sqrt((p[i]*p[i] - p0*p0)*(p[i]*p[i] - p0*p0)/p0/p0 + 4.*m*m);
		P[i] -= dtau*((-p[ir]+2*p[i]-p[il])/a/m/omega/omega + a*(p[i]*p[i] - p0*p0)*p[i]/p0/p0/denomin);
    }
}
*/
void update_moment(double dtau)
{
    for(int i=0; i<N; i++)
    {
        int ir=(i+1)%N;
        int il=(i-1+N)%N;
        double denomin = sqrt((p[i]*p[i] - p0*p0)*(p[i]*p[i] - p0*p0)/p0/p0 + 4.0*B*B*A*A);
		P[i] -= dtau*((-p[ir]+2*p[i]-p[il])/a/d + a*A*(p[i]*p[i] - p0*p0)*p[i]/p0/p0/denomin);
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
	for(int t = 0; t<T_md-1; t++) {
		update_phi(dtau);
		update_moment(dtau);
	}
	update_phi(dtau);
	update_moment(dtau/2);
}

double perform_sweeps()
{
    double accepted=0;//it is 0 or 1
    for(int i=0;i<N;i++)
        oldp[i]=p[i];
    
    for(int i=0;i<N;i++)
        P[i]=my_normal_double(gen);
    double H_old=action();
    for(int i=0;i<N;i++)
        H_old+=P[i]*P[i]/2.0;

    run_md();

    double H_new=action();
    for(int i=0;i<N;i++)
        H_new+=P[i]*P[i]/2.0;

    double prob_acc,gamma;
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
    return accepted;
}

void print_traj(FILE* out_traj)
{
	for (int i=0; i < N; i++)
		fprintf(out_traj,"%.10lf ",p[i]);
	fprintf(out_traj,"\n");    
}

double aver_T()
{
    double res=0;
    for(int i=0;i<N;i++)
        res+=A*sqrt( B*B*A*A + (p[i]*p[i]-p0*p0)*(p[i]*p[i]-p0*p0)/(4*p0*p0) );
    return res/N;
}

double aver_p_dot()
{
    double res=0;
    int ir;
    for(int i=0;i<N;i++)
    {
        ir=(i+1)%N;
        res+= (p[i]-p[ir])*(p[i]-p[ir]);
    }
    res=res/(2*a*a*d);
    return res/N;
}

void print_log(FILE* out_log)
{
    double aver_E_val,aver_T_val,aver_p_dot_val;
    aver_T_val=aver_T();
    aver_p_dot_val=aver_p_dot();
    aver_E_val=0.5/a+aver_T_val-aver_p_dot_val;
	fprintf(out_log,"%.15lf, %.15lf, %.15lf\n",aver_E_val,aver_T_val,aver_p_dot_val);
}



int main(int argc, char *argv[])
{
    struct timeval start, end;
	gettimeofday(&start, NULL);
	srand(start.tv_usec);

	//histogram parameters
	const double p_range=3*p0;
	
	//traj range for plotter
	const double traj_p_range=3*p0;

	//open files for output
	FILE *out_gen_des;//lists simulation parameters
	out_gen_des=fopen("out_gen_des.txt","w");
	FILE *out_p_traj;
	out_p_traj=fopen("out_p_traj.txt","w");
    FILE *out_log;
	out_log=fopen("out_log.txt","w");

	//print general simulation description to file
	fprintf(out_gen_des,"N_spots,%d\n",N);
	fprintf(out_gen_des,"N_waiting_trajectories,%d\n",N_waiting_trajectories);
	fprintf(out_gen_des,"N_sample_trajectories,%d\n",N_sample_trajectories);
	fprintf(out_gen_des,"N_steps_per_traj,%d\n",N_steps_per_traj);
	fprintf(out_gen_des,"a,%.8lf\n",a);
	fprintf(out_gen_des,"beta,%.8lf\n",beta);
	fprintf(out_gen_des,"v_fermi,%.8lf\n",v_fermi);
	fprintf(out_gen_des,"m,%.8lf\n",m);
	fprintf(out_gen_des,"omega,%.8lf\n",omega);
	fprintf(out_gen_des,"p_b,%.8lf\n",p0);
	fprintf(out_gen_des,"p_range,%.8lf\n",p_range);
	fprintf(out_gen_des,"traj_p_range,%.8lf\n",traj_p_range);
	fprintf(out_gen_des,"A,%.8lf\n",A);
	fprintf(out_gen_des,"B,%.8lf\n",B);
	fprintf(out_gen_des,"d,%.8lf\n",d);


    //init
    for(int i=0;i<N;i++)
        p[i]=oldp[i]=0.0;
    
    double accepted=0;
    //termo
    for(int i=0;i<N_waiting_trajectories;i++)
    {
        for(int j=0; j<N_steps_per_traj; j++)
            accepted+=perform_sweeps();
        if (i%20==0)
            printf("acc_rate=%.2lf\n",accepted/N_steps_per_traj*100);
        accepted=0;
    }
    //sampling
    for(int i=0;i<N_sample_trajectories;i++)
    {
        for(int j=0; j<N_steps_per_traj; j++)
            accepted+=perform_sweeps();
        if (i%3000==0)
            printf("i=%d\n",i);     
        accepted=0;
        if (vid_traj>0)
        {
		    print_traj(out_p_traj);
            vid_traj--;
        }
        print_log(out_log);
    }



	fclose(out_p_traj);
    fclose(out_log);
    fclose(out_gen_des);
	//check for errors and print report
	printf("===launch status report===\n");
	
	gettimeofday(&end, NULL);
	double total_time=((end.tv_sec  - start.tv_sec) * 1000000u + 
        end.tv_usec - start.tv_usec) / 1.e6;//in seconds
	printf("TOTAL TIME: %.1lf seconds (%.1lf minutes)\n",total_time,total_time/60);

    FILE *out_dummy;//to signal governor script that this core is done
	out_dummy=fopen("out_dummy.txt","w");
    fprintf(out_dummy,"p_range,%.8lf\n",p_range);
    fclose(out_dummy);
	printf("===CPP CODE FINISHED WORKING===\n");
    return 0;
}
