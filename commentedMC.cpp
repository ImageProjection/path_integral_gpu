#define _CRT_SECURE_NO_WARNINGS

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <map>

//=======================================================================================
//
//                            MERSENNE GENERATOR
//

/* Period parameters */
#define MERSENNE_N 624
#define MERSENNE_M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned int mt[MERSENNE_N]; /* the array for the state vector  */
int mersenne_i = -1; /*  < 0 means mt[N] is not initialized */
double mersenne_array[MERSENNE_N];

/* initializing the array with a NONZERO seed */
void seed_mersenne(long seed)
{
	int mti;
	mt[0] = seed & 0xffffffffUL;
	for(mti = 1; mti<MERSENNE_N; mti++) {
		mt[mti] =
			(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
	mersenne_i = 0;
}

double mersenne_generate()
{
	register unsigned int y;
	register int kk;
	static unsigned int mag01[2] = { 0x0, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if(mersenne_i < 0) {  /* if sgenrand() has not been called, */
		printf("DUMMY: you did not seed the generator!\n");
		exit(0);
	}

	/* generate MERSENNE_N words at one time */

	for(kk = 0; kk<MERSENNE_N-MERSENNE_M; kk++) {
		y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
		mt[kk] = mt[kk+MERSENNE_M] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	for(; kk<MERSENNE_N-1; kk++) {
		y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
		mt[kk] = mt[kk+(MERSENNE_M-MERSENNE_N)] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	y = (mt[MERSENNE_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
	mt[MERSENNE_N-1] = mt[MERSENNE_M-1] ^ (y >> 1) ^ mag01[y & 0x1];

	for(kk = 0; kk<MERSENNE_N; kk++) {
		y = mt[kk];
		y ^= TEMPERING_SHIFT_U(y);
		y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
		y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
		y ^= TEMPERING_SHIFT_L(y);
		mersenne_array[kk] = (double)y * 2.3283064365386963e-10;  /* reals: interval [0,1) */
	}

	mersenne_i = MERSENNE_N;
	return (mersenne_array[--mersenne_i]);
	/* return y; */ /* for integer generation */
}

#define mersenne() ( mersenne_i > 0 ? mersenne_array[--mersenne_i] : mersenne_generate() )

double BoxMuller2(double mean, double stdDev) {

	static bool hasSpare = false;
	static double spare;

	if(hasSpare) {
		hasSpare = false;
		return mean + stdDev * spare;
	}

	hasSpare = true;
	static double u, v, s;
	do {
		u = mersenne() * 2.0 - 1.0;
		v = mersenne() * 2.0 - 1.0;
		s = u * u + v * v;
	} while((s >= 1.0) || (s == 0.0));

	s = sqrt(-2.0 * log(s) / s);
	spare = v * s;
	return mean + stdDev * u * s;
}

/////////////////////////////////////////////////////////////////////////////////////////

//#define TEST

const double omega = 3.;
const double mass = 0.1;
const double p0 = 2.;
const int N = 490;	//500
const double beta = 0.01*N;
#ifdef TEST
const int LEN = 1000;
#else
const int LEN = 1000000*20;//number of sample trajectories
#endif
const double dtau = 1.e-2;//molecular or langevin epsilon
const int Ntau = 20;//T_molec, T_lang==1

double p[N], oldp[N];//p_traj_new, p_traj_old
int acc;
double dH;

const double hist_d = 0.5;		//for x
//const double hist_d = 0.02;	//for p
int hist_tot;
std::map<int, int> hist;

double x[N];
void makex()
{
	double delta = beta/N;
	x[0] = 0.;
	for(int i = 0; i<N-1; i++)
		x[i+1] = x[i] + delta*p[i]/mass;
}

void puthist(const double *arr)
{
	for(int i = 0; i<N; i++) {
		int j = fabs(arr[i])/hist_d;
		int idx = arr[i]>0 ? j : -(j+1);
		hist[idx]++;
	}
	hist_tot++;
}

void init()//set both p_traj to 0
{
	for(int i = 0; i<N; i++)
		oldp[i] = p[i] = 0.;
}

double action()//twin-peaks action
{
	double act = 0.;
	double delta = beta/N;//time step
	for(int i = 0; i<N; i++) {
		int ir = (i+1)%N;
		act += (p[ir]-p[i])*(p[ir]-p[i])/2./delta/mass/omega/omega;
		/*if(ir == 0)
			act += (-p[ir]-p[i])*(-p[ir]-p[i])/2./delta/mass/omega/omega;
		else
			act += (p[ir]-p[i])*(p[ir]-p[i])/2./delta/mass/omega/omega;*/
		//act += delta*sqrt(p[i]*p[i] + mass*mass);
		act += delta*sqrt((p[i]*p[i] - p0*p0)*(p[i]*p[i] - p0*p0)/4./p0/p0 + mass*mass);
	}
	return act;
}

double P[N];//h_pi_vect

void update_moment(double dtau)//propagates h_pi_vect with eps==dtau h_pi_vect-=dtau*S_der
{
	double delta = beta/N;
	for(int i = 0; i<N; i++) {
		int ir = (i+1)%N;
		int il = (N+i-1)%N;
		//P[i] -= dtau*((-p[ir]+2*p[i]-p[il])/delta/mass/omega/omega + delta*p[i]/sqrt(p[i]*p[i]+mass*mass));
		double tt = sqrt((p[i]*p[i] - p0*p0)*(p[i]*p[i] - p0*p0)/p0/p0 + 4.*mass*mass);
		P[i] -= dtau*((-p[ir]+2*p[i]-p[il])/delta/mass/omega/omega + delta*(p[i]*p[i] - p0*p0)*p[i]/p0/p0/tt);
		/*if(ir == 0)
			P[i] -= dtau*((p[ir]+2*p[i]-p[il])/delta/mass/omega/omega + delta*(p[i]*p[i] - p0*p0)*p[i]/p0/p0/tt);
		else if(il == N-1)
			P[i] -= dtau*((-p[ir]+2*p[i]+p[il])/delta/mass/omega/omega + delta*(p[i]*p[i] - p0*p0)*p[i]/p0/p0/tt);
		else
			P[i] -= dtau*((-p[ir]+2*p[i]-p[il])/delta/mass/omega/omega + delta*(p[i]*p[i] - p0*p0)*p[i]/p0/p0/tt);*/
	}
}

void update_phi(double dtau)//md molecular update for traj
{
	for(int i = 0; i<N; i++)
		p[i] += dtau*P[i];
}

void makeChange()//performs one cycle of 20 md updates and 1 langevin //!!!makes half stepes for md!!!
{
	update_moment(dtau/2);
	for(int t = 0; t<Ntau-1; t++) {
		update_phi(dtau);
		update_moment(dtau);
	}
	update_phi(dtau);
	update_moment(dtau/2);
}

bool findAccept(double expdh)//performs metrofork and tells if we should accept
{
	if(expdh >= 1.)
		return true;
	return mersenne() < expdh;
}

void dostep()//perform_sweeps equivalent
{
	for(int i = 0; i<N; i++)
		oldp[i] = p[i];//backup so can revert after metrofork
	//
	for(int i = 0; i<N; i++)
		P[i] = BoxMuller2(0., 1.);//init momentum with N(0,1)
	double oldH = action();
	for(int i = 0; i<N; i++)
		oldH += P[i]*P[i]/2.;
	//
	makeChange();
	//
	double H = action();
	for(int i = 0; i<N; i++)
		H += P[i]*P[i]/2.;
	//
	dH = H-oldH;
	bool curacc = findAccept(exp(-dH));
	if(curacc)
		acc++;
	else
		for(int i = 0; i<N; i++)
			p[i] = oldp[i];
}

double aveP()
{
	double res = 0.;
	for(int i = 0; i<N; i++)
		res += p[i]*p[i];
	return sqrt(res)/N;
}

double E1, E2, E1s, E2s;
int obs_tot;

void evalObs()
{
	double delta = beta/N;
	double E1_ = 0., E2_ = 0.;
	for(int i = 0; i<N; i++) {
		int ir = (i+1)%N;
		E1_ += (p[ir]-p[i])*(p[ir]-p[i])/2./delta/delta/mass/omega/omega;
		E2_ += sqrt((p[i]*p[i] - p0*p0)*(p[i]*p[i] - p0*p0)/4./p0/p0 + mass*mass);
	}
	E1 += E1_;
	E1s += E1_*E1_;
	E2 += E2_;
	E2s += E2_*E2_;
	obs_tot++;
}

int main()
{
	E1 = E2 = E1s = E2s = 0.;
	obs_tot = 0;
	seed_mersenne(time(NULL));
	init();
	FILE *plog = fopen("log.txt", "w");
	for(int i = 1; i<=LEN; i++) {
		dostep();
#ifdef TEST
		if(i%10==0)
			printf("%4d\t%e\t%.2f\t%.2f%%\n", i, aveP(), dH, 100.*acc/i);
#else
		if(i%100==0)
			fprintf(plog, "%4d\t%e\t%e\t%.2f\t%.2f%%\n", i, aveP(), action(), dH, 100.*acc/i);
		if(i%100000==0)
			printf("%.1fM\n", i/1000000.);
		if(i>= 2*100000 && i%10000==0) {
			//makex();
			//puthist(x);
			evalObs();
		}
#endif
	}
	fclose(plog);
	//
	FILE *pconf = fopen("conf.txt", "w");
	for(int i = 0; i<N; i++)
		fprintf(pconf, "%.14e\n", p[i]);
	fclose(pconf);
	//
#ifdef TEST
	/*double pmin = p[0], pmax = p[0];
	for(int i=1; i<N; i++) {
		if(p[i] < pmin)
			pmin = p[i];
		if(p[i] > pmax)
			pmax = p[i];
	}
	printf("p in [%e ; %e]\n", pmin, pmax);*/
#else
	/*FILE *phist = fopen("hist.txt", "w");
	for(const auto &pa : hist)
		fprintf(phist, "%f\t%f\n", hist_d*(pa.first+0.5), (double)pa.second/hist_tot);
	fclose(phist);*/
	char buf[100];
	sprintf(buf, "obs_%d.txt", N);
	FILE *pobs = fopen(buf, "w");
	double E1ave = E1/obs_tot;
	double E2ave = E2/obs_tot;
	double dE1 = sqrt(fabs(E1s/obs_tot - E1ave*E1ave));
	double dE2 = sqrt(fabs(E2s/obs_tot - E2ave*E2ave));
	fprintf(pobs, "%e +- %e\n%e +- %e\n", E1ave/N, dE1/N, E2ave/N, dE2/N);
	fclose(pobs);
	printf("mea = %d\n", obs_tot);
#endif
	return 0;
}
