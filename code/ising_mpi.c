#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>

void init ( int ** F, int L, gsl_rng * r ) {
  int i,j;

  /* Visit each position (i,j) in the lattice */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      /* 2*x-1, where x is randomly 0,1 */
      F[i][j]=2*(int)gsl_rng_uniform_int(r,2)-1;
    }
  }
}

int main(int argc, char** argv)
{

  int ** F;       /* The 2D array of spins */
  int L = 1000;     /* The sidelength of the array */
  int N;          /* The total number of spins = L*L */
  double T = 1.0; /* Dimensionless temperature = (T*k)/J */

  /* Run parameters */
  int nCycles = 1000000; /* number of MC cycles to run; one cycle is N 
			    consecutive attempted spin flips */
  int fSamp = 1000;      /* Frequency with which samples are taken */

  /* Computational variables */
  int nSamp;      /* Number of samples taken */
  int de;         /* energy change due to flipping a spin */
  double b;       /* Boltzman factor */
  double x;       /* random number */
  int i,j,a,c;    /* loop counters */

  /* Observables */
  double s=0.0, ssum=0.0;    /* average magnetization */
  double e=0.0, esum=0.0;    /* average energy per spin */

  /* This line creates a random number generator
     of the "Mersenne Twister" type, which is much
     better than the default random number generator. */
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-L")) L=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-s")) Seed = (unsigned long)atoi(argv[++i]);
  }
  
  /* Output some initial information */
  /*fprintf(stdout,"# command: ");
  for (i=0;i<argc;i++) fprintf(stdout,"%s ",argv[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"# ISING simulation, NVT Metropolis Monte Carlo\n");
  fprintf(stdout,"# L = %i, T = %.3lf, nCycles %i, fSamp %i, Seed %u\n",
	  L,T,nCycles,fSamp,Seed);*/

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* Compute the number of spins */
  N=L*L;

  /* Allocate memory for the system */
  F=(int**)malloc(L*sizeof(int*));
  for (i=0;i<L;i++) F[i]=(int*)malloc(L*sizeof(int));

  /* Generate an initial state */
  init(F,L,r);

  /* For computational efficiency, convert T to reciprocal T */
  T=1.0/T;

  s = 0.0;
  e = 0.0;
  nSamp = 0;

  int ierr;
  int size, rank;
  int dims[] = {2,2};
  int periodic[] = {1,1};
  int coords[2];

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm cartcomm;

  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 0, &cartcomm);
  ierr = MPI_Comm_rank(cartcomm, &rank);
  ierr = MPI_Cart_coords(cartcomm, rank, 2, coords);

  fprintf("rank: %d, coords: %d,%d", rank, coords[0], coords[1]);



  return 0;
}