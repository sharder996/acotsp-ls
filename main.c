/*
 * \file    main.c
 * \brief   Driver program for ACO-LS
 */

#include "acotsp.h"
#include <assert.h>
#include <stdbool.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "util.h"


/* ----- CONSTANTS ----- */
#define BUFFER_LENGTH     128


/* ----- FUNCTION PROTOTYPES ----- */
void parse_args(int argc, char *argv[]);
ACO_instance* init_ACO();
ACO_instance* reinit_ACO(ACO_instance *aco);
void init_dist();
void init_pher();
void init_ants();
void destory_ACO();
void print_settings(ACO_instance *aco, FILE *stream);
void print_report(ACO_instance *aco, FILE *stream);
void print_stats(ACO_instance* aco, FILE *stream);


/* ----- STATIC VARIABLES ----- */
// ACO parameters
static int g = -1;
static int N = -1;
static double alpha = -1;
static double beta = -1;
static double rho = -1;
static double Q = -1;
static double q_0 = -1;
static double tau = -1;
static int search_iter = 0;
static int seed = -1;

static int n_iterations = 1;

// Input/output options/filenames
static char in_filename[BUFFER_LENGTH];
static FILE *out_fp = NULL;

// ACO instance
static ACO_instance *aco;


/*
 * \brief       Parse arguments
 *
 * \param[in]   argc: Number of arguments
 * \param[in]   argv: Argument vector
 */
int main(int argc, char *argv[])
{
  int i;

  parse_args(argc, argv);

  /* Initialization */
  aco = init_ACO();

  /* Seed RNG */
  if (seed < 0)
  {
    seed = time(NULL);
    srand(seed);
  }
  else
  {
    srand(seed);
  }
  /* Seed RNG */

  /* End Initialization */

  if (n_iterations > 1)
  {
    print_settings(aco, stdout);

    for (i = 0; i < n_iterations; i++)
    {
      /* Run ACO */
      run_ACO(aco);

      /* Run ACO */

      print_stats(aco, out_fp);
      reinit_ACO(aco);
    }
  }
  else
  {
    print_settings(aco, out_fp);

    /* Run ACO */
    run_ACO(aco);
    /* Run ACO */

    print_report(aco, out_fp);
  }

  destory_ACO();

  exit(EXIT_SUCCESS);
}

/*
 * \brief       Parse arguments
 *
 * \param[in]   argc: Number of arguments
 * \param[in]   argv: Argument vector
 */
void parse_args(int argc, char *argv[])
{
  char opt, *eptr;
  char out_filename[BUFFER_LENGTH];

  while((opt = getopt(argc, argv, "a:b:g:N:n:i:l:o:Q:q:r:s:t:h")) != -1)
  {
    switch (opt)
    {
      case 'a':
        alpha = strtod(optarg, &eptr);
        break;
      case 'b':
        beta = strtod(optarg, &eptr);
        break;
      case 'g':
        g = atoi(optarg);
        break;
      case 'l':
        search_iter = atoi(optarg);
        assert(search_iter > 0);
        break;
      case 'n':
        n_iterations = atoi(optarg);
        break;
      case 'N':
        N = atoi(optarg);
        break;
      case 'i':
        strncpy(in_filename, optarg, BUFFER_LENGTH);
        break;
      case 'o':
        strncpy(out_filename, optarg, BUFFER_LENGTH);
        out_fp = fopen(out_filename, "a");
        break;
      case 'Q':
        Q = strtod(optarg, &eptr);
        break;
      case 'q':
        q_0 = strtod(optarg, &eptr);
        assert(q_0 <= 1.0);
        break;
      case 'r':
        rho = strtod(optarg, &eptr);
        break;
      case 's':
        seed = atoi(optarg);
        assert(seed >= 0);
        break;
      case 't':
        tau = strtod(optarg, &eptr);
        break;
      case 'h':
      default:
        printf("Usage: %s [-i filename]\n", argv[0]);
        printf("ACO parameters:\n");
        printf("\t-a: weight of pheromone on decision\n");
        printf("\t-b: weight of heuristic data on decision\n");
        printf("\t-g: max number of generations\n");
        printf("\t-l: use local search (2-opt) number of iterations\n");
        printf("\t-N: number of ants\n");
        printf("\t-Q: amount of pheromone to be deposited\n");
        printf("\t-q: degree of random choice\n");
        printf("\t-r: rate of pheromone evaporation\n");
        printf("\t-s: seed for random choice\n");
        printf("\t-t: initial pheromone levels\n");
        printf("\t-n: number of ACO runs\n");
        printf("\t-h: help\n");
        printf("Input/Output:\n");
        printf("\t-i: specify tsp file to use as input/source file\n");
        printf("\t-o: specify file to output results to\n");

        exit(EXIT_SUCCESS);
    }
  }
}

/*
 * \brief       Create ACO instance from tsp file
 *
 * \param[in]   filename: TSP filename
 * 
 * return       Pointer to ACO_instance struct
 */
ACO_instance* init_ACO()
{
  FILE *fp;
  char buffer[BUFFER_LENGTH];
  char *endptr;
  bool header = true;
  int n_cities, i = 0;

  aco = malloc(sizeof(ACO_instance));
  fp = fopen(in_filename, "r");
  while (header && fgets(buffer, sizeof(buffer), fp))
  {
    if (strstr(buffer, "DIMENSION"))
    {
      while (!(buffer[i] >= '0' && buffer[i] <= '9'))
        i++;
      n_cities = (int)strtol(&buffer[i], &endptr, 10);
    }
    else if (strstr(buffer, "NODE_COORD_SECTION"))
    {
      fgets(buffer, sizeof(buffer), fp);
      break;
    }
  }

  aco->num_cities = n_cities;
  aco->cities = malloc(n_cities*sizeof(ACO_city));
  i = 0;

  while (i < n_cities)
  {
    (void)strtol(buffer, &endptr, 10);
    aco->cities[i].x = (int)strtol(endptr, &endptr, 10);
    aco->cities[i].y = (int)strtol(endptr, &endptr, 10);
    i++;
    fgets(buffer, sizeof(buffer), fp);
  }

  fclose(fp);

  if (out_fp == NULL)
    out_fp = stdout;

  if (N == -1)
    aco->N = n_cities;
  else
    aco->N = N;
  
  if (g == -1)
    aco->g = DEFAULT_G;
  else
    aco->g = g;

  aco->search_iter = search_iter;
  
  if (alpha == -1)
    aco->alpha = DEFAULT_ALPHA;
  else
    aco->alpha = alpha;
  
  if (beta == -1)
    aco->beta = DEFAULT_BETA;
  else
    aco->beta = beta;
  
  if (rho == -1)
    aco->rho = DEFAULT_RHO;
  else
    aco->rho = rho;

  if (Q == -1)
    aco->Q = DEFAULT_Q;
  else
    aco->Q = Q;

  if (q_0 == -1)
    aco->q_0 = DEFAULT_Q_0;
  else
    aco->q_0 = q_0;
  
  if (tau == -1)
    aco->tau = DEFAULT_TAU;
  else
    aco->tau = tau;

  aco->best = malloc(sizeof(ACO_tour));
  aco->best->tour = malloc(n_cities*sizeof(int));
  aco->best->distance = 0.0;
  aco->repeat = 0;

  aco->mean = malloc(aco->g*sizeof(double));
  aco->std_dev = malloc(aco->g*sizeof(double));
  aco->best_per_gen = malloc(aco->g*sizeof(double));
  aco->ACO_time = 0;
  aco->search_time = 0;
  aco->com_gen = 0;

  init_dist();
  init_ants();
  init_pher();

  return aco;
}

/*
 * \brief       Reinitialize ACO instance from tsp file
 *
 * \param[in]   aco: Pointer to existing aco instance
 * 
 * return       Pointer to ACO_instance struct
 */
ACO_instance* reinit_ACO(ACO_instance *aco)
{
  int i;
  int n = aco->num_cities;

  aco->best->distance = 0.0;
  aco->repeat = 0;

  aco->ACO_time = 0;
  aco->search_time = 0;
  aco->com_gen = 0;

  for (i = 0; i < n*n; i++)
  {
    aco->pheromone[i] = aco->tau;
  }

  return aco;
}

/*
 * \brief       Initiatize ACO distance matrix
 */
void init_dist()
{
  int i, j;
  int n = aco->num_cities;
  aco->distance = malloc(n*n*sizeof(double));

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      aco->distance[i*n+j] = 0.0;
      if (i != j)
      {
        aco->distance[i*n+j] = aco->distance[j*n+i] = euclid_dist(aco->cities[i].x, aco->cities[i].y, aco->cities[j].x, aco->cities[j].y);
      }
    }
  }
}

/*
 * \brief       Initiatize ACO pheromone matrix
 */
void init_pher()
{
  int i;
  int n = aco->num_cities;
  aco->pheromone = malloc(n*n*sizeof(double));

  for (i = 0; i < n*n; i++)
  {
    aco->pheromone[i] = aco->tau;
  }
}

/*
 * \brief       Initiatize ACO ants
 */
void init_ants()
{
  int i;
  int n = aco->N;
  aco->ants = malloc(n*sizeof(ACO_ant));

  for (i = 0; i < n; i++)
  {
    aco->ants[i].visited = malloc(aco->num_cities*sizeof(int));
    aco->ants[i].tour = malloc(aco->num_cities*sizeof(int));
  }
}

/*
 * \brief       Free memory allocated for ACO instance
 */
void destory_ACO()
{
  int i;

  free(aco->cities);

  for (i = 0; i < aco->N; i++)
  {
    free(aco->ants[i].visited);
    free(aco->ants[i].tour);
  }

  free(aco->ants);
  free(aco->best->tour);
  free(aco->best);
  free(aco->mean);
  free(aco->std_dev);
  free(aco->best_per_gen);

  free(aco);
}

/*
 * \brief       Prints the settings of the passed ACO instance
 *
 * \param[in]   aco: Instance of ACO
 * \param[in]   stream: Destination to print the report to
 */
void print_settings(ACO_instance *aco, FILE *stream)
{
  fprintf(stream, "----------------------------------------------\n");
  fprintf(stream, " ACO initialized with the following settings: \n");
  fprintf(stream, "----------------------------------------------\n");

  fprintf(stream, "%12s %s\n", "File: ", in_filename);
  fprintf(stream, "%12s %d\n", "Num cities: ", aco->num_cities);
  fprintf(stream, "%12s %d\n", "Num ants: ", aco->N);
  fprintf(stream, "%12s %d\n", "Max gens: ", aco->g);
  fprintf(stream, "%12s %.2f\n", "Alpha: ", aco->alpha);
  fprintf(stream, "%12s %.2f\n", "Beta: ", aco->beta);
  fprintf(stream, "%12s %.2f\n", "Rho: ", aco->rho);
  fprintf(stream, "%12s %.2f\n", "Q: ", aco->Q);
  fprintf(stream, "%12s %.2f\n", "q_0: ", aco->q_0);
  fprintf(stream, "%12s %.2f\n", "Tau: ", aco->tau);
  fprintf(stream, "%12s %d\n", "LS iters: ", aco->search_iter);
  fprintf(stream, "%12s %d\n", "Rand seed: ", seed);

  fprintf(stream, "----------------------------------------------\n\n");
}

/*
 * \brief       Prints a report of the ACO results
 *
 * \param[in]   aco: Instance of ACO
 * \param[in]   stream: Destination to print the report to
 */
void print_report(ACO_instance* aco, FILE *stream)
{
  int i;
  char *row[] = { "Gen", "Best", "Avg", "Std Dev" };

  fprintf(stream, "                    Results                   \n");
  fprintf(stream, "-----------------------------------------------\n");
  fprintf(stream, " %*s | %*s | %*s | %*s\n", -3, row[0], -12, row[1], -12, row[2], -11, row[3]);
  fprintf(stream, "-----------------------------------------------\n");
  for (i = 0; i < aco->com_gen; i++)
  {
    fprintf(stream, " %3d | %10.5f | %10.5f | %10.5f \n", i+1, aco->best_per_gen[i], aco->mean[i], aco->std_dev[i]);
  }
  fprintf(stream, "-----------------------------------------------\n");

  fprintf(stream, "Last generation of improvement: %d\n", aco->com_gen-TERM_COND);

  if (aco->search_iter)
  {
    fprintf(stream, "Local search execution %.2f (secs)\n", aco->search_time / CLOCKS_PER_SEC);
  }
  fprintf(stream, "Completed in %.2f (secs)\n\n", aco->ACO_time / CLOCKS_PER_SEC);
}

/*
 * \brief       Prints a shorter report of the ACO results
 *
 * \param[in]   aco: Instance of ACO
 * \param[in]   stream: Destination to print the report to
 */
void print_stats(ACO_instance* aco, FILE *stream)
{
  // fprintf(stream, "%10.5f, %d, %.2f\n", aco->best->distance, aco->com_gen-TERM_COND, aco->ACO_time / CLOCKS_PER_SEC);
  fprintf(stream, "%10.5f\t%d\t%.2f\n", aco->best->distance, aco->com_gen-TERM_COND, aco->ACO_time / CLOCKS_PER_SEC);
}
