/*
 * \file    acotsp.h
 * \brief   ACO for TSP implementation
 */


#ifndef ACOTSP_H
#define ACOTSP_H


#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "util.h"


/* ----- CONSTANTS ----- */
#define DEFAULT_G         1000
#define DEFAULT_ALPHA     0.5
#define DEFAULT_BETA      0.5
#define DEFAULT_RHO       0.95
#define DEFAULT_Q         2
#define DEFAULT_Q_0       0.5
#define DEFAULT_TAU       0.5
#define TERM_COND         25


/* ----- TYPE DEFINITIONS ----- */
typedef struct 
{
  double x;                    /* x-coordinate of city */
  double y;                    /* y-coordinate of city */
} ACO_city;

typedef struct 
{
  double distance;          /* Distance of tour */
  int *tour;                /* Ordering of path of tour */
} ACO_tour;

typedef struct 
{
  int city;                 /* Current position of ant */
  int next_city;            /* Next city in ant's tour */
  int *visited;             /* List of seen cities */
  int *tour;                /* Path of ant's tour */
  int path_index;           /* Current index in tour */
  double tour_distance;     /* Distance of tour */
} ACO_ant;

typedef struct 
{
  /* ACO data structures */
  ACO_city *cities;         /* List of cities */
  ACO_ant *ants;            /* List of ants */
  ACO_tour *best;           /* Best tour found */
  int num_cities;           /* Number of cities in space */
  double *distance;         /* Distance matrix */
  double *pheromone;        /* Pheromone matrix */

  /* Settings */
  int search_iter;          /* Flag for using local search */
  int N;                    /* Number of ants */
  int g;                    /* Number of generations */
  double alpha;             /* Weight of pheromone */
  double beta;              /* Weight of heuristic data */
  double rho;               /* Rate of pheromone evaporation */
  double Q;                 /* Amount of pheromone to be deposited */
  double q_0;               /* Degree of random choice */
  double tau;               /* Initial amount of pheromone */
  int repeat;               /* Number of generations with no change */

  /* Stats */
  int com_gen;              /* Number of completed generations */
  double *mean;             /* Mean length of tours per generation */
  double *std_dev;          /* Standard deviation of length of tours per generation */
  double *best_per_gen;     /* Shortest tour per generation */
  double ACO_time;          /* Time spent on ACO algorithm */
  double search_time;       /* Time spent on local search */
} ACO_instance;


/* ----- FUNCTION PROTOTYPES ----- */
void run_ACO(ACO_instance *aco);

void reset_ants(ACO_instance* aco);

void step_ants(ACO_instance* aco);

void update_pheromone(ACO_instance *aco);

void update_best(ACO_instance *aco);

double prob_product(ACO_instance *aco, int from, int to);

double calc_dist(ACO_instance* aco, int ant_index);

bool valid_tour(ACO_ant* ant, int n);

void two_opt_search(ACO_instance *aco, int ant_index);

#endif /* ACOTSP_H */
