/*
 * \file    acotsp.c
 * \brief   ACO for TSP implementation
 */


#include "acotsp.h"


/* ----- FUNCTION PROTOTYPES ----- */
void reset_ant(ACO_ant* ant, int n, int pos);
int* generate_random_perm(int n);
int next_city(ACO_instance* aco, int ant_index);
double mean_dist(ACO_instance *aco);
double std_deviation_dist(ACO_instance *aco, double mean);


/*
 * \brief       Run ACO based on the parameters given in the passed ACO instance
 *
 * \param[in]   aco: Instance of ACO
 */
void run_ACO(ACO_instance *aco)
{
  int i, j;
  clock_t start_time, end_time, s_start_time, s_end_time;
  double prev_dist = 0.;

  start_time = clock();

  reset_ants(aco);

  for (i = 0; i < aco->g; i++)
  {
    for (j = 0; j < aco->num_cities; j++)
    {
      step_ants(aco);
    }


    // Apply local search
    if (aco->search_iter)
    {
      s_start_time = clock();
      for (j = 0; j < aco->N; j++)
        two_opt_search(aco, j);
      s_end_time = clock();
      aco->search_time += s_end_time - s_start_time;
    }
    
    
    update_pheromone(aco);
    update_best(aco);

    // Keep track of number of generations with no improvement
    if (fabs(prev_dist - aco->best->distance) < 0.00001)
      aco->repeat++;
    else
      aco->repeat = 0;
    prev_dist = aco->best->distance;
    

    // Record stats
    aco->best_per_gen[i] = aco->best->distance;
    aco->com_gen++;
    aco->mean[i] = mean_dist(aco);
    aco->std_dev[i] = std_deviation_dist(aco, aco->mean[i]);

    // Exit on termination condition
    if (aco->repeat >= TERM_COND)
      break;
    
    reset_ants(aco);
  }

  end_time = clock();
  aco->ACO_time = end_time - start_time;
}

/*
 * \brief       Calculate total distance of ant's tour
 *
 * \param[in]   aco: ACO instance
 * \param[in]   ant_index: Index of ant to calculate distance of
 * 
 * return       Euclidean distance of tour
 */
double calc_dist(ACO_instance* aco, int ant_index)
{
  int i, n = aco->num_cities;
  double sum = 0;

  for (i = 0; i < n - 1; i++)
  {
    sum += aco->distance[aco->ants[ant_index].tour[i]*n+aco->ants[ant_index].tour[i+1]];
  }
  sum += aco->distance[aco->ants[ant_index].tour[n-1]*n+aco->ants[ant_index].tour[0]];

  return sum;
}

/*
 * \brief       Validate a tour
 *
 * \param[in]   ant: Pointer to ant
 * \param[in]   n: Number of cities in tour
 * 
 * return       Boolean value of valid tour
 */
bool valid_tour(ACO_ant* ant, int n)
{
  int i;
  bool *seen, retval = true;

  seen = malloc(n * sizeof(bool));
  for (i = 0; i < n; i++)
  {
    seen[i] = false;
  }

  for (i = 0; i < n && retval; i++)
  {
    retval = !seen[ant->tour[i]];
    seen[ant->tour[i]] = true;
  }

  free(seen);
  return retval;
}

/*
 * \brief       Reset the position and previous tour of each ant
 *
 * \param[in]   aco: Pointer to ACO instance
 */
void reset_ants(ACO_instance* aco)
{
  int i;

  for (i = 0; i < aco->N; i++)
  {
    reset_ant(&aco->ants[i], aco->num_cities, i);
  }
}

/*
 * \brief       Helper method that resets a single ant
 *
 * \param[in]   aco: Pointer to ant
 * \param[in]   n: Number of cities in tour
 * \param[in]   i: Initial position of ant
 */
void reset_ant(ACO_ant* ant, int n, int pos)
{
  int i;

  ant->city = rand() % n;
  ant->path_index = 1;
  ant->tour_distance = 0.0;

  for (i = 0; i < n; i++)
  {
    ant->visited[i] = 0;
    ant->tour[i] = -1;
  }
  ant->visited[ant->city] = 1;
  ant->tour[0] = ant->city;
}

/*
 * \brief       Advance each ant to next city in its tour
 *
 * \param[in]   aco: Pointer to ACO instance
 */
void step_ants(ACO_instance* aco)
{
  int i, n = aco->num_cities;

  for (i = 0; i < aco->N; i++)
  {
    if (aco->ants[i].path_index < n)
    {
      aco->ants[i].next_city = next_city(aco, i);
      aco->ants[i].tour_distance += aco->distance[aco->ants[i].city*n+aco->ants[i].next_city];
      aco->ants[i].tour[aco->ants[i].path_index++] = aco->ants[i].next_city;
      aco->ants[i].visited[aco->ants[i].next_city] = 1;
      aco->ants[i].city = aco->ants[i].next_city;

      if (aco->ants[i].path_index == n)
      {
        aco->ants[i].tour_distance += aco->distance[aco->ants[i].tour[n-1]*n+aco->ants[i].tour[0]];
      }
    }
  }
}

/*
 * \brief       Choose the next city in a tour
 *
 * \param[in]   aco: Pointer to aco instance
 * \param[in]   ant_index: Index of ant moving to next city
 * 
 */
int next_city(ACO_instance* aco, int ant_index)
{
  int i;
  double r, sum_prob = 0.0, part_sum = 0.;
  
  r = (double)rand() / (double)RAND_MAX * aco->q_0;

  for (i = 0; i < aco->num_cities; i++)
  {
    if (!aco->ants[ant_index].visited[i])
    {
      sum_prob += prob_product(aco, aco->ants[ant_index].city, i);
    }
  }

  r *= sum_prob;
  for (i = 0; i < aco->num_cities; i++)
  {
    if (!aco->ants[ant_index].visited[i])
    {
      part_sum += prob_product(aco, aco->ants[ant_index].city, i);
      if (part_sum > r)
      {
        break;
      }
    }
  }

  return i;
}

/*
 * \brief       Update pheromone matrix from new tours
 *
 * \param[in]   aco: Pointer to aco instance
 */
void update_pheromone(ACO_instance *aco)
{
  int i, j, from, to, n = aco->num_cities;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (i != j)
      {
        aco->pheromone[i*n+j] *= 1.0 - aco->rho;
        if (aco->pheromone[i*n+j] < 0.0)
        {
          aco->pheromone[i*n+j] = 1.0 / n;
        }
      }
    }
  }
  
  for (i = 0; i < aco->N; i++)
  {
    for (j = 0; j < n; j++)
    {
      from = aco->ants[i].tour[j];

      if (j < n - 1)
      {
        to = aco->ants[i].tour[j+1];
      }
      else
      {
        to = aco->ants[i].tour[0];
      }

      // aco->pheromone[from*n+to] += aco->Q / aco->ants[i]->tour_distance;
      aco->pheromone[from*n+to] += aco->Q;
      aco->pheromone[to*n+from] = aco->pheromone[from*n+to];
    }
  }
}

/*
 * \brief       Calculate the pheromone-distance product
 *
 * \param[in]   aco: Pointer to instance of ACO
 * \param[in]   from: Index of starting city
 * \param[in]   y1: Index of destination city
 * 
 * return       Product
 */
double prob_product(ACO_instance *aco, int from, int to)
{
  int n = aco->num_cities;
  double a = pow(aco->pheromone[from*n+to], aco->alpha);
  double b = pow((1.0 / (aco->distance[from*n+to] + 0.1)), aco->beta);
  return a * b;
}

/*
 * \brief       Compare tours and save the overall best one
 *
 * \param[in]   aco: Pointer to aco instance
 */
void update_best(ACO_instance *aco)
{
  int i, j;
  for (i = 0; i < aco->N; i++)
  {
    if (aco->ants[i].tour_distance < aco->best->distance || aco->best->distance == 0.0)
    {
      aco->best->distance = aco->ants[i].tour_distance;
      for (j = 0; j < aco->num_cities; j++)
      {
        aco->best->tour[j] = aco->ants[i].tour[j];
      }
    }
  }
}

/*
 * \brief       Use 2-opt search to search for better tour
 *
 * \param[in]   aco: Pointer to instance of aco
 * \param[in]   ant_index: Index of ant to perform search on
 */
void two_opt_search(ACO_instance *aco, int ant_index)
{
  int i, j, k, l, temp;
  double new_dist;
  int c1, s_c1, c2, s_c2, p_c1;
  int *rnd;
  int *pos;
  int n_improves = 0;
  int n = aco->num_cities;


  rnd = generate_random_perm(n);
  pos = malloc(n * sizeof(int));
  for (i = 0; i < n; i++)
  {
    pos[aco->ants[ant_index].tour[i]] = i;
  }

  for (i = 0; i < aco->search_iter; i++)
  {
    c1 = rnd[i];
    p_c1 = aco->ants[ant_index].tour[(pos[c1]+1)%n];

    for (j = 0; j < n; j++)
    {
      s_c1 = aco->ants[ant_index].tour[(pos[c1]+n-1)%n];
      c2 = rand() % (n-1);
      s_c2 = aco->ants[ant_index].tour[(pos[c2]+1)%n];
      if (c2 == c1 || c2 == s_c1 || c2 == p_c1)
      {
        continue;
      }
      
      new_dist = aco->ants[ant_index].tour_distance - aco->distance[c1*n+s_c1] - aco->distance[c2*n+s_c2];
      new_dist += aco->distance[s_c1*n+c2] + aco->distance[c1*n+s_c2];

      if ((int)new_dist < (int)aco->ants[ant_index].tour_distance)
      {
        k = pos[c1]; l = pos[c2];
        if (k > l)
        {
          temp = k; k = l+1; l = temp-1;
        }
        
        while (k < l)
        {
          temp = aco->ants[ant_index].tour[k];
          pos[temp] = l;
          pos[aco->ants[ant_index].tour[l]] = k;
          aco->ants[ant_index].tour[k] = aco->ants[ant_index].tour[l];
          aco->ants[ant_index].tour[l] = temp;
          k++; l--;
        }

        aco->ants[ant_index].tour_distance = calc_dist(aco, ant_index);
        n_improves++;
      }
    }
  }
  
  free(rnd);
  free(pos);
}

/*
 * \brief       Calculate mean of values in array
 *
 * \param[in]   aco: ACO Instance
 * 
 * return       Mean of values
 */
double mean_dist(ACO_instance *aco)
{
  int i;
  double mean = 0;

  for (i = 0; i < aco->N; i++)
  {
    mean += aco->ants[i].tour_distance;
  }
  mean /= aco->N;

  return mean;
}

/*
 * \brief       Calculate the standard deviation from an array of values
 *
 * \param[in]   aco: ACO instance
 * \param[in]   mean: Mean of values in array
 * 
 * return       Standard deviation of values
 */
double std_deviation_dist(ACO_instance *aco, double mean)
{
  int i;
  double std_dev = 0;

  for (i = 0; i < aco->N; i++)
  {
    std_dev += (aco->ants[i].tour_distance - mean) * (aco->ants[i].tour_distance - mean);
  }

  return sqrt(std_dev / (aco->N-1));
}
