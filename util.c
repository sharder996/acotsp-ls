/*
 * \file    util.c
 * \brief   Additional functions for working with ACO TSP instances
 */


#include "util.h"



/*
 * \brief       Calculate Euclidean distance between two points on a 2D plane
 *
 * \param[in]   x1: x-coordinate of first point
 * \param[in]   y1: y-coordinate of first point
 * \param[in]   x2: x-coordinate of second point
 * \param[in]   y2: y-coordinate of second point
 * 
 * return       Euclidean distance between points
 */
double euclid_dist(double x1, double y1, double x2, double y2)
{
  return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}

/*
 * \brief       Create an array containing a random permutation from 0 to n-1
 *
 * \param[in]   n: Size of array
 * 
 */
int* generate_random_perm(int n)
{
  int i, node, temp, tot_assigned = 0;
  double rnd;
  int *r;

  r = malloc(n * sizeof(int));

  for (i = 0; i < n; i++)
  {
    r[i] = i;
  }

  for (i = 0; i < n; i++)
  {
    rnd = (double)rand() / (double)RAND_MAX;
    node = (int)(rnd * (n - tot_assigned));
    temp = r[i];
    r[i] = r[i+node];
    r[i+node] = temp;
    tot_assigned++;
  }

  return r;
}

/*
 * \brief       Calculate mean of values in array
 *
 * \param[in]   nums: array of values
 * \param[in]   len: length of array
 * 
 * return       Mean of values
 */
double mean(double *nums, int len)
{
  int i;
  double mean = 0;

  for (i = 0; i < len; i++)
  {
    mean += nums[i];
  }
  mean /= len;

  return mean;
}

/*
 * \brief       Calculate the standard deviation from an array of values
 *
 * \param[in]   nums: array of values
 * \param[in]   len: length of array
 * \param[in]   mean: mean of values in array
 * 
 * return       Standard deviation of values
 */
double std_deviation(double *nums, int len, double mean)
{
  int i;
  double std_dev = 0;

  for (i = 0; i < len; i++)
  {
    std_dev += (nums[i] - mean) * (nums[i] - mean);
  }

  return sqrt(std_dev / (len-1));
}
