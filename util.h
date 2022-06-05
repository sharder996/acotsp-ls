/*
 * \file    util.h
 * \brief   Additional functions for working with ACO TSP instances
 */


#ifndef UTIL_H
#define UTIL_H


#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* ----- FUNCTION PROTOTYPES ----- */
double euclid_dist(double x1, double y1, double x2, double y2);

int* generate_random_perm(int n);

double mean(double *nums, int len);

double std_deviation(double *nums, int len, double mean);

#endif /* UTIL_H */
