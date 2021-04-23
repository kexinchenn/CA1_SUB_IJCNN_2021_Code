/*
filtfilt.cpp
A C++ implementation of MATLAB filtfilt.m.
Written by Kris Carlson @ UCI.
11-04-2014
Adapted from:
FILTER.C
An ANSI C implementation of MATLAB FILTER.M (built-in)
Written by Chen Yangquan <elecyq@nus.edu.sg>
1998-11-11.
*/
#include <cstdio>
#include <cstdlib>


void filter(int ord, float *a, float *b, int np, float *x, float *y);
// filtfilt will call the filter function twice
void filtfilt(int ord, float *a, float *b, int np, float *x, float *y, int size);
void custom(float *x, float *y, int size);
void getFiltered(int ord, float *a, float *b, int np, float *x, float *y);