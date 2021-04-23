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
#include <filtfilt.h>
#include <iostream>

// filter function does half of filtfilt
void filter(int ord, float *a, float *b, int np, float *x, float *y)
{
        int i,j;
	y[0]=b[0]*x[0];
	for (i=1;i<ord+1;i++)
	{
        y[i]=0.0;
        for (j=0;j<i+1;j++)
        	y[i]=y[i]+b[j]*x[i-j];
        for (j=0;j<i;j++)
        	y[i]=y[i]-a[j+1]*y[i-j-1];
	}
/* end of initial part */
for (i=ord+1;i<np+1;i++)
{
	y[i]=0.0;
        for (j=0;j<ord+1;j++)
        y[i]=y[i]+b[j]*x[i-j];
        for (j=0;j<ord;j++)
        y[i]=y[i]-a[j+1]*y[i-j-1];
}

} /* end of filter */


// filtfilt will call the filter function twice
void filtfilt(int ord, float *a, float *b, int np, float *x, float *y, int nBins)
{
	float initial_x[nBins];

	for (int i = 0; i < nBins; i++) {
		initial_x[i] = x[i];
	}

	int padding = 40;
	int initial_np = np + padding * 2;
	int start = padding+1;
	int stop = padding+np;

	np = np+padding*2;

	float padded_x[np];

	float padded_y[np];

	for (int i = 0; i < np; i++) {
		padded_x[i] = 1.0;
	}

	for (int i = 0; i < nBins; i++) {
		padded_x[start+i] = initial_x[i];
	}

	// set leading and trailing padding values to first and last input values respectively
	// (reduces error on edges)
	for (int i = 0; i <= (start-1); i++) {
		padded_x[i] = padded_x[i] * initial_x[0];
	}
	for (int i = stop; i < nBins; i++) {
		padded_x[i] = padded_x[i] * initial_x[(initial_np-1)];
	}
	
	// first filter pass
	getFiltered(ord,a,b,np,padded_x,padded_y);
	
	// second filter pass
	getFiltered(ord,a,b,np,padded_x,padded_y);

	for (int i = 0; i < nBins; i++) {
		y[i] = padded_x[start+i];
	}
} /* end of filtfilt */

void custom(float *x, float *y, int nBins) {
	int ord = 13;
	float a[14] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float b[14] = {6.07588281736937e-09, 1.48671950679309e-06, 0.000133830225050040, 0.00443184838826559, 
		0.0539909662247990, 0.241970723226673, 0.398942278270510, 0.241970723226673, 0.0539909662247990,
		0.00443184838826559, 0.000133830225050040, 1.48671950679309e-06, 6.07588281736937e-09, 9.13472035957211e-12};

	filter(ord,a,b,nBins,x,y);
	filtfilt(ord, a, b, nBins, x, y, nBins);

}

void getFiltered(int ord, float *a, float *b, int np, float *x, float *y) {
	//float temp_y[np];

	filter(ord,a,b,np,x,y);

	// for (int i = 1; i < (np+1); i++) {
	// 	temp_y[i-1] = y[i];
	// }

	// for (int i = 0; i < np; i++) {
	// 	y[i] = temp_y[i];
	// }

	// flip array
	for (int i = 0; i < np; i++) {
		x[i] = y[np-(i+1)];
	}
}