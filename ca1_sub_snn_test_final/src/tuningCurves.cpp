#include <iostream>
#include <cmath>
#include <vector>

#include "tuningCurves.h"

using namespace std;

void linspace(float start, float end, int size, float* out_arr) {
	float step = (end - start) / (size - 1);

	for (int i =0; i < size; i ++) {
		out_arr[i] = start + step * i;
	}
}

float norm(float* in_array, int nItems) {

    float sumSq = 0;
    float output;

    for (int i = 0; i < nItems; i++) {
        sumSq += powf(in_array[i], 2.0);
    }
    output = sqrt(sumSq);
    return output;
}

vector<vector<float> > gauss2d(vector<float> orig_arr_x, vector<float> orig_arr_y, vector<float> prefs_x, vector<float> prefs_y, int numCurves, 
    int numBins, float sigma_x, float r_max) {

    vector<vector<float> > out_arr(numCurves, vector<float>(numBins, 0.0));
    float numerator;
    float denom;

    for (int bin = 0; bin < numBins; bin ++) {
        for (int neur = 0; neur < numCurves; neur++) {
            float diff[2];
            diff[0] = orig_arr_x[bin] - prefs_x[neur];
            diff[1] = orig_arr_y[bin] - prefs_y[neur];

            numerator = powf(norm(diff, 2), 2);
            denom = 2.0f * powf(sigma_x,2);

            out_arr[neur][bin] = r_max * exp(-(numerator/denom));
        }
    }
    return out_arr;
}

vector<vector<float> > gauss2dXY(vector<float> orig_arr_x, vector<float> orig_arr_y, float* prefs_x, float* prefs_y, int numCurvesX, int numCurvesY, 
    int numBins, float sigma_x, float r_max) {

    int numCurves = numCurvesX * numCurvesY;
    float numerator;
    float denom;
    vector<vector<float> > out_arr(numCurves, vector<float>(numBins, 0.0));

    for (int bin = 0; bin < numBins; bin ++) {
        for (int y = 0; y < numCurvesY; y++) {
            for (int x = 0; x < numCurvesX; x++) {
                float diff[2];

                diff[0] = orig_arr_x[bin] - prefs_x[x];
                diff[1] = orig_arr_y[bin] - prefs_y[y];

                numerator = powf(norm(diff, 2), 2);
                denom = 2.0f * powf(sigma_x,2);

                out_arr[y*numCurvesX+x][bin] = r_max * exp(-(numerator/denom));
            }
        }
    }
    return out_arr;
}

vector<vector<float> > gauss1d(vector<float> orig_arr, float* prefs, int numCurves, int numBins, float sigma, float r_max) {

    vector<vector<float> > out_arr(numCurves, vector<float>(numBins, 0.0));

    float base;
    float result;
    float exponent = 2.0;

    for (int bin = 0; bin < numBins; bin ++){
        for (int neur = 0; neur < numCurves; neur++) {
            base = (orig_arr[bin] - prefs[neur])/sigma;
            result = powf(base, exponent);
            out_arr[neur][bin] = r_max * exp(-(0.5) * result);
        }
    }
    return out_arr;
}

vector<vector<float> > cosine(vector<float> orig_arr, float* prefs, int numCurves, int numBins, float r_max, float r_min) {

    vector<vector<float> > out_arr(numCurves, vector<float>(numBins, 0.0));

    float fr;
    for (int bin = 0; bin < numBins; bin ++){
        for (int neur = 0; neur < numCurves; neur++) {
            // fr = r_min + (r_max - r_min) * cos(orig_arr[bin] - prefs[neur]);
            // if (fr < 0) {
            //     fr = 0;
            // }
            if ( (orig_arr[bin] >= (prefs[neur] - 0.7854) ) &&
                (orig_arr[bin] <= (prefs[neur] + 0.7854) ) ) {
                fr = r_min + (r_max - r_min) * cos(orig_arr[bin] - prefs[neur]);
            } else {
                fr = 0;
            }
            out_arr[neur][bin] = fr;
        }
    }
    return out_arr;
}