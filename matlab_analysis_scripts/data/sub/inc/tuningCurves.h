#include <cmath>
#include <vector>

using namespace std;

void linspace(float start, float end, int size, float* out_arr);

float norm(float * a, int n);

// vector<vector<float> > gauss2d(vector<float> orig_arr_x, vector<float> orig_arr_y, float* prefs_x, float* prefs_y, int numCurves, 
//     int numBins, float sigma_x, float r_max);
vector<vector<float> > gauss2d(vector<float> orig_arr_x, vector<float> orig_arr_y, float* prefs_x, float* prefs_y, int numCurvesX, int numCurvesY, 
    int numBins, float sigma_x, float r_max);
vector<vector<float> > gauss1d(vector<float> orig_arr, float* prefs, int numCurves, int numBins, float sigma, float r_max);
vector<vector<float> >cosine(vector<float> orig_arr, float* prefs, int numCurves, int numBins, float r_max, float r_min);