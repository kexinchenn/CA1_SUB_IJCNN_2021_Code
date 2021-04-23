// corrcoef.cpp : Defines the entry point for the console application.
// remark
//#include <atltime.h>
//#include <iostream>
//#include <cstdlib>
//#include <stdlib.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <iostream>
using namespace std;

// int main(void)
// {
//         float a[9] = {1,2,3,4,5,4,3,2,1};
//         float b[9] = {1,2,3,4,5,4,3,1,1};
//         float c[9] = {5,4,3,2,1,2,3,4,5};
//         float coef = 0;
//         coef = corrCoef(a,c,9);
//         cout<<"corrCorf(a,b): "<<coef<<endl;
//         coef = corrCoef(a,b,9);
//         cout<<"corrCorf(a,c): "<<coef<<endl;
//         coef = corrCoef(a,a,9);
//         cout<<"corrCorf(a,a): "<<coef<<endl;
//         int x;
//         cin >> x;
//         return 0;
// }
//*****************************************************************/
//function name:        corrCoef
//arguments:            x[] (in) - rendom vector
//                                      y[] (in) - rendom vector, same length as x
//                                      length (in) - length of thr vectors
//return value :        result of corrCoef
//*****************************************************************/
float corrCoef(float* x,float* y, int length)
{
    const float TINY=1.0e-20;//tiny number to avoide division with zero
    float res = 0;
    int j;
    float syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
    /**first loop to find the mean of x and y ***************************/
        for (j=0;j<length;j++) 
        {
        ax += x[j];
        ay += y[j];
    }
    ax /= length;
    ay /= length;
        /**secound loop for standard deviation ******************************/
    for (j=0;j<length;j++) 
        {     
        sxx += (x[j]-ax)*(x[j]-ax);
        syy += (y[j]-ay)*(y[j]-ay);
        sxy += (x[j]-ax)*(y[j]-ay);
    }
        /********************************************************************/
    res=sxy/(sqrt(sxx*syy)+TINY);
        return res;
}
