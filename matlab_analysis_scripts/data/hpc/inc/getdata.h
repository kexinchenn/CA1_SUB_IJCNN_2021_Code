
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <sstream>
#include <vector>

using namespace std;


vector<vector<float> > getXYData (int, int, bool, int, int);

void getRealTestFRs(int, int, float**);

vector<float> getAngVel (int, int, bool, int, int);

vector<float> getLinVel (int, int, bool, int, int);

vector<float> getHeadDir (int, int, bool, int, int);

vector<float> getOcc (int rec, int posTrajType, bool testing, int trial, int bins);

void getCoordPrefs (float*, float*, int);
void getProgPrefs(float*, int);
void getLVPrefs(float*, int);
void getAVPrefs(float*, int);

// int getTrTrials(int, int);
// int getTeTrials(int, int);

vector<vector<float> > loadData(std::string, int, int);

string getFileName(string dataDir, string inputVar, int path);
