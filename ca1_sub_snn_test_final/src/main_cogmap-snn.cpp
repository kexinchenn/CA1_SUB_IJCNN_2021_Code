/*! \brief 
 *
 */
#include <carlsim.h>
#include "PTI.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
// #include <string.h>
#include <vector>
 #include <limits>
#include <limits.h>
#include <time.h>
#include <algorithm>
#include <fstream>
#include <bits/stdc++.h> 
#include <sys/stat.h> 
#include <sys/types.h> 

#include "filtfilt.h" // C++ implementation of Matlab filtfilt function
#include "getdata.h"
#include "corrcoef.h"
 #include "getCombination.h"
// #include <covariance.h>
#include "tuningCurves.h"

using namespace std;

class CogmapExperiment : public Experiment {
public:

	const LoggerMode verbosity;
	const SimMode simMode;
	CogmapExperiment(const SimMode simMode, const LoggerMode verbosity): simMode(simMode), verbosity(verbosity) {}

	void run(const ParameterInstances &parameters, ostream &outputStream) const {

		// time_t seed = time(NULL);
		// time_t seed = 1611767486; // CA1_1_indi
		time_t seed = 1611539389; // SUB_one
		// time_t seed = 1611097203; // SUB_conn_prob_5_indi
		// time_t seed = 1609285682; // CA1_pos_all_2_indi
		// time_t seed = 1609618455; // SUB_pos_all_2_indi
		// time_t seed = 1609864806; // SUB_pos_all_3_indi
		// time_t seed = 1609896419; // CA1_pos_all_3_indi
		// srand(seed);

		bool loadFlag = 1;
		bool writeSimData = 1;
		bool forceMatch = 1;
		bool inputOneGroup = 1;

		string net_folder = "/data/cogmap_ca1_sub/SUB_5_indi/";
		if (loadFlag == 0) {
			if (writeSimData) {
				if (mkdir(net_folder.c_str(), 0777) == -1) {
					 cerr << "Error :  " << strerror(errno) << endl; 
					 return;
				}
			}
		}
		string sim_name = net_folder+"sim_direct_driver_sub.dat";

		// bool LESION = 0;

		bool lesionAV = 1;
		bool lesionHD = 0;
		bool lesionPos = 1;
		bool lesionLV = 1;
		bool lesionExc = 0;
		bool lesionInh = 0;

		bool testing;

		// int i = 7;
		// int run = 10;
		string dataRoot = "./data/";

		time_t timerTrial_start,timerTrial_end;
		time_t globalTimer_start, globalTimer_end;

		time(&globalTimer_start);

		// DATA LOADING VARIABLES
		const int totalTrTrials = 20;
		const int totalTeTrials = 20;
		const int repTrialsTr = 5;
		const int repTrialsTe = 5;
		// const int totalTrTrials = 2;
		// const int totalTeTrials = 2;
		// const int repTrialsTr = 2;
		// const int repTrialsTe = 2;

		// SUB
		const int numNeuronsAll = 382;
		const int numNeurons = 382;
		const int numRecordings = 49;
		// // HPC
		// const int numNeuronsAll = 322;
		// const int numNeurons = 295;
		// const int numRecordings = 32;

		const int numPaths = 6;

		int r_min = 0;

		int INP_ANG = 12;
		int INP_HD = 8;
		int INP_POS_X = 25;
		int INP_POS_Y = 18;
		int INP_POS = INP_POS_X * INP_POS_Y;
		// int INP_POS = 394;
		// int INP_POS = 398;
		int INP_LV = 12;
		int INP_N = INP_ANG+INP_HD+INP_POS+INP_LV;

		int INP_POS_START = INP_ANG+INP_HD;
		int INP_POS_END = INP_ANG+INP_HD+INP_POS;

		int BINS_IN = 140;
		int BINS_OUT = 197;
		int NUM_IN = 4;
		int NUM_OUT = 2;
		int BINS_ALL = BINS_IN * NUM_IN + BINS_OUT * NUM_OUT;

		float base;
		float result;
		float r_max = 60;
		float r_max_ang_lin = 40;
		// float coef_thresh = .5;
		float maxFR_thresh = 100.0;

		//POSITION ON TRACK PARAMETERS:
		float xCoordPref[INP_POS_X];
		float yCoordPref[INP_POS_Y];
		// vector<float> xCoordPref(INP_POS, 0.0);
		// vector<float> yCoordPref(INP_POS, 0.0);

		// float sigma_x = 15.0;
		float sigma_x = 25.0;

		float sigmaAngVel = 6.0;
		float angVelPref[INP_ANG];

		float sigmaLinVel = 20.0;
		float linVelPref[INP_LV];
		
		// float headDirPref[INP_HD];
		float headDirPref[INP_HD] = {-2.3562, -1.5708, -0.7854, 0, 0.7854, 1.5708, 2.3562, 3.1416};

		// get preferred velocity, direction, pos
		const float pi = 3.14;
		linspace(-40, 40, INP_ANG, angVelPref);
		linspace(0, 200, INP_LV, linVelPref);
		// linspace(-pi, pi, INP_HD, headDirPref);
		linspace(0, 700, INP_POS_X, xCoordPref);
		linspace(0, 500, INP_POS_Y, yCoordPref);

		// read in trial count
		string dataDir = (dataRoot + "data_SUB_train_test_20201216/");
		// string dataDir = (dataRoot + "data_HPC_train_test_20201204/");

		// number of trials for each path in each recording (numRecording x numPath)
		string trialCountFileTrain = (dataDir + "trial_count_train.csv");
		vector<vector<float> > trialCountTrain; 
		trialCountTrain = loadData(trialCountFileTrain, numRecordings, numPaths);

		string trialCountFileTest = (dataDir + "trial_count_test.csv");
		vector<vector<float> > trialCountTest; 
		trialCountTest = loadData(trialCountFileTest, numRecordings, numPaths);				

		// // template XY position to create input place cells (2 x numPosXY)
		// string templatePosXYFile = (dataDir + "allPathTemplateXYDownsampledBy5.csv");
		// // string templatePosXYFile = (dataDir + "pathTemplateXYDownsampledBy4.csv");
		// vector<vector<float> > templatePosXY;
		// templatePosXY = loadData(templatePosXYFile, 2, INP_POS);

		// ofstream templateFile;
		// templateFile.open(net_folder+"templateXY.csv", ofstream::out | ofstream::app);

		// xCoordPref = templatePosXY[0];
		// yCoordPref = templatePosXY[1];	
		// for (int xyInd = 0; xyInd < INP_POS; xyInd++) {
		// 	xCoordPref[xyInd] = templatePosXY[0][xyInd];
		// 	yCoordPref[xyInd] = templatePosXY[1][xyInd];
		// 	templateFile << xCoordPref[xyInd] << ",";
			
		// }
		// templateFile << endl;
		// for (int xyInd = 0; xyInd < INP_POS; xyInd++) {
		// 	templateFile << yCoordPref[xyInd] << ",";
		// }
		// templateFile << endl;
		/* ----------------------------- READ IN REAL NEURON DATA ------------------------------------------*/

		float **trueTeFR;
		trueTeFR = new float*[numNeurons];
		for (int n = 0; n < numNeurons; n++){
			trueTeFR[n] = new float[BINS_ALL];
		}
		
		getRealTestFRs(dataDir, numNeurons, numNeuronsAll, numPaths, trueTeFR);

		// -------------------------------- DONE READING IN ALL DATA ------------------------------------------------

		// CARLSIM VARIABLES:
		int NUM_N = 800;
		int numExc = 640;
		int numInh = 160;

		// Izhikevich parameters
		const float REG_IZH[] = { 0.02f, 0.2f, -65.0f, 8.0f };
		const float FAST_IZH[] = { 0.1f, 0.2f, -65.0f, 2.0f };
		const float COND_tAMPA=5.0, COND_tNMDA=150.0, COND_tGABAa=6.0, COND_tGABAb=150.0;
		#define HOMEO_FACTOR (0.1) //original
		// #define HOMEO_AVERAGE_TIME_SCALE (1.0) // original: begin with this value from Jay's
		float HOMEO_AVERAGE_TIME_SCALE_EXC;
		float HOMEO_AVERAGE_TIME_SCALE_INH;

		const int indiNum = parameters.getNumInstances();
		// if (indiNum > 1) {
		// 	writeSimData = 0;
		// }
		int runTimeSec = 0;
		int runTimeMs = 15;

		float fitness[indiNum];

		int exc[indiNum];
		int inh[indiNum];

		SpikeMonitor* excMonitor[indiNum];
		SpikeMonitor* inhMonitor[indiNum];

		ConnectionMonitor* inpToExcCM[indiNum];
		ConnectionMonitor* inpToInhCM[indiNum];
		ConnectionMonitor* excToExcCM[indiNum];
		ConnectionMonitor* inhToExcCM[indiNum];

		ConnectionMonitor* inpAVToExcCM[indiNum];
		ConnectionMonitor* inpHDToExcCM[indiNum];
		ConnectionMonitor* inpPosToExcCM[indiNum];
		ConnectionMonitor* inpLVToExcCM[indiNum];

		// EXPERIMENT VARIABLES: 
		int trialNumTr;
		int recNum;

		int exc2D;
		int inh2D;
		int inp2D;

		int count = 0;
		
		// 4d array to store smoothed exc FR (ind-rec-trial-nExc)
		float *****smoothedFRs = new float****[indiNum];
		for (int i = 0; i < indiNum; i ++) {
			smoothedFRs[i] = new float***[totalTeTrials];
			for (int r = 0; r < totalTeTrials; r++) {
				smoothedFRs[i][r] = new float**[repTrialsTe];
				for (int j = 0; j < repTrialsTe; j++) {
					smoothedFRs[i][r][j] = new float*[numExc];
					for (int n = 0; n < numExc; n++) {
						smoothedFRs[i][r][j][n] = new float[BINS_ALL];
						for (int b = 0; b < BINS_ALL; b++) {
							smoothedFRs[i][r][j][n][b] = 0.0;
						}
					}
				}
			}
		}
		
		// 4d array to store inp data (rec-trial-nTypesInp)
		float ****inpsDataAll = new float***[totalTeTrials];
		for (int r = 0; r < totalTeTrials; r++) {
			inpsDataAll[r] = new float**[repTrialsTe];
			for (int j = 0; j < repTrialsTe; j++) {
				inpsDataAll[r][j] = new float*[5];
				for (int n = 0; n < 5; n++) {
					inpsDataAll[r][j][n] = new float[BINS_ALL];
					for (int b = 0; b < BINS_ALL; b++) {
						inpsDataAll[r][j][n][b] = 0.0;
					}
				}
			}
		}

		// 4d array to store inp FR (rec-trial-nInpNeuron)
		float ****inpsSmthFRs = new float***[totalTeTrials];
		for (int r = 0; r < totalTeTrials; r++) {
			inpsSmthFRs[r] = new float**[repTrialsTe];
			for (int j = 0; j < repTrialsTe; j++) {
				inpsSmthFRs[r][j] = new float*[INP_N];
				for (int n = 0; n < INP_N; n++) {
					inpsSmthFRs[r][j][n] = new float[BINS_ALL];
					for (int b = 0; b < BINS_ALL; b++) {
						inpsSmthFRs[r][j][n][b] = 0.0;
					}
				}
			}
		}

		// 2d array to store mean correlation for each neuron pair (nNeur-nExc)
		float **avgCoefficient;
		avgCoefficient = new float*[numNeurons];
		for (int n = 0; n < numNeurons; n++) {
			avgCoefficient[n] = new float[numExc];
		}

		// 3d array to store mean FR for each exc neuron (indi-nExc-BinsAll)
		float ***smoothedFRsAVG = new float**[indiNum];
		for (int i = 0; i < indiNum; i ++) {
			smoothedFRsAVG[i] = new float*[numExc];
			for (int n = 0; n < numExc; n++) {
				smoothedFRsAVG[i][n] = new float[BINS_ALL];
				for (int b = 0; b < BINS_ALL; b++) {
					smoothedFRsAVG[i][n][b] = 0.0;
				}
			}
		}

		// 3d array to store mean FR for even trials (indi-nExc-BinsAll)
		float ***smoothedFRsAVG_EVEN = new float**[indiNum];
		for (int i = 0; i < indiNum; i ++) {
			smoothedFRsAVG_EVEN[i] = new float*[numExc];
			for (int n = 0; n < numExc; n++) {
				smoothedFRsAVG_EVEN[i][n] = new float[BINS_ALL];
				for (int b = 0; b < BINS_ALL; b++) {
					smoothedFRsAVG_EVEN[i][n][b] = 0.0;
				}
			}
		}
		// 3d array to store mean FR for odd trials (indi-nExc-BinsAll)
		float ***smoothedFRsAVG_ODD = new float**[indiNum];
		for (int i = 0; i < indiNum; i ++) {
			smoothedFRsAVG_ODD[i] = new float*[numExc];
			for (int n = 0; n < numExc; n++) {
				smoothedFRsAVG_ODD[i][n] = new float[BINS_ALL];
				for (int b = 0; b < BINS_ALL; b++) {
					smoothedFRsAVG_ODD[i][n][b] = 0.0;
				}
			}
		}

		float **inpSmthFRsAVG = new float*[INP_N];
		for (int n = 0; n < INP_N; n++) {
			inpSmthFRsAVG[n] = new float[BINS_ALL];
			for (int b = 0; b < BINS_ALL; b++) {
				inpSmthFRsAVG[n][b] = 0.0;
			}
		}

		ofstream lesionFitScores;
		CARLsim* network;

		if (loadFlag == 0) {
			network = new CARLsim("direct_driver", simMode, verbosity, seed);
		}
		else {
			/* Seeds for each evolved network, necessary for loading. Uncomment the seed corresponding to the network being loaded. 'NETWORK X' refers to the generation number associated with the network being loaded (for example, NETWORK 36 indicates indirect_network36.dat). The corresponding evolutionary run number is listed afterward (so, for network 36, the evolutionary run is 10). Lastly, the best fit individual is listed (in the case of network 36, i = 7). */
				
			//time_t seed = 1446878707; // NETWORK 30 - 1 - 9
			//time_t seed = 1456539732; // NETWORK 32 - 2 - 6
			//time_t seed = 1457264465; // NETWORK 49 - 3 - 1
			//time_t seed = 1459526736; // NETWORK 9 - 4 - 5
			//time_t seed = 1460532400; // NETWORK 32 - 5 - 6
			//time_t seed = 1461052506; // NETWORK 43 - 6 - 4
			//time_t seed = 1461413861; // NETWORK 48 - 7 - 1
			//time_t seed = 1461512243; // NETWORK 7 - 8 - 8
			//time_t seed = 1462071562; // NETWORK 31 - 9 - 7
			// seed = 1462301880; // NETWORK 36 - 10 - 7
			network = new CARLsim("loadSim",  simMode, verbosity, seed);
		}

		float inpToExcMaxWght, inpToInhMaxWght, excToExcMaxWght, inhToExcMaxWght;
		float inpAVToExcMaxWght, inpLVToExcMaxWght, inpHDToExcMaxWght, inpPosToExcMaxWght;
		float inAVpToInhMaxWght, inpLVToInhMaxWght, inpHDToInhMaxWght, inpPosToInhMaxWght;
		
		float inpAVConnProb, inpLVConnProb, inpHDConnProb, inpPosConnProb;

		float excBaseFR, inhBaseFR;
		float ALPHA_LTP_EE, ALPHA_LTD_EE, TAU_LTP_EE, TAU_LTD_EE;
		float ALPHA_LTP_EI, ALPHA_LTD_EI, TAU_LTP_EI, TAU_LTD_EI;
		float ALPHA_LTP_IE, ALPHA_LTD_IE, TAU_LTP_IE, TAU_LTD_IE;

		int inputs;
		int gInAngVel;
		int gInLinVel;
		int gInHD;
		int gInPos;
		float initWtScale = 0.5;

		short int inpToExc[indiNum];
		short int inpToInh[indiNum];
		short int inAngVelToExc[indiNum];
		short int inLinVelToExc[indiNum];
		short int inHDToExc[indiNum];
		short int inPosToExc[indiNum];

		short int inAngVelToInh[indiNum];
		short int inLinVelToInh[indiNum];
		short int inHDToInh[indiNum];
		short int inPosToInh[indiNum];

		short int inhToExc[indiNum];
		short int excToExc[indiNum];
		if (inputOneGroup) {
			inputs = network->createSpikeGeneratorGroup("inputs", INP_N, EXCITATORY_NEURON);
		} else {
			gInAngVel = network->createSpikeGeneratorGroup("inAngVel", INP_ANG, EXCITATORY_POISSON);
			gInLinVel = network->createSpikeGeneratorGroup("inLinVel", INP_LV, EXCITATORY_POISSON);
			gInHD = network->createSpikeGeneratorGroup("inHD", INP_HD, EXCITATORY_POISSON);
			gInPos = network->createSpikeGeneratorGroup("inPos", INP_POS, EXCITATORY_POISSON);
		}

		for (int i = 0; i < indiNum; i++){

			ALPHA_LTP_EE = parameters.getParameter(i,0);
			ALPHA_LTD_EE = parameters.getParameter(i,3);
			TAU_LTP_EE = parameters.getParameter(i,6);
			TAU_LTD_EE = parameters.getParameter(i,7);

			ALPHA_LTP_EI = parameters.getParameter(i,1);
			ALPHA_LTD_EI = parameters.getParameter(i,4);
			TAU_LTP_EI = parameters.getParameter(i,8);
			TAU_LTD_EI = parameters.getParameter(i,9);

			ALPHA_LTP_IE = parameters.getParameter(i,5);
			ALPHA_LTD_IE = parameters.getParameter(i,2);
			TAU_LTP_IE = parameters.getParameter(i,10);
			TAU_LTD_IE = parameters.getParameter(i,11);

			excBaseFR = parameters.getParameter(i,19);
			if (inputOneGroup) {
				inhBaseFR = parameters.getParameter(i,12);

				inpToExcMaxWght = parameters.getParameter(i,13);
				inpToInhMaxWght = parameters.getParameter(i,14);
				excToExcMaxWght = parameters.getParameter(i,15);
				inhToExcMaxWght = parameters.getParameter(i,16);

				HOMEO_AVERAGE_TIME_SCALE_EXC = parameters.getParameter(i,17);
				HOMEO_AVERAGE_TIME_SCALE_INH = parameters.getParameter(i,18);
			} else {
				inhBaseFR = parameters.getParameter(i,12);
				inpAVToExcMaxWght = parameters.getParameter(i,13);
				inpLVToExcMaxWght = parameters.getParameter(i,13);
				inpHDToExcMaxWght = parameters.getParameter(i,13);
				inpPosToExcMaxWght = parameters.getParameter(i,13);

				inAVpToInhMaxWght = parameters.getParameter(i,14);
				inpLVToInhMaxWght = parameters.getParameter(i,14);
				inpHDToInhMaxWght = parameters.getParameter(i,14);
				inpPosToInhMaxWght = parameters.getParameter(i,14);

				excToExcMaxWght = parameters.getParameter(i,15);
				inhToExcMaxWght = parameters.getParameter(i,16);

				HOMEO_AVERAGE_TIME_SCALE_EXC = parameters.getParameter(i,17);
				HOMEO_AVERAGE_TIME_SCALE_INH = parameters.getParameter(i,18);

				inpAVConnProb = parameters.getParameter(i,20);
				inpLVConnProb = parameters.getParameter(i,21);
				inpHDConnProb = parameters.getParameter(i,22);
				inpPosConnProb  = parameters.getParameter(i,23);
			}

			exc[i] = network->createGroup("excitatory", numExc, EXCITATORY_NEURON);
			network->setNeuronParameters(exc[i], REG_IZH[0], REG_IZH[1], REG_IZH[2], REG_IZH[3]);
			inh[i] = network->createGroup("inhibitory", numInh, INHIBITORY_NEURON);
			network->setNeuronParameters(inh[i], FAST_IZH[0], FAST_IZH[1], FAST_IZH[2], FAST_IZH[3]);
			network->setConductances(true,COND_tAMPA,COND_tNMDA,COND_tGABAa,COND_tGABAb);

			if (inputOneGroup) {
				inpToExc[i] = network->connect(inputs, exc[i], "random", 
					RangeWeight(0.0f, inpToExcMaxWght*initWtScale, inpToExcMaxWght), 0.1f, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
				inpToInh[i] = network->connect(inputs, inh[i], "random", 
					RangeWeight(0.0f, inpToInhMaxWght*initWtScale, inpToInhMaxWght), 0.1f, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
			} else {
				// connections from input groups to exc group
				inAngVelToExc[i] = network->connect(gInAngVel, exc[i], "random", 
					RangeWeight(0.0f,inpAVToExcMaxWght*initWtScale,inpAVToExcMaxWght), inpAVConnProb, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
				inLinVelToExc[i] = network->connect(gInLinVel, exc[i], "random", 
					RangeWeight(0.0f,inpLVToExcMaxWght*initWtScale,inpLVToExcMaxWght), inpLVConnProb, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
				inHDToExc[i] = network->connect(gInHD, exc[i], "random", 
					RangeWeight(0.0f,inpHDToExcMaxWght*initWtScale,inpHDToExcMaxWght), inpHDConnProb, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
				inPosToExc[i] = network->connect(gInPos, exc[i], "random", 
					RangeWeight(0.0f,inpPosToExcMaxWght*initWtScale,inpPosToExcMaxWght), inpPosConnProb, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);

				// connections from input groups to inh group
				inAngVelToInh[i] = network->connect(gInAngVel, inh[i], "random", 
					RangeWeight(0.0f,inAVpToInhMaxWght*initWtScale,inAVpToInhMaxWght), inpAVConnProb, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
				inLinVelToInh[i] = network->connect(gInLinVel, inh[i], "random", 
					RangeWeight(0.0f,inpLVToInhMaxWght*initWtScale,inpLVToInhMaxWght), inpLVConnProb, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
				inHDToInh[i] = network->connect(gInHD, inh[i], "random", 
					RangeWeight(0.0f,inpHDToInhMaxWght*initWtScale,inpHDToInhMaxWght), inpHDConnProb, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
				inPosToInh[i] = network->connect(gInPos, inh[i], "random", 
					RangeWeight(0.0f,inpPosToInhMaxWght*initWtScale,inpPosToInhMaxWght), inpPosConnProb, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
			}
			excToExc[i] = network->connect(exc[i], exc[i], "random", 
				RangeWeight(0.0f, excToExcMaxWght*initWtScale, excToExcMaxWght), 0.1f, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
			inhToExc[i] = network->connect(inh[i], exc[i], "random", 
				RangeWeight(0.0f, inhToExcMaxWght*initWtScale, inhToExcMaxWght), 0.1f, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);

			network->setESTDP(exc[i], true, STANDARD, ExpCurve(ALPHA_LTP_EE, TAU_LTP_EE, ALPHA_LTD_EE, TAU_LTD_EE));
			network->setISTDP(exc[i], true, STANDARD, ExpCurve(ALPHA_LTP_IE, TAU_LTP_IE, TAU_LTP_IE, TAU_LTD_IE));
			network->setESTDP(inh[i], true, STANDARD, ExpCurve(ALPHA_LTP_EI, TAU_LTP_EI, ALPHA_LTD_EI, TAU_LTD_EI));
			network->setHomeostasis(exc[i], true, HOMEO_FACTOR, HOMEO_AVERAGE_TIME_SCALE_EXC);
			network->setHomeostasis(inh[i], true, HOMEO_FACTOR, HOMEO_AVERAGE_TIME_SCALE_INH);

			network->setHomeoBaseFiringRate(exc[i], excBaseFR, 0);
			network->setHomeoBaseFiringRate(inh[i], inhBaseFR, 0);
		}			

		time_t trialTimer_start, trialTimer_end;

		FILE* simFid;
		if (loadFlag == 1) {
			simFid = NULL;
			simFid = fopen(sim_name.c_str(), "rb");
			cout << "Loading network" << endl;
			network->loadSimulation(simFid);
		}

		network->setupNetwork();

		if (loadFlag == 1) {
			fclose(simFid);
		}

		string ch2;
		stringstream strs2;
		string simdataDir;
		// string net_folder = "./";
		string driver_str = "_Driver/";
		string recstr = "recordings.csv";
		string repstr = "repetitions.csv";

		if (loadFlag == 0) {
			if (writeSimData) {
				for (int i = 0; i < indiNum; i++) {
					strs2 << i;
					ch2 = strs2.str();
					strs2.str("");
					simdataDir = (net_folder+ch2+driver_str);

					if (mkdir(simdataDir.c_str(), 0777) == -1) {
						 cerr << "Error :  " << strerror(errno) << endl; 
						 return;
					}
				}
			}
		}

		string conn = "conn_";
		string connInpToExc = "inp_exc";
		string connInpToInh = "inp_inh";
		string connExcToExc = "exc_exc";
		string connInhToExc = "inh_exc";

		string conn1 = "inAngVel_exc";
		string conn2 = "inHD_exc";
		string conn3 = "inPos_exc";
		string conn4 = "inLinVel_exc";

		string ind = "results/spk_";
		string ind1 = "excitatory_";
		string ind2 = "inhibitory_";
		string ind3 = "excInps_";
		string e = ".dat";

		if (loadFlag == 0) {
			for (int i = 0; i < indiNum; i++) {
				stringstream strs1;
				strs1 << i;
				string ch1 = strs1.str();

				if (writeSimData) {
					if (inputOneGroup) {
						string allConnInpToExc = (net_folder+ch1+driver_str+conn+connInpToExc+e);
						inpToExcCM[i] = network->setConnectionMonitor(inputs, exc[i], allConnInpToExc.c_str());
						inpToExcCM[i]->setUpdateTimeIntervalSec(-1);

						string allConnInpToInh = (net_folder+ch1+driver_str+conn+connInpToInh+e);
						inpToInhCM[i] = network->setConnectionMonitor(inputs, inh[i], allConnInpToInh.c_str());
						inpToInhCM[i]->setUpdateTimeIntervalSec(-1);
					} else {
						string allConnInpAVToExc = (net_folder+ch1+driver_str+conn+conn1+e);
						inpAVToExcCM[i] = network->setConnectionMonitor(gInAngVel, exc[i], allConnInpAVToExc.c_str());
						inpAVToExcCM[i]->setUpdateTimeIntervalSec(-1);

						string allConnInpHDToExc = (net_folder+ch1+driver_str+conn+conn2+e);
						inpHDToExcCM[i] = network->setConnectionMonitor(gInHD, exc[i], allConnInpHDToExc.c_str());
						inpHDToExcCM[i]->setUpdateTimeIntervalSec(-1);

						string allConnInpPosToExc = (net_folder+ch1+driver_str+conn+conn3+e);
						inpPosToExcCM[i] = network->setConnectionMonitor(gInPos, exc[i], allConnInpPosToExc.c_str());
						inpPosToExcCM[i]->setUpdateTimeIntervalSec(-1);

						string allConnInpLVToExc = (net_folder+ch1+driver_str+conn+conn4+e);
						inpLVToExcCM[i] = network->setConnectionMonitor(gInLinVel, exc[i], allConnInpLVToExc.c_str());
						inpLVToExcCM[i]->setUpdateTimeIntervalSec(-1);
					}

					string allConnExcToExc = (net_folder+ch1+driver_str+conn+connExcToExc+e);
					excToExcCM[i] = network->setConnectionMonitor(exc[i], exc[i], allConnExcToExc.c_str());
					excToExcCM[i]->setUpdateTimeIntervalSec(-1);

					string allConnInhToExc = (net_folder+ch1+driver_str+conn+connInhToExc+e);
					inhToExcCM[i] = network->setConnectionMonitor(inh[i], exc[i], allConnInhToExc.c_str());
					inhToExcCM[i]->setUpdateTimeIntervalSec(-1);

				} else {
					if (inputOneGroup) {
						inpToExcCM[i] = network->setConnectionMonitor(inputs, exc[i], "NULL");
						inpToExcCM[i]->setUpdateTimeIntervalSec(-1);

						inpToInhCM[i] = network->setConnectionMonitor(inputs, inh[i], "NULL");
						inpToInhCM[i]->setUpdateTimeIntervalSec(-1);
					} else {
						string allConnInpAVToExc = (conn+conn1+e);
						inpAVToExcCM[i] = network->setConnectionMonitor(gInAngVel, exc[i], "NULL");
						inpAVToExcCM[i]->setUpdateTimeIntervalSec(-1);

						string allConnInpHDToExc = (conn+conn2+e);
						inpHDToExcCM[i] = network->setConnectionMonitor(gInHD, exc[i], "NULL");
						inpHDToExcCM[i]->setUpdateTimeIntervalSec(-1);

						string allConnInpPosToExc = (conn+conn3+e);
						inpPosToExcCM[i] = network->setConnectionMonitor(gInPos, exc[i], "NULL");
						inpPosToExcCM[i]->setUpdateTimeIntervalSec(-1);

						string allConnInpLVToExc = (conn+conn4+e);
						inpLVToExcCM[i] = network->setConnectionMonitor(gInLinVel, exc[i], "NULL");
						inpLVToExcCM[i]->setUpdateTimeIntervalSec(-1);
					}
					excToExcCM[i] = network->setConnectionMonitor(exc[i], exc[i], "NULL");
					excToExcCM[i]->setUpdateTimeIntervalSec(-1);

					inhToExcCM[i] = network->setConnectionMonitor(inh[i], exc[i], "NULL");
					inhToExcCM[i]->setUpdateTimeIntervalSec(-1);
				}
			}
		}
		for (int i = 0; i < indiNum; i++) {
			// string allExc = (ind+ind1+e);
			excMonitor[i] = network->setSpikeMonitor(exc[i], "NULL");
			// string allInh = (ind+ind2+e);
			inhMonitor[i] = network->setSpikeMonitor(inh[i], "NULL");
		}

		vector<vector<float> > posXY;
		vector<float> angVel;
		vector<float> linVel;
		vector<float> headDir;
		vector<float> occ;

		vector<vector<float> > coordsFR;
		vector<vector<float> > angVelFR;
		vector<vector<float> > linVelFR;
		vector<vector<float> > headDirFR;
		int numBins;
		int startBinInd;

		int numTeTrials;

		int mod;
		int trialNumTe;

		// string rec_name = (net_folder+ch2+driver_str+recstr);
		// string rep_name = (net_folder+ch2+driver_str+repstr);
		string rec_name = (net_folder+recstr);
		string rep_name = (net_folder+repstr);

		ofstream writeRecordings;
		ofstream writeRepetitions;

		int recordings[totalTeTrials] = {};
		int repetitions[totalTeTrials][numPaths][repTrialsTe] = {};

		if (forceMatch == 0) {
			writeRecordings.open(rec_name.c_str());
			writeRepetitions.open(rep_name.c_str());
		}
		else {
			ifstream getRecordings;
			getRecordings.open(rec_name.c_str());

			for (int total = 0; total < totalTeTrials; total++) {
				getRecordings >> recordings[total];
				getRecordings.get();
			}
			getRecordings.close();

			ifstream getRepetitions;
			getRepetitions.open(rep_name.c_str());
			for (int total = 0; total < totalTeTrials; total++) {
				for (int r = 0; r < numPaths; r++) {
					for (int z = 0; z < repTrialsTe; z++) {
						getRepetitions >> repetitions[total][r][z];
						// cout << "TRIAL NUM IS: " << repetitions[total][r][z] << endl;
						getRepetitions.get();
					}
				}
			}
			getRepetitions.close();
		}

		/********************************************************************************
									Training
		********************************************************************************/
		if (loadFlag == 0) {
			testing = 0;
			for (int total = 0; total < totalTrTrials; total++) {
				recNum = (rand() % numRecordings);
				time(&trialTimer_start);
				for (int reps = 0; reps < numPaths; reps++) {
					if (reps < NUM_IN) {
						numBins = BINS_IN;
						startBinInd = reps * BINS_IN;
					} else {
						numBins = BINS_OUT;
						startBinInd = BINS_IN * NUM_IN + BINS_OUT * (reps-NUM_IN);
					}
					for (int z = 0; z < repTrialsTe; z++) {
						numTeTrials = trialCountTrain[recNum][reps];							

						if (numTeTrials > 0) {
							trialNumTe = (rand() % numTeTrials);
						}
						else {
							trialNumTe = 0;
						}

						posXY = getXYData(dataDir, recNum, reps, testing, trialNumTe, numBins);						
						coordsFR = gauss2dXY(posXY[0], posXY[1], xCoordPref, yCoordPref, INP_POS_X, INP_POS_Y, numBins, sigma_x, r_max);						
						// coordsFR = gauss2d(posXY[0], posXY[1], xCoordPref, yCoordPref, INP_POS, numBins, sigma_x, r_max);	
						angVel = getAngVel(dataDir, recNum, reps, testing, trialNumTe, numBins);
						angVelFR = gauss1d(angVel, angVelPref, INP_ANG, numBins, sigmaAngVel, r_max_ang_lin);
						linVel = getLinVel(dataDir, recNum, reps, testing, trialNumTe, numBins);
						linVelFR = gauss1d(linVel, linVelPref, INP_LV, numBins, sigmaLinVel, r_max_ang_lin);
						headDir = getHeadDir(dataDir, recNum, reps, testing, trialNumTe, numBins);
						headDirFR = cosine(headDir, headDirPref, INP_HD, numBins, r_max, r_min);
						// occ = getOcc(dataDir, recNum, reps, testing, trialNumTe, numBins);
						time(&timerTrial_start);

						PoissonRate inpIn(INP_N);
						PoissonRate inpAV(INP_ANG);
						PoissonRate inpLV(INP_LV);
						PoissonRate inpHD(INP_HD);
						PoissonRate inpPos(INP_POS);

						// BEGIN BINS LOOP
						//time(&timerBin_start);
						for (int b = 0; b < numBins; b++){
							if (inputOneGroup) {
								count = 0;
								for (int j = 0; j < INP_ANG; j++) {
									inpIn.setRate(count, angVelFR[j][b]);
									count ++;
								}
								for (int j = 0; j < INP_HD; j++) {
									inpIn.setRate(count, headDirFR[j][b]);
									count ++;
								}
								for (int j = 0; j < INP_POS; j++) {
									inpIn.setRate(count, coordsFR[j][b]);
									count ++;
								}
								for (int j = 0; j < INP_LV; j ++){
									inpIn.setRate(count, linVelFR[j][b]);
									count ++;
								}

								network->setSpikeRate(inputs, &inpIn);
							} else {
								for (int j = 0; j < INP_ANG; j++) {
									inpAV.setRate(j, angVelFR[j][b]);
								}
								for (int j = 0; j < INP_HD; j++) {
									inpHD.setRate(j, headDirFR[j][b]);
								}
								for (int j = 0; j < INP_POS; j++) {
									inpPos.setRate(j, coordsFR[j][b]);
								}
								for (int j = 0; j < INP_LV; j ++){
									inpLV.setRate(j, linVelFR[j][b]);
								}

								network->setSpikeRate(gInAngVel, &inpAV);
								network->setSpikeRate(gInLinVel, &inpLV);
								network->setSpikeRate(gInHD, &inpHD);
								network->setSpikeRate(gInPos, &inpPos);
							} // end if
							network->runNetwork(runTimeSec,runTimeMs);
						} // end bins
					} // end rep trials
				} // end path
			} // end rec
			if (writeSimData) {
				network->saveSimulation(sim_name, true);
			}
		} 

		network->startTesting();

		// vector< vector< float > > weights;

		if (loadFlag == 0) {
			for (int i = 0; i < indiNum; i++) {
				if (inputOneGroup) {
					inpToExcCM[i]->takeSnapshot();
					inpToInhCM[i]->takeSnapshot();
				} else {
					inpAVToExcCM[i]->takeSnapshot();
					inpHDToExcCM[i]->takeSnapshot();
					inpPosToExcCM[i]->takeSnapshot();
					inpLVToExcCM[i]->takeSnapshot();
				}
				excToExcCM[i]->takeSnapshot();
				inhToExcCM[i]->takeSnapshot();
			}
		}

		for (int i = 0; i < indiNum; i++) {
			if (lesionAV == 1) {
				for (int preId = 0; preId < INP_ANG; preId++) {
					if (inputOneGroup) {
						float w = 0.0;
						// lesion inp to exc
						for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
							network->setWeight(inpToExc[i], preId, postId, w); // between 0 and 1
						}
						// lesion inp to inh
						for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
							network->setWeight(inpToInh[i], preId, postId, w); // between 0 and 1
						}
					} else {
						float w = 0.0;
						// lesion inp to exc
						for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
							network->setWeight(inAngVelToExc[i], preId, postId, w); // between 0 and 1
						}
						// lesion inp to inh
						for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
							network->setWeight(inAngVelToInh[i], preId, postId, w); // between 0 and 1
						}
					}
				}			
			}
			if (lesionHD == 1) {
				for (int preId = INP_ANG; preId < INP_POS_START; preId++) {
					if (inputOneGroup) {
						float w = 0.0;
						for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
							network->setWeight(inpToExc[i], preId, postId, w); // between 0 and 1
						}
						for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
							network->setWeight(inpToInh[i], preId, postId, w); // between 0 and 1
						}
					} else {
						float w = 0.0;
						for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
							network->setWeight(inHDToExc[i], preId-INP_ANG, postId, w); // between 0 and 1
						}
						for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
							network->setWeight(inHDToInh[i], preId-INP_ANG, postId, w); // between 0 and 1
						}
					}
				}			
			}

			if (lesionPos == 1) {
				for (int preId = INP_POS_START; preId < INP_POS_END; preId++) {
					if (inputOneGroup) {
						float w = 0.0;
						for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
							network->setWeight(inpToExc[i], preId, postId, w); // between 0 and 1
						}
						for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
							network->setWeight(inpToInh[i], preId, postId, w); // between 0 and 1
						}
					} else {
						float w = 0.0;
						for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
							network->setWeight(inPosToExc[i], preId-INP_POS_START, postId, w); // between 0 and 1
						}
						for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
							network->setWeight(inPosToInh[i], preId-INP_POS_START, postId, w); // between 0 and 1
						}
					}
				}			
			}

			if (lesionLV == 1) {
				for (int preId = INP_POS_END; preId < INP_N; preId++) {
					if (inputOneGroup) {
						float w = 0.0;
						for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
							network->setWeight(inpToExc[i], preId, postId, w); // between 0 and 1
						}
						for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
							network->setWeight(inpToInh[i], preId, postId, w); // between 0 and 1
						}
					} else {
						float w = 0.0;
						for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
							network->setWeight(inLinVelToExc[i], preId-INP_POS_END, postId, w); // between 0 and 1
						}
						for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
							network->setWeight(inLinVelToInh[i], preId-INP_POS_END, postId, w); // between 0 and 1
						}
					}
				}			
			}	

			if (lesionExc == 1) {
				for (int preId = 0; preId < network->getGroupNumNeurons(exc[i]); preId++) {
					for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
						if (preId != postId) {
							float w = 0.0;
							network->setWeight(excToExc[i], preId, postId, w); // between 0 and 1
						}
					}
				}		
			}

			if (lesionInh == 1) {
				for (int preId = 0; preId < network->getGroupNumNeurons(inh[i]); preId++) {
					for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
						float w = 0.0;
						network->setWeight(inhToExc[i], preId, postId, w); // between 0 and 1
					}
				}		
			}
		}
		/********************************************************************************
									Testing
		********************************************************************************/		
		// Loop over every recording
		testing = 1;

		for (int total = 0; total < totalTeTrials; total++) {
			if (forceMatch == 0) {
				recNum = (rand() % numRecordings);
				writeRecordings << recNum << ",";
			}
			else {
				recNum = recordings[total];
				//cout << recNum << endl;
			}

			time(&trialTimer_start);
			for (int reps = 0; reps < numPaths; reps++) {
				if (reps < NUM_IN) {
					numBins = BINS_IN;
					startBinInd = reps * BINS_IN;
				} else {
					numBins = BINS_OUT;
					startBinInd = BINS_IN * NUM_IN + BINS_OUT * (reps-NUM_IN);
				}
				for (int z = 0; z < repTrialsTe; z++) {
					if (forceMatch == 0) {
						numTeTrials = trialCountTest[recNum][reps];

						if (numTeTrials > 0) {
							trialNumTe = (rand() % numTeTrials);
						}
						else {
							trialNumTe = 0;
						}
						writeRepetitions << trialNumTe << ",";					
					}
					else {
						trialNumTe = repetitions[total][reps][z];
						//cout << trialNumTe << endl;
					}

					posXY = getXYData(dataDir, recNum, reps, testing, trialNumTe, numBins);						
					coordsFR = gauss2dXY(posXY[0], posXY[1], xCoordPref, yCoordPref, INP_POS_X, INP_POS_Y, numBins, sigma_x, r_max);						
					// coordsFR = gauss2d(posXY[0], posXY[1], xCoordPref, yCoordPref, INP_POS, numBins, sigma_x, r_max);	
					angVel = getAngVel(dataDir, recNum, reps, testing, trialNumTe, numBins);
					angVelFR = gauss1d(angVel, angVelPref, INP_ANG, numBins, sigmaAngVel, r_max_ang_lin);
					linVel = getLinVel(dataDir, recNum, reps, testing, trialNumTe, numBins);
					linVelFR = gauss1d(linVel, linVelPref, INP_LV, numBins, sigmaLinVel, r_max_ang_lin);
					headDir = getHeadDir(dataDir, recNum, reps, testing, trialNumTe, numBins);
					headDirFR = cosine(headDir, headDirPref, INP_HD, numBins, r_max, r_min);
					// occ = getOcc(dataDir, recNum, reps, testing, trialNumTe, numBins);					

					time(&timerTrial_start);
					// BEGIN BINS LOOP
					//time(&timerBin_start);
					vector<vector<vector<float> > > excFRs(indiNum, vector<vector<float> >(numBins, vector<float>(numExc, 0)));

					for (int b = 0; b < numBins; b++){
						inpsDataAll[total][z][0][b+startBinInd] = angVel[b];
						inpsDataAll[total][z][1][b+startBinInd] = headDir[b];
						inpsDataAll[total][z][2][b+startBinInd] = posXY[0][b];
						inpsDataAll[total][z][3][b+startBinInd] = posXY[1][b];
						inpsDataAll[total][z][4][b+startBinInd] = linVel[b];
						// inpsSmthFRs[total][z][5][b+startBinInd] = occ[b];
						PoissonRate inpIn(INP_N);
						PoissonRate inpAV(INP_ANG);
						PoissonRate inpLV(INP_LV);
						PoissonRate inpHD(INP_HD);
						PoissonRate inpPos(INP_POS);

						if (inputOneGroup) {
							count = 0;
							for (int j = 0; j < INP_ANG; j++) {
								inpIn.setRate(count, angVelFR[j][b]);
								inpsSmthFRs[total][z][count][b+startBinInd] = angVelFR[j][b];
								count ++;
							}
							for (int j = 0; j < INP_HD; j++) {
								inpIn.setRate(count, headDirFR[j][b]);
								inpsSmthFRs[total][z][count][b+startBinInd] = headDirFR[j][b];
								count ++;
							}
							for (int j = 0; j < INP_POS; j++) {
								inpIn.setRate(count, coordsFR[j][b]);
								inpsSmthFRs[total][z][count][b+startBinInd] = coordsFR[j][b];
								count ++;
							}
							for (int j = 0; j < INP_LV; j ++){
								inpIn.setRate(count, linVelFR[j][b]);
								inpsSmthFRs[total][z][count][b+startBinInd] = linVelFR[j][b];
								count ++;
							}
							network->setSpikeRate(inputs, &inpIn);
						} else {
							count = 0;
							for (int j = 0; j < INP_ANG; j++) {
								inpsSmthFRs[total][z][count][b+startBinInd] = angVelFR[j][b];
								inpAV.setRate(j, angVelFR[j][b]);
								count ++;

							}
							for (int j = 0; j < INP_HD; j++) {
								inpsSmthFRs[total][z][count][b+startBinInd] = headDirFR[j][b];
								inpHD.setRate(j, headDirFR[j][b]);
								count ++;

							}
							for (int j = 0; j < INP_POS; j++) {
								inpsSmthFRs[total][z][count][b+startBinInd] = coordsFR[j][b];
								inpPos.setRate(j, coordsFR[j][b]);
								count ++;

							}
							for (int j = 0; j < INP_LV; j ++){
								inpsSmthFRs[total][z][count][b+startBinInd] = linVelFR[j][b];
								inpLV.setRate(j, linVelFR[j][b]);
								count ++;

							}

							network->setSpikeRate(gInAngVel, &inpAV);
							network->setSpikeRate(gInLinVel, &inpLV);
							network->setSpikeRate(gInHD, &inpHD);
							network->setSpikeRate(gInPos, &inpPos);
						}

						for (int i = 0; i < indiNum; i++) {
							excMonitor[i]->startRecording();
							inhMonitor[i]->startRecording();
						}

						// runTimeMs = occ[b] / 60 * 1000;
						// runTimeSec = int(floor(runTimeMs / 1000));
						// runTimeMs = runTimeMs % 1000;

						network->runNetwork(runTimeSec,runTimeMs);

						for (int i = 0; i < indiNum; i++) {
							excMonitor[i]->stopRecording();
							inhMonitor[i]->stopRecording();

							long int recTime = excMonitor[i]->getRecordingTotalTime();
							// For this config ID, store the spikes for this BIN and TRIAL.
							for(int index = 0; index < numExc; index++){
								// Grab the spike vectors for this CONFIG ID and BIN (on current trial).
								exc2D = 0;
								//exc_temp[index] = 0.0;
								exc2D = excMonitor[i]->getNeuronNumSpikes(index);
								// if (recTime > 0) {
									excFRs[i][b][index] = (exc2D * 1000.0) / recTime;
								// } else {
									// excFRs[i][b][index] = 0.0;
								// }
							}
							// for(int index = 0; index < numInh; index++){
							// 	// Grab the spike vectors for this CONFIG ID and BIN (on current trial).
							// 	inh2D = 0;
							// 	//exc_temp[index] = 0.0;
							// 	inh2D = inhMonitor[i]->getNeuronNumSpikes(index);
							// 	// if (i == 0 && total == 0 && trial == 0) {
							// 	// 	cout << " Bin " << b << ", Neuron " << index << " spikes: " << exc2D << endl;
							// 	// }
							// 	inhFRs[b][index] = (inh2D * 1000.0) / recTime;
							// }
							excMonitor[i]->clear();
							inhMonitor[i]->clear();

							// If we are on the last bin, then we can evaluate the fitness for a given neuron.
							// If we are on the last bin, we can smooth and store the data for each neuron.
							if (b == (numBins-1)) { 
								for (int i = 0; i < indiNum; i++) {
									for (int n = 0; n < numExc; n++) {

										// Create a temporary vector for smoothing.
										float* temp_vec_in = new float[numBins];
										float* temp_vec_out = new float[numBins];
										for (int l = 0; l < numBins; l++) {
											// if (n < (NUM_N*0.8)) {
											temp_vec_in[l] = excFRs[i][l][n];
											// }
											// else {
											// 	temp_vec_in[l] = inhFRs[l][n];
											// }
											temp_vec_out[l] = 0.0;
										}

										// Smooth the data.
										custom(temp_vec_in, temp_vec_out, numBins);

										count = 0;
										// Store the smoothed data.
										for (int l = 0; l < numBins; l++) {
											smoothedFRs[i][total][z][n][l+startBinInd] = temp_vec_out[l];
										}
										delete[] temp_vec_in;
										delete[] temp_vec_out;
									}
								} // end indi
							} // Only if last bin
						} // End bin loop
					} // end pos/traj loop
				} // end trial loop
				if (forceMatch == 0) {
					writeRepetitions << endl;
				}
				time(&trialTimer_end);
				// cout << "Finished testing trial " << total+1 << "/" << totalTeTrials << " with time " 
					// << difftime(trialTimer_end, trialTimer_start) << endl;
				time(&trialTimer_start);
			} 

		} // END TESTING
		if (forceMatch == 0) {
			writeRecordings.close();
			writeRepetitions.close();
		}

		network->stopTesting();

		// BEGIN FITNESS CALC
		// float bestCorr = -1.0;
		vector<vector<float> > bestCorr(indiNum, vector<float>(numNeurons, -1.0));
		int temp_index = -2;
		float coefficient = 0.0;

		float avgSmooth[indiNum] = {};
		float inpAvgSmooth = 0.0;
		float avgSmoothR[indiNum] = {};
		float inpAvgSmthR = 0.0;
		float avgSmoothOdd[indiNum] = {};
		float avgSmoothROdd[indiNum] = {};
		float avgSmoothEven[indiNum] = {};
		float avgSmoothREven[indiNum] = {};

		// BEGIN FITNESS CALC

		// Calculate averaged firing rates (across trials).
		// Divide into even and odd pools for reconstruction analyses.
		int numOdd = 0;
		int numEven = 0;

		for (int n = 0; n < numExc; n++) {
			for (int b = 0; b < BINS_ALL; b++) {
				for (int r = 0; r < totalTeTrials; r++) {
					for (int z = 0; z < repTrialsTe; z++) {
						if (n < INP_N) {
							inpAvgSmooth += inpsSmthFRs[r][z][n][b];
						}
						for (int i = 0; i < indiNum; i ++) {
							avgSmooth[i] += smoothedFRs[i][r][z][n][b];
						}
						if ( (z+1) % 2 == 0) {
							for (int i = 0; i < indiNum; i ++) {
								avgSmoothEven[i] += smoothedFRs[i][r][z][n][b];
							}
							numEven ++;
						}
						else {
							for (int i = 0; i < indiNum; i ++) {
								avgSmoothOdd[i] += smoothedFRs[i][r][z][n][b];
							}
							numOdd ++;
						}
					}
					if (n < INP_N) {
						inpAvgSmthR += (inpAvgSmooth / repTrialsTe);
					}
					for (int i = 0; i < indiNum; i ++) {
						avgSmoothR[i] += (avgSmooth[i] / repTrialsTe);
						avgSmoothREven[i] += (avgSmoothEven[i] / numEven);
						avgSmoothROdd[i] += (avgSmoothOdd[i] / numOdd);
						avgSmooth[i] = 0.0;
						avgSmoothOdd[i] = 0.0;
						avgSmoothEven[i] = 0.0;
					}
					inpAvgSmooth = 0.0;
					numEven = 0;
					numOdd = 0;
				}
				if (n < INP_N) {
					inpSmthFRsAVG[n][b] = (inpAvgSmthR / totalTeTrials);
				}
				for (int i = 0; i < indiNum; i ++) {
					smoothedFRsAVG[i][n][b] = (avgSmoothR[i] / totalTeTrials);
					smoothedFRsAVG_EVEN[i][n][b] = (avgSmoothREven[i] / totalTeTrials);
					smoothedFRsAVG_ODD[i][n][b] = (avgSmoothROdd[i] / totalTeTrials);
					avgSmoothR[i] = 0.0;
					avgSmoothREven[i] = 0.0;
					avgSmoothROdd[i] = 0.0;
				}
				inpAvgSmthR = 0.0;
			}
		}

		// float bin_sum[indiNum][numNeurons];
		float bin_sum_max = 0.0;
		float avgFR[numExc];
		float maxFR[indiNum];

		for (int i = 0; i < indiNum; i++) {
			maxFR[i] = 0.0;
			// for (int n = 0; n < numNeurons; n++) {
			// 	bin_sum[i][n] = 0.0;
			// }	
		}					

		for (int n = 0; n < numExc; n++) {
			avgFR[n] = 0.0;
		}

		string match_string = "MatchList";
		string match_stats_string = "MatchStats";
		string end = ".csv";
		// string ind_str = "Indi";

		// int chosen[indiNum][numNeurons] = {}; // which exc neuron is matched to this real neuron
		// int match[indiNum][numExc] = {}; // which real neuron is matched to this exc neuron
		vector<vector<int> > chosen(indiNum, vector<int>(numNeurons, -2));
		vector<vector<int> > match(indiNum, vector<int>(numExc, -2));

		vector<vector<float> > matchedCoef(indiNum, vector<float>(numNeurons, -1.0));

		if (forceMatch == true) {
			for (int i = 0; i < indiNum; i++) {
				strs2 << i;
				ch2 = strs2.str();
				strs2.str("");

				ifstream matchList;
				string name = (net_folder+ch2+driver_str+match_string+end);
				matchList.open(name.c_str());

				coefficient = 0;

				for (int n = 0; n < numNeurons; n++) {
					chosen[i][n] = 0;
					matchList >> chosen[i][n];
					chosen[i][n] = chosen[i][n] - 1;
					//cout << "Match to neuron " << n << " is " << chosen[n] << endl;
					matchedCoef[i][n] = corrCoef(trueTeFR[n], smoothedFRsAVG[i][chosen[i][n]], BINS_ALL);
					coefficient += matchedCoef[i][n];
				}
				fitness[i] = coefficient;
			}
				// cout << coefficient << endl;
		}
		else {
			// Get the avg corr coeff across all six paths between each neuron
			for (int i = 0; i < indiNum; i++) {
				for (int n = 0; n < numNeurons; n++) {
					chosen[i][n] = -2;
				}
				for (int n = 0; n < numExc; n++) {
					match[i][n] = (n+1);
				}
				coefficient = 0;

				// pairwise correlation for all real neurons and exc neurons
				for (int n = 0; n < numNeurons; n++) {
					for (int f = 0; f < numExc; f++) {
						avgCoefficient[n][f] = corrCoef(smoothedFRsAVG[i][f], trueTeFR[n], BINS_ALL);
					}
				}

				for (int n = 0; n < numNeurons; n++) {
					// if this real neuron is not matched yet
					if (chosen[i][n] == -2) {
						// go through corr between this real neuron and all exc neurons
						for (int f = 0; f < numExc; f++) {
							// find max corr in pool of exc neurons not yet matched
							if (match[i][f] != 0) {
								if (avgCoefficient[n][f] > bestCorr[i][n]) {
									bestCorr[i][n] = avgCoefficient[n][f];
									temp_index = f;
								}
							}
						}
						coefficient += bestCorr[i][n];

						chosen[i][n] = temp_index;
						if (temp_index != -2.0) {
							match[i][temp_index] = 0;
						}
						temp_index = -2.0;
					}
				}
				fitness[i] = coefficient;
			}
		}

		string odd_str = "Odd";
		string even_str = "Even";

		string fit_all_str = "fitness";
		string name_fit_all = (net_folder+fit_all_str+end);
		ofstream writeFitAll;
		writeFitAll.open(name_fit_all.c_str());

		// If on the last trial, calculate the fitness for this configuration.
		for (int i = 0; i < indiNum; i++) {
			// Calculate avg FR 
			for (int n = 0; n < numExc; n++) {
				bin_sum_max = 0;
				for (int b = 0; b < BINS_ALL; b++) {
					bin_sum_max += smoothedFRsAVG[i][n][b];
				}
				avgFR[n] = (bin_sum_max/BINS_ALL);

				if (avgFR[n] > maxFR[i]) {
					maxFR[i] = avgFR[n];
				}
				avgFR[n] = 0.0;
			}
			if (maxFR[i] > maxFR_thresh) {
				fitness[i] = fitness[i] - (maxFR[i]-maxFR_thresh);
			}

			ofstream writeAVGFRs;
			ofstream writeInpAVGFRs;
			ofstream writeAVGFRsEven;
			ofstream writeAVGFRsOdd;

			strs2 << i;
			ch2 = strs2.str();
			strs2.str("");
			
			string inp_fr_str = "InpFR";
			string lesion_str = "";

			if (lesionAV == 1) {
				lesion_str = lesion_str + "AV_";
			}
			if (lesionHD == 1) {
				lesion_str = lesion_str + "HD_";
			}
			if (lesionPos == 1) {
				lesion_str = lesion_str + "Pos_";
			}
			if (lesionLV == 1) {
				lesion_str = lesion_str + "LV_";
			}
			if (lesionExc == 1) {
				lesion_str = lesion_str + "Exc_";
			}
			if (lesionInh == 1) {
				lesion_str = lesion_str + "Inh_";
			}
			if (lesionAV == 0 && lesionHD == 0 && lesionPos == 0 && lesionLV == 0 && lesionExc == 0 && lesionInh == 0 && forceMatch == 0) {
				lesion_str = "Unlesioned_";
			}
			if (lesionAV == 0 && lesionHD == 0 && lesionPos == 0 && lesionLV == 0 && lesionExc == 0 && lesionInh == 0 && forceMatch == 1) {
				lesion_str = "Unlesioned_ForcedMatch_";
			}

			lesion_str = lesion_str.substr(0, lesion_str.length()-1)+"/";
			string lesion_folder = net_folder+ch2+driver_str+lesion_str;
			if (writeSimData) {
				if (mkdir(lesion_folder.c_str(), 0777) == -1) {
					 cerr << "Error :  " << strerror(errno) << endl; 
					 return;
				}
			}

			string fit_str = "FitnessScores";
			string name_fitscore = (net_folder+ch2+driver_str+lesion_str+fit_str+end);
			lesionFitScores.open(name_fitscore.c_str());

			lesionFitScores << seed << endl;
			lesionFitScores << maxFR[i] << endl;
			lesionFitScores << fitness[i] << endl;
			lesionFitScores.close();

			writeFitAll << lesion_str.substr(0, lesion_str.length()-1)+": " << fitness[i] << endl;

			outputStream << fitness[i] << endl;
			string mean_fr_str = "MeanFR";

			string name = (net_folder+ch2+driver_str+lesion_str+mean_fr_str+end);
			string name_even = (net_folder+ch2+driver_str+lesion_str+even_str+end);
			string name_odd = (net_folder+ch2+driver_str+lesion_str+odd_str+end);
			string inpName = (net_folder+ch2+driver_str+lesion_str+inp_fr_str+end);

			writeAVGFRs.open(name.c_str());
			writeInpAVGFRs.open(inpName.c_str());
			writeAVGFRsEven.open(name_even.c_str());
			writeAVGFRsOdd.open(name_odd.c_str());

			for (int n = 0; n < numExc; n++) {
				for (int b = 0; b < BINS_ALL; b++) {
					if (n < INP_N) {
						writeInpAVGFRs << inpSmthFRsAVG[n][b] << ",";
					}
					writeAVGFRs << smoothedFRsAVG[i][n][b] << ",";
					writeAVGFRsEven << smoothedFRsAVG_EVEN[i][n][b] << ",";
					writeAVGFRsOdd << smoothedFRsAVG_ODD[i][n][b] << ",";
				}
				if (n < INP_N) {
					writeInpAVGFRs << endl;
				}
				writeAVGFRs << endl;
				writeAVGFRsEven << endl;
				writeAVGFRsOdd << endl;
			}
			writeInpAVGFRs.close();
			writeAVGFRs.close();
			writeAVGFRsEven.close();
			writeAVGFRsOdd.close();

			ofstream writeFRs[totalTeTrials*repTrialsTe];
			// ofstream writeFRs;
			string rec_str = "Rec";
			string trial_str = "Trial";
			int t = 0;
			for (int r = 0; r < totalTeTrials; r++) {
				for (int z = 0; z < repTrialsTe; z++) {
					stringstream tx;
					tx << t;
					string ty = tx.str();
					name = (net_folder+ch2+driver_str+lesion_str+trial_str+ty+end);
					writeFRs[t].open(name.c_str());
					/** write input data (nTypesInput * numBinsAll) **/
					for (int b = 0; b < BINS_ALL; b++) {
						writeFRs[t] << inpsDataAll[r][z][0][b] << ",";
					}
					writeFRs[t] << endl;
					for (int b = 0; b < BINS_ALL; b++) {
						writeFRs[t] << inpsDataAll[r][z][1][b] << ",";
					}
					writeFRs[t] << endl;	
					for (int b = 0; b < BINS_ALL; b++) {
						writeFRs[t] << inpsDataAll[r][z][2][b] << ",";
					}
					writeFRs[t] << endl;	
					for (int b = 0; b < BINS_ALL; b++) {
						writeFRs[t] << inpsDataAll[r][z][3][b] << ",";
					}
					writeFRs[t] << endl;			
					for (int b = 0; b < BINS_ALL; b++) {
						writeFRs[t] << inpsDataAll[r][z][4][b] << ",";
					}
					writeFRs[t] << endl;
					/** finish writing input data **/ 	
					// write input firing rates (INP_N * numBinsAll)
					for (int n = 0; n < INP_N; n++) {
						for (int b= 0; b < BINS_ALL; b++) {
							writeFRs[t] << inpsSmthFRs[r][z][n][b] << ",";
						}
						writeFRs[t] << endl;
					}
					// write exc firing rates (nExc * numBinsAll)
					for (int n = 0; n < numExc; n++) {
						for (int b= 0; b < BINS_ALL; b++) {
							writeFRs[t] << smoothedFRs[i][r][z][n][b] << ",";
						}
						writeFRs[t] << endl;
					}
					writeFRs[t].close();
					t++;
				}
			}
			// write matched neurons and correlation scores
			if (forceMatch == true) {
				ofstream match_stats;
				name = (net_folder+ch2+driver_str+lesion_str+match_stats_string+end);
				match_stats.open(name.c_str());
				for (int n = 0; n < numNeurons; n++) {
					match_stats << (chosen[i][n] + 1) << ", " << matchedCoef[i][n] << endl;
				}	
				match_stats.close();
			} else {
				ofstream match_list;
				string name = (net_folder+ch2+driver_str+match_string+end);
				match_list.open(name.c_str());
				for (int n = 0; n < numNeurons; n++) {
					match_list << (chosen[i][n] + 1) << endl;
				}			
				match_list.close();

				ofstream match_stats;
				name = (net_folder+ch2+driver_str+lesion_str+match_stats_string+end);
				match_stats.open(name.c_str());
				for (int n = 0; n < numNeurons; n++) {
					match_stats << (chosen[i][n] + 1) << ", " << bestCorr[i][n] << endl;
				}	
				match_stats.close();
			}

		}

		if (loadFlag == 0) {
			string realFRstr = "realFR";
			string realfr_name = (net_folder+realFRstr+end);

		 ofstream fileRealFR;
		 fileRealFR.open(realfr_name.c_str());

			for (int n = 0; n < numNeurons; n++) {
		 	for (int b = 0; b < BINS_ALL; b++) {
					fileRealFR << trueTeFR[n][b] << ",";
				}
				fileRealFR << endl;
			}
		}
		// cout << endl << "End generation" << endl << endl;	

		// delete inpIn;

		// delete network;

		// //////////////////////////////////////--- DATA DELETION ---///////////////////////////////////////////////////////////////////////

		// for (int r = 0; r < totalTeTrials; r++) {
		// 	for (int j = 0; j < repTrialsTe; j++) {
		// 		for (int n = 0; n < numExc; n++) {
		// 			delete[] smoothedFRs[r][j][n];
		// 		}
		// 		delete[] smoothedFRs[r][j];
		// 	}
		// 	delete[] smoothedFRs[r];
		// }
		// delete[] smoothedFRs;

		// for (int r = 0; r < totalTeTrials; r++) {
		// 	for (int j = 0; j < repTrialsTe; j++) {
		// 		for (int n = 0; n <INP_N; n++) {
		// 			delete[] inpsSmthFRs[r][j][n];
		// 		}
		// 		delete[] inpsSmthFRs[r][j];
		// 	}
		// 	delete[] inpsSmthFRs[r];
		// }
		// delete[] inpsSmthFRs;
		
		// for (int n = 0; n < numExc; n++) {
		// 	delete[] smoothedFRsAVG[i][n];
		// }
		// delete[] smoothedFRsAVG;

		// for (int n = 0; n < INP_N; n++) {
		// 	delete[] inpSmthFRsAVG[n];
		// }
		// delete[] inpSmthFRsAVG;

		// for (int n = 0; n < numExc; n++) {
		// 	delete[] smoothedFRsAVG_EVEN[i][n];
		// }
		// delete[] smoothedFRsAVG_EVEN;

		// for (int n = 0; n < numExc; n++) {
		// 	delete[] smoothedFRsAVG_ODD[i][n];
		// }
		// delete[] smoothedFRsAVG_ODD;

		// for (int n = 0; n < numNeurons; n++) {
		// 	delete[] avgCoefficient[n];
		// }
		// delete[] avgCoefficient;

		// for (int n = 0; n < numNeurons; n ++){
		// 	delete[] trueTeFR[n];
		// }
		// delete[] trueTeFR;	
	}
};

int main(int argc, char* argv[]) {

	const SimMode simMode = CPU_MODE;
  	const LoggerMode verbosity = SILENT;

	const CogmapExperiment experiment(simMode, verbosity);
	const PTI pti(argc, argv, cout, cin);

	pti.runExperiment(experiment);

	return 0;
}
