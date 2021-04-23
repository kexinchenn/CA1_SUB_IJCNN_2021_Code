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

// C++ implementation of Matlab filtfilt function
#include "filtfilt.h"
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

	void run(const ParameterInstances &parameters, std::ostream &outputStream) const {
		bool loadFlag = 0;
		bool writeSimData = 0;

		bool LESION = 0;

		bool lesionAV = 0;
		bool lesionHD = 0;
		bool lesionPos = 0;
		bool lesionLV = 0;
		bool lesionExc = 0;
		bool lesionInh = 0;

		bool testing;
		bool forceMatch = 0;

		// int i = 7;
		int run = 10;

		srand(time(NULL));
		string dataRoot = "./data/";

		time_t timerTrial_start,timerTrial_end;
		time_t globalTimer_start, globalTimer_end;

		time(&globalTimer_start);

		// DATA LOADING VARIABLES
		const int totalTrTrials = 20;
		const int totalTeTrials = 20;
		const int repTrialsTr = 5;
		const int repTrialsTe = 5;
		// const int totalTrTrials = 1;
		// const int totalTeTrials = 1;
		// const int repTrialsTr = 1;
		// const int repTrialsTe = 1;

		const int numNeurons = 295;
		const int numRecordings = 32;
		// const int numNeurons = 1;
		// const int numRecordings = 1;
		const int numPaths = 6;

		int r_min = 0;

		int INP_ANG = 12;
		int INP_HD = 8;
		int INP_POS_X = 25;
		int INP_POS_Y = 18;
		int INP_POS = INP_POS_X * INP_POS_Y;
		// int INP_POS = 394;
		int INP_LV = 12;
		int INP_N = INP_ANG+INP_HD+INP_POS+INP_LV;
		// int INP_PROG = 25;

		// int INP_2 = 20;
		// int INP_3 = 410;
		// int INP_4 = 572;
	
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
		// float modifierMSE = 1;

		//POSITION ON TRACK PARAMETERS:
		float xCoordPref[INP_POS];
		float yCoordPref[INP_POS];
		// vector<float> xCoordPref;
		// vector<float> yCoordPref;

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
		string dataDir = (dataRoot + "data_HPC_train_test_20201204/");

		// number of trials for each path in each recording (numRecording x numPath)
		string trialCountFileTrain = (dataDir + "trial_count_train.csv");
		vector<vector<float> > trialCountTrain; 
		trialCountTrain = loadData(trialCountFileTrain, numRecordings, numPaths);

		string trialCountFileTest = (dataDir + "trial_count_test.csv");
		vector<vector<float> > trialCountTest; 
		trialCountTest = loadData(trialCountFileTest, numRecordings, numPaths);

		// template XY position to create input place cells (2 x numPosXY)
		// string templatePosXYFile = (dataDir + "allPathTemplateXYDownsampledBy5.csv");
		// vector<vector<float> > templatePosXY; 
		// templatePosXY = loadData(templatePosXYFile, 2, INP_POS);
		// xCoordPref = templatePosXY[0];
		// yCoordPref = templatePosXY[1];

		// ofstream templateFile;
		// templateFile.open("templateXY.csv", ofstream::out | ofstream::app);

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
		getRealTestFRs(numNeurons, numPaths, trueTeFR);

		// -------------------------------- DONE READING IN ALL DATA ------------------------------------------------

		int order;

		float r = ((float) rand() / (RAND_MAX)) * 1;
		if (r < 0.5) {
			order = 0;
		}
		else {
			order = 2;
		}		
			
		int repstart1;
		int repstart2;
		int repend1;
		int repend2; 

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
		if (indiNum > 1) {
			writeSimData = 0;
		}
		int runTimeSec = 0;
		int runTimeMs = 15;

		float fitness[indiNum];

		int exc[indiNum];
		int inh[indiNum];

		SpikeMonitor* excMonitor[indiNum];
		SpikeMonitor* inhMonitor[indiNum];

		ConnectionMonitor* inpToExc;
		ConnectionMonitor* inpToInh;
		ConnectionMonitor* inhToExc;
		ConnectionMonitor* excToExc;

		// EXPERIMENT VARIABLES: 
		int trialNumTr;
		int recNum;

		int exc2D;
		int inh2D;
		int inp2D;

		int count = 0;
		
		int ithGPU = 0;

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
		

		float ****inpsSmthFRs = new float***[totalTeTrials];
		for (int r = 0; r < totalTeTrials; r++) {
			inpsSmthFRs[r] = new float**[repTrialsTe];
			for (int j = 0; j < repTrialsTe; j++) {
				inpsSmthFRs[r][j] = new float*[5];
				for (int n = 0; n < 5; n++) {
					inpsSmthFRs[r][j][n] = new float[BINS_ALL];
					for (int b = 0; b < BINS_ALL; b++) {
						inpsSmthFRs[r][j][n][b] = 0.0;
					}
				}
			}
		}

		float **avgCoefficient;
		avgCoefficient = new float*[numNeurons];
		for (int n = 0; n < numNeurons; n++) {
			avgCoefficient[n] = new float[numExc];
		}

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

		int p_max = 1;
		int p = 0;

		if (forceMatch && LESION) {
			p_max = 32;
			p = 31;
		}
		// while (p < p_max) {	
		// cout << p << endl;
		// cout << p_max << endl;
		
		if (forceMatch && LESION) {
			getCombination(p, lesionAV, lesionHD, lesionPos, lesionLV, lesionExc, lesionInh);
		}
		p ++;
		// cout << "Lesion AV is: " << lesionAV << endl;
		// cout << "lesion HD is: " << lesionHD << endl;
		// cout << "Lesion position is: " << lesionPos << endl;
		// cout << "Lesion LV is: " << lesionLV << endl;
		// cout << "Lesion exc is: " << lesionExc << endl;
		// cout << "Lesion inh is: " << lesionInh << endl;	

		CARLsim* network;

		if (loadFlag == 0) {
			network = new CARLsim("direct_driver", simMode, verbosity, time(NULL));
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
			time_t seed = 1462301880; // NETWORK 36 - 10 - 7
			network = new CARLsim("loadSim", simMode, verbosity, time(NULL));
		}

		float inpToExcMaxWght, inpToInhMaxWght, excToExcMaxWght, inhToExcMaxWght;
		float excBaseFR, inhBaseFR;
		float ALPHA_LTP_EE, ALPHA_LTD_EE, TAU_LTP_EE, TAU_LTD_EE;
		float ALPHA_LTP_EI, ALPHA_LTD_EI, TAU_LTP_EI, TAU_LTD_EI;
		float ALPHA_LTP_IE, ALPHA_LTD_IE, TAU_LTP_IE, TAU_LTD_IE;
	
		// int inpExc[indiNum];
		int inputs;
		float initWtScale = 0.5;

		short int inpExcToExc[indiNum];
		// short int inpExcToInh[indiNum];
		// short int excToInpExc[ind
		short int inpExcToInhib[indiNum];
		short int inhibToExcit[indiNum];
		short int excitToExcit[indiNum];
		inputs = network->createSpikeGeneratorGroup("inputs", INP_N, EXCITATORY_NEURON);

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

			inhBaseFR = parameters.getParameter(i,12);
			excBaseFR = parameters.getParameter(i,19);

			inpToExcMaxWght = parameters.getParameter(i,13);
			inpToInhMaxWght = parameters.getParameter(i,14);
			excToExcMaxWght = parameters.getParameter(i,15);
			inhToExcMaxWght = parameters.getParameter(i,16);

			HOMEO_AVERAGE_TIME_SCALE_EXC = parameters.getParameter(i,17);
			HOMEO_AVERAGE_TIME_SCALE_INH = parameters.getParameter(i,18);

			exc[i] = network->createGroup("excitatory", numExc, EXCITATORY_NEURON);
			network->setNeuronParameters(exc[i], REG_IZH[0], REG_IZH[1], REG_IZH[2], REG_IZH[3]);
			inh[i] = network->createGroup("inhibitory", numInh, INHIBITORY_NEURON);
			network->setNeuronParameters(inh[i], FAST_IZH[0], FAST_IZH[1], FAST_IZH[2], FAST_IZH[3]);
			network->setConductances(true,COND_tAMPA,COND_tNMDA,COND_tGABAa,COND_tGABAb);

			inpExcToExc[i] = network->connect(inputs, exc[i], "random", 
				RangeWeight(0.0f, inpToExcMaxWght*initWtScale, inpToExcMaxWght), 0.1f, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
			excitToExcit[i] = network->connect(exc[i], exc[i], "random", 
				RangeWeight(0.0f, excToExcMaxWght*initWtScale, excToExcMaxWght), 0.1f, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
			inpExcToInhib[i] = network->connect(inputs, inh[i], "random", 
				RangeWeight(0.0f, inpToInhMaxWght*initWtScale, inpToInhMaxWght), 0.1f, RangeDelay(1), RadiusRF(-1), SYN_PLASTIC);
			inhibToExcit[i] = network->connect(inh[i], exc[i], "random", 
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
		// // OBVIOUSLY THIS NEEDS TO BE CLEANED UP
		if (loadFlag == 1) {
			simFid = NULL;
			simFid = fopen("sim_direct_driver.dat", "rb");
			// cout << "Loading previous network" << std::endl;
			network->loadSimulation(simFid);
		}

		network->setupNetwork();

		if (loadFlag == 1) {
			fclose(simFid);
		}

		// int lesion_start;
		// int lesion_end;
		
		// if (lesionAV == true) {
		// 	lesion_start = 0;
		// 	lesion_end = 12;
		// }
		// else if (lesionHD == true) {
		// 	lesion_start = 12;
		// 	lesion_end = 20;
		// }
		// else if (lesionPos == true) {
		// 	lesion_start = 20;
		// 	lesion_end = 410;
		// }
		// else if (lesionLV == true) {
		// 	lesion_start = 410;
		// 	lesion_end = 417;
		// }

		// if (lesionAV == 1) {
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inputs); preId++) {
		// 		if (preId >= 0 && preId < 12) {
		// 			for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
		// 				float w = 0.0;
		// 				// cout << w << std::endl;
		// 				network->setWeight(inpExcToExc[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inputs); preId++) {
		// 		if (preId >= 0 && preId < 12) {
		// 			for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
		// 				float w = 0.0;
		// 				network->setWeight(inpExcToInhib[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}			
		// }
		// if (lesionHD == 1) {
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inputs); preId++) {
		// 		if (preId >= 12 && preId < 20) {
		// 			for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
		// 				float w = 0.0;
		// 				// cout << w << std::endl;
		// 				network->setWeight(inpExcToExc[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inputs); preId++) {
		// 		if (preId >= 12 && preId < 20) {
		// 			for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
		// 				float w = 0.0;
		// 				network->setWeight(inpExcToInhib[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}			
		// }
		// if (lesionPos == 1) {
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inputs); preId++) {
		// 		if (preId >= 20 && preId < 410) {
		// 			for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
		// 				float w = 0.0;
		// 				// cout << w << std::endl;
		// 				network->setWeight(inpExcToExc[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inputs); preId++) {
		// 		if (preId >= 20 && preId < 410) {
		// 			for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
		// 				float w = 0.0;
		// 				network->setWeight(inpExcToInhib[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}			
		// }
		// if (lesionLV == 1) {
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inputs); preId++) {
		// 		if (preId >= 410 && preId < 417) {
		// 			for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
		// 				float w = 0.0;
		// 				// cout << w << std::endl;
		// 				network->setWeight(inpExcToExc[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inputs); preId++) {
		// 		if (preId >= 410 && preId < 417) {
		// 			for (int postId = 0; postId < network->getGroupNumNeurons(inh[i]); postId++) {
		// 				float w = 0.0;
		// 				network->setWeight(inpExcToInhib[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}			
		// }									
		// if (lesionExc == 1) {
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(exc[i]); preId++) {
		// 		for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
		// 			if (preId != postId) {
		// 				float w = 0.0;
		// 				network->setWeight(excitToExcit[i], preId, postId, w); // between 0 and 1
		// 			}
		// 		}
		// 	}		
		// }
		// if (lesionInh == 1) {
		// 	for (int preId = 0; preId < network->getGroupNumNeurons(inh[i]); preId++) {
		// 		for (int postId = 0; postId < network->getGroupNumNeurons(exc[i]); postId++) {
		// 			float w = 0.0;
		// 			network->setWeight(inhibToExcit[i], preId, postId, w); // between 0 and 1
		// 		}
		// 	}		
		// }

		string conn = "results/conn_";
		string conn1 = "inp_exc";
		string conn2 = "inp_inh";
		// string conn3 = "inp_inh_";
		string conn3 = "inh_exc";
		string conn4 = "exc_exc";

		string ind = "results/spk_";
		string ind1 = "excitatory_";
		string ind2 = "inhibitory_";
		string ind3 = "excInps_";
		string e = ".dat";

		for (int i = 0; i < indiNum; i++) {
			// stringstream strs1;
			// strs1 << i;
			// string = strs1.str();
			// string allConnInpToExc = (conn+conn1+e);
			// inpToExc = network->setConnectionMonitor(inputs, exc[i], allConnInpToExc.c_str());
			// inpToExc->setUpdateTimeIntervalSec(-1);

			// string allConnInpToInh = (conn+conn2+e);
			// inpToInh = network->setConnectionMonitor(inputs, inh[i], allConnInpToInh.c_str());
			// inpToInh->setUpdateTimeIntervalSec(-1);

			// string allConnInhToExc = (conn+conn3+e);
			// inhToExc = network->setConnectionMonitor(inh[i], exc[i], allConnInhToExc.c_str());
			// inhToExc->setUpdateTimeIntervalSec(-1);

			// string allConnExcToExc = (conn+conn4+e);
			// excToExc = network->setConnectionMonitor(exc[i], exc[i], allConnExcToExc.c_str());
			// excToExc->setUpdateTimeIntervalSec(-1);

			string allExc = (ind+ind1+e);
			excMonitor[i] = network->setSpikeMonitor(exc[i], "NULL");
			string allInh = (ind+ind2+e);
			inhMonitor[i] = network->setSpikeMonitor(inh[i], "NULL");

		}

		vector<vector<float> > posXY;
		// vector<float> posY;
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

		r = ((float) rand() / (RAND_MAX)) * 1;
		if (r < 0.5) {
			order = 0;
		}
		else {
			order = 2;
		}

		int numTeTrials;

		int mod;
		int trialNumTe;

		stringstream strs2;
		strs2 << run;
		string ch2 = strs2.str();
		string net_folder = "./";
		string driver_str = "_Driver/";
		string recstr = "recordings.csv";
		string repstr = "repetitions.csv";

		string rec_name = (net_folder+ch2+driver_str+recstr);
		string rep_name = (net_folder+ch2+driver_str+repstr);

		string simDataDir = (net_folder+ch2+driver_str);
		string makeDirCommand = "mkdir -p ";

		if (writeSimData) {
			if (mkdir(simDataDir.c_str(), 0777) == -1) {
				 cerr << "Error :  " << strerror(errno) << endl; 
				 return;
			}
		}

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
				for (int r = 0; r < 4; r++) {
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
		testing = 0;
		for (int total = 0; total < totalTrTrials; total++) {
			if (forceMatch == 0) {
				recNum = (rand() % numRecordings) +0;
				// writeRecordings << recNum << ",";
			}
			else {
				recNum = recordings[total];
				//cout << recNum << endl;
			}

			time(&trialTimer_start);
			for (int reps = 0; reps < numPaths; reps++) {
				if (reps < 4) {
					numBins = BINS_IN;
					startBinInd = reps * BINS_IN;
				} else {
					numBins = BINS_OUT;
					startBinInd = BINS_IN * 4 + BINS_OUT * (reps-4);
				}
				for (int z = 0; z < repTrialsTe; z++) {
					if (forceMatch == 0) {
						numTeTrials = trialCountTrain[recNum][reps];							

						if (numTeTrials > 0) {
							trialNumTe = (rand() % numTeTrials);
						}
						else {
							trialNumTe = 0;
						}
						// writeRepetitions << trialNumTe << ",";					
					}
					else {
						trialNumTe = repetitions[total][reps][z];
						//cout << trialNumTe << endl;
					}

					posXY = getXYData(recNum, reps, testing, trialNumTe, numBins);						
					coordsFR = gauss2dXY(posXY[0], posXY[1], xCoordPref, yCoordPref, INP_POS_X, INP_POS_Y, numBins, sigma_x, r_max);						
					// coordsFR = gauss2d(posXY[0], posXY[1], xCoordPref, yCoordPref, INP_POS, numBins, sigma_x, r_max);						
					angVel = getAngVel(recNum, reps, testing, trialNumTe, numBins);
					angVelFR = gauss1d(angVel, angVelPref, INP_ANG, numBins, sigmaAngVel, r_max_ang_lin);
					linVel = getLinVel(recNum, reps, testing, trialNumTe, numBins);
					linVelFR = gauss1d(linVel, linVelPref, INP_LV, numBins, sigmaLinVel, r_max_ang_lin);
					headDir = getHeadDir(recNum, reps, testing, trialNumTe, numBins);
					headDirFR = cosine(headDir, headDirPref, INP_HD, numBins, r_max, r_min);
					// occ = getOcc(recNum, reps, testing, trialNumTe, numBins);
					time(&timerTrial_start);

					// BEGIN BINS LOOP
					//time(&timerBin_start);
					for (int b = 0; b < numBins; b++){
						PoissonRate inpIn(INP_N);
						count = 0;
						for (int j = 0; j < INP_ANG; j++) {
							inpIn.setRate(count, angVelFR[j][b]);
							count ++;
						}
						// count = 0;
						for (int j = 0; j < INP_HD; j++) {
							inpIn.setRate(count, headDirFR[j][b]);
							count ++;
						}
						// count = 0;
						for (int j = 0; j < INP_POS; j++) {
							inpIn.setRate(count, coordsFR[j][b]);
							count ++;
						}
						// count = 0;
						for (int j = 0; j < INP_LV; j ++){
							inpIn.setRate(count, linVelFR[j][b]);
							count ++;
						}
						// Start recording...
						network->setSpikeRate(inputs, &inpIn);

						// runTimeMs = occ[b] / 60 * 1000;
						// runTimeSec = int(floor(runTimeMs / 1000));
						// runTimeMs = runTimeMs % 1000;

						network->runNetwork(runTimeSec,runTimeMs);
					}
				}
			}
		}
		network->startTesting();

		/********************************************************************************
									Testing
		********************************************************************************/		
		// Loop over every recording
		testing = 1;
		for (int total = 0; total < totalTeTrials; total++) {
			if (forceMatch == 0) {
				recNum = (rand() % numRecordings) +0;
				writeRecordings << recNum << ",";
			}
			else {
				recNum = recordings[total];
				//cout << recNum << endl;
			}

			time(&trialTimer_start);
			for (int reps = 0; reps < numPaths; reps++) {
				if (reps < 4) {
					numBins = BINS_IN;
					startBinInd = reps * BINS_IN;
				} else {
					numBins = BINS_OUT;
					startBinInd = BINS_IN * 4 + BINS_OUT * (reps-4);
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

					posXY = getXYData(recNum, reps, testing, trialNumTe, numBins);						
					coordsFR = gauss2dXY(posXY[0], posXY[1], xCoordPref, yCoordPref, INP_POS_X, INP_POS_Y, numBins, sigma_x, r_max);						
					// coordsFR = gauss2d(posXY[0], posXY[1], xCoordPref, yCoordPref, INP_POS, numBins, sigma_x, r_max);						
					angVel = getAngVel(recNum, reps, testing, trialNumTe, numBins);
					angVelFR = gauss1d(angVel, angVelPref, INP_ANG, numBins, sigmaAngVel, r_max_ang_lin);
					linVel = getLinVel(recNum, reps, testing, trialNumTe, numBins);
					linVelFR = gauss1d(linVel, linVelPref, INP_LV, numBins, sigmaLinVel, r_max_ang_lin);
					headDir = getHeadDir(recNum, reps, testing, trialNumTe, numBins);
					headDirFR = cosine(headDir, headDirPref, INP_HD, numBins, r_max, r_min);
					// occ = getOcc(recNum, reps, testing, trialNumTe, numBins);

					time(&timerTrial_start);
					// BEGIN BINS LOOP
					//time(&timerBin_start);
					vector<vector<vector<float> > > excFRs(indiNum, vector<vector<float> >(numBins, vector<float>(numExc, 0)));

					for (int b = 0; b < numBins; b++){
						PoissonRate inpIn(INP_N);
						inpsSmthFRs[total][z][0][b+startBinInd] = angVel[b];
						inpsSmthFRs[total][z][1][b+startBinInd] = headDir[b];
						inpsSmthFRs[total][z][2][b+startBinInd] = posXY[0][b];
						inpsSmthFRs[total][z][3][b+startBinInd] = posXY[1][b];
						inpsSmthFRs[total][z][4][b+startBinInd] = linVel[b];
						// inpsSmthFRs[total][z][5][b+startBinInd] = occ[b];

						count = 0;
						for (int j = 0; j < INP_ANG; j++) {
							inpIn.setRate(count, angVelFR[j][b]);
							count ++;
						}
						// count = 0;
						for (int j = 0; j < INP_HD; j++) {
							inpIn.setRate(count, headDirFR[j][b]);
							count ++;
						}
						// count = 0;
						for (int j = 0; j < INP_POS; j++) {
							inpIn.setRate(count, coordsFR[j][b]);
							count ++;
						}
						// count = 0;
						for (int j = 0; j < INP_LV; j ++){
							inpIn.setRate(count, linVelFR[j][b]);
							count ++;
						}
						// count = 0;
						// for (int j = INP_4; j < INP_N; j++){
						// 	if (reps == 0 || reps == 2) { // 0 = ALRL; 2 = BLRL
						// 		inpIn.setRate(j, progFR_LRL[b][count]);
						// 		count ++;
						// 	}
						// 	else if (reps == 1 || reps == 3) { // 1 = ARLR; 3 = BRLR
						// 		inpIn.setRate(j, progFR_RLR[b][count]);
						// 		count ++;									
						// 	}
						// }

						// Start recording...
						network->setSpikeRate(inputs, &inpIn);

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
								excFRs[i][b][index] = (exc2D * 1000.0) / recTime;
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
							// inhMonitor[i]->clear();
							// }

							// If we are on the last bin, then we can evaluate the fitness for a given neuron.
												// If we are on the last bin, we can smooth and store the data for each neuron.
							if (b == (numBins-1)) { //&& (z == 19) ) {
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
					// if (forceMatch == 0) {
					// 	writeRepetitions << endl;
					// }
				} // end trial loop
				if (forceMatch == 0) {
					writeRepetitions << endl;
				}
				time(&trialTimer_end);
				// cout << "Finished testing trial " << total+1 << "/" << totalTeTrials << " with time " 
					// << difftime(trialTimer_end, trialTimer_start) << endl;
				time(&trialTimer_start);
			} 

			if (forceMatch == 0) {
				writeRecordings.close();
				writeRepetitions.close();
			}
		} // END TESTING

		network->stopTesting();

		//TAKING SNAPSHOTS!
		// cout << "Taking weight snapshots" << endl;
		// inpToExc->takeSnapshot();
		// inpToInh->takeSnapshot();
		// inhToExc->takeSnapshot();
		// excToExc->takeSnapshot();
		// cout << "Finished snapshots" << endl;

		// BEGIN FITNESS CALC
		// float bestCorr = -1.0;
		vector<vector<float> > bestCorr(indiNum, vector<float>(numNeurons, -1.0));
		int temp_index = -2;
		float coefficient = 0.0;

		int chosen[numNeurons];
		int match[numExc];

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
							for (int inpType = 0; inpType < 5; inpType++) {
								inpAvgSmooth += inpsSmthFRs[r][z][inpType][b];
							}
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

		float bin_sum[indiNum][numNeurons];
		float bin_sum_max = 0.0;
		float avgFR[numExc];
		float maxFR[indiNum];

		float MSE_N[indiNum];
		float avgMSE[indiNum];
		for (int i = 0; i < indiNum; i++) {
			MSE_N[i] = 0.0;
			avgMSE[i] = 0.0;
			maxFR[i] = 0.0;
			for (int n = 0; n < numNeurons; n++) {
				bin_sum[i][n] = 0.0;
			}	
		}					

		for (int n = 0; n < numExc; n++) {
			avgFR[n] = 0.0;
		}

		string match_string = "_Driver/MatchList_Indi";
		string end = ".csv";

		if (forceMatch == true) {
			ifstream matchList;
			string name = (net_folder+ch2+match_string+end);
			matchList.open(name.c_str());
			int match[numNeurons];
			for (int n = 0; n < numNeurons; n++) {
				match[n] = 0;
				matchList >> match[n];
				match[n] = match[n] - 1;
				//cout << "Match to neuron " << n << " is " << match[n] << endl;
			}
			for (int i = 0; i < indiNum; i++) {
				for (int n = 0; n < numNeurons; n++) {
					//cout << "Match to neuron " << n << " is " << match[n] << endl;
					coefficient += corrCoef(trueTeFR[n], smoothedFRsAVG[i][match[n]], BINS_ALL);
				}
			}
				// cout << coefficient << endl;
		}
		else {
			// Get the avg corr coeff across all four positions between each neuron
			for (int i = 0; i < indiNum; i++) {
				for (int n = 0; n < numNeurons; n++) {
					chosen[n] = -2;
				}
				for (int n = 0; n < numExc; n++) {
					match[n] = (n+1);
				}
				coefficient = 0;

				for (int n = 0; n < numNeurons; n++) {
					for (int f = 0; f < numExc; f++) {
						avgCoefficient[n][f] = corrCoef(smoothedFRsAVG[i][f], trueTeFR[n], BINS_ALL);
					}
				}

				// for (int i = 0; i < indiNum; i++) {
				// bestCorr = -1.0;
				for (int n = 0; n < numNeurons; n++) {
					if (chosen[n] == -2) {
						for (int f = 0; f < numExc; f++) {
							if (match[f] != 0) {
								if (avgCoefficient[n][f] > bestCorr[i][n]) {
									bestCorr[i][n] = avgCoefficient[n][f];
									temp_index = f;
								}
							}
						}
						coefficient += bestCorr[i][n];

						chosen[n] = temp_index;
						if (temp_index != -2.0) {
							match[temp_index] = 0;
						}

						// bin_sum[i][n] = 0.0;
						// for (int b = 0; b < BINS_ALL; b++) {
						// 	float diff = (smoothedFRsAVG[temp_index][b] - trueTeFR[n][b]);
						// 	float sq_diff = powf(diff, 2.0);
						// 	bin_sum[i][n] += sq_diff;
						// }
						// MSE_N[i] = sqrt( (bin_sum[i][n]/BINS_ALL) );

						// bestCorr = -1.0;
						temp_index = -2.0;
					}
				}
				//	MSE_N[i] = (MSE_N[i] / numNeurons);
				fitness[i] = coefficient;

			}
							
			// Calculate avg FR (within trial)
			for (int i = 0; i < indiNum; i++) {
				for (int n = 0; n < numExc; n++) {
					bin_sum_max = 0;
					for (int b = 0; b < BINS_ALL; b++) {
						bin_sum_max += smoothedFRsAVG[i][n][b];
					}
					avgFR[n] = (bin_sum_max/BINS_ALL);

					// bin_sum_max = 0;
					// for (int r = 0; r < totalTeTrials; r++) {
					// 	for (int z = 0; z < repTrialsTe; z++) {
					// 		for (int b = 0; b < BINS_ALL; b++) {
					// 			if (smoothedFRs[i][r][z][n][b] > bin_sum_max) {
					// 				bin_sum_max = smoothedFRs[i][r][z][n][b];
					// 			}
					// 		}
					// 	}
					// }
					// avgFR[n] = bin_sum_max;

					if (avgFR[n] > maxFR[i]) {
						maxFR[i] = avgFR[n];
					}
					avgFR[n] = 0.0;
				}
			}
		}

		// float bestFitness = -10000.0;
		// int punishFR = -500;

		// If on the last trial, calculate the fitness for this configuration.
		for (int i = 0; i < indiNum; i++) {
			// fitness[i] = coefficient;
			if (maxFR[i] > maxFR_thresh) {
				fitness[i] = fitness[i] - (maxFR[i]-maxFR_thresh);
			}

			lesionFitScores.open("UnlesionedFitnessScores.csv", std::ofstream::out | std::ofstream::app);
			lesionFitScores << maxFR[i] << endl;
			lesionFitScores << fitness[i] << endl;
			lesionFitScores.close();

			outputStream << fitness[i] << endl;

			// if (fitness > bestFitness) {
			// 	bestFitness = fitness;
			// 	i = i;
			// }

			//cout << "Sum of partial matches: " << coefficient[i] << endl;
	//		cout << "Avg MSE = " << -1*(MSE_N[i]*modifierMSE) << endl;
			//cout << "Max Avg FR: " << maxFR[i] << endl;
			// cout << "FITNESS TOTAL: " << fitness << endl << endl;

			ofstream writeAVGFRs;
			ofstream writeInpAVGFRs;
			ofstream writeAVGFRsEven;
			ofstream writeAVGFRsOdd;
			
			string inp_fr_str = "InpFR_";
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

			string odd_str = "_Odd";
			string even_str = "_Even";
			string ind_str = "Indi";

			string name = (net_folder+ch2+driver_str+lesion_str+ind_str+end);
			string name_even = (net_folder+ch2+driver_str+lesion_str+ind_str+even_str+end);
			string name_odd = (net_folder+ch2+driver_str+lesion_str+ind_str+odd_str+end);
			string inpName = (net_folder+ch2+driver_str+inp_fr_str+lesion_str+ind_str+end);

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

			// ofstream writeFRs[totalTeTrials*repTrialsTe];
			ofstream writeFRs;
			string rec_str = "_Rec";
			string trial_str = "_Trial";
			int t = 0;
			for (int r = 0; r < totalTeTrials; r++) {
				for (int z = 0; z < repTrialsTe; z++) {
					stringstream tx;
					tx << t;
					string ty = tx.str();
					name = (net_folder+ch2+driver_str+lesion_str+ind_str+trial_str+ty+end);
					writeFRs.open(name.c_str());
					// for (int b = 0; b < BINS_ALL; b++) {
					// 	writeFRs[t] << inpsSmthFRs[r][z][0][b] << ",";
					// }
					// writeFRs[t] << endl;
					// for (int b = 0; b < BINS_ALL; b++) {
					// 	writeFRs[t] << inpsSmthFRs[r][z][1][b] << ",";
					// }
					// writeFRs[t] << endl;	
					// for (int b = 0; b < BINS_ALL; b++) {
					// 	writeFRs[t] << inpsSmthFRs[r][z][2][b] << ",";
					// }
					// writeFRs[t] << endl;	
					// for (int b = 0; b < BINS_ALL; b++) {
					// 	writeFRs[t] << inpsSmthFRs[r][z][3][b] << ",";
					// }
					// writeFRs[t] << endl;			
					// for (int b = 0; b < BINS_ALL; b++) {
					// 	writeFRs[t] << inpsSmthFRs[r][z][4][b] << ",";
					// }
					// writeFRs[t] << endl;	
					for (int i = 0; i < indiNum; i ++) {																										
						for (int n = 0; n < numExc; n++) {
							for (int b= 0; b < BINS_ALL; b++) {
								writeFRs << smoothedFRs[i][r][z][n][b] << ",";
							}
							writeFRs << endl;
						}
					}
					writeFRs.close();
					t++;
				}
			}

			if (forceMatch == 0) {
				ofstream match_list;
				name = (net_folder+ch2+match_string+end);
				match_list.open(name.c_str());
				for (int n = 0; n < numNeurons; n++) {
					match_list << (chosen[n] + 1) << ", " << bestCorr[i][n] << endl;
				}			
				match_list.close();
			}

		}
	 //    std::ofstream fileRealFR;
	 //    fileRealFR.open("RealMeanFR.csv", std::ofstream::out | std::ofstream::app);

		// for (int n = 0; n < numNeurons; n++) {
	 //    	for (int b = 0; b < BINS_ALL; b++) {
		// 		fileRealFR << trueTeFR[n][b] << ",";
		// 	}
		// 	fileRealFR << endl;
		// }
		// cout << endl << "End generation" << endl << endl;	
		// network->stopTesting();

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
	const PTI pti(argc, argv, std::cout, std::cin);

	pti.runExperiment(experiment);

	return 0;
}