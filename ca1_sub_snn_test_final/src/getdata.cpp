#include <getdata.h>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

vector<vector<float> > getXYData (string dataDir, int rec, int posTrajType, bool testing, int trial, int bins) {
    vector<vector<float> > out_arr(2, vector<float>(bins, 0.0f));

    std::ofstream fileErrors;
    fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);
    std::string line, csvItem;

    // std::string filex2("Rec");
    std::string end(".csv");
    std::string stimx("Stim_");
    std::string stimy("Stim_");

    std::string combType;
    std::string trialType;

    stringstream strs2;
    strs2 << (posTrajType+1);
    combType = "Path_"+strs2.str()+"_";
    if (testing == 0) {
        trialType = "Tr";
    }
    else {
        trialType = "Te";
    }

    std::string filex(dataDir+"posX/posX");
    std::string filey(dataDir+"posY/posY");

    std::ifstream xCoord;
    std::ifstream yCoord;

    stringstream strs1;
    strs1 << (rec+1);
    string ch1 = strs1.str();

    std::string tempx = (filex+trialType+stimx+combType+ch1+end);
    std::string tempy = (filey+trialType+stimy+combType+ch1+end);

    xCoord.open(tempx.c_str());
    if (xCoord.fail()){
        fileErrors << tempx.c_str() << std::endl;
        fileErrors << trialType << " X Coord " << posTrajType << " " << (rec+1) << " failed" << std::endl;
    }

    yCoord.open(tempy.c_str());
    if (yCoord.fail()){
        fileErrors << tempy.c_str() << std::endl;
        fileErrors << trialType << " Y Coord " << posTrajType << " " << (rec+1) << " failed" << std::endl;
    }

    if (posTrajType >= 0 && posTrajType < 6) {
        int c = 0;
        float dummy;
        for (int l = 0; getline(xCoord, line) && l <= trial; l++) {
            if (l == trial) {
                istringstream myline(line);
                while(getline(myline, csvItem, ',')) {
                    std::istringstream str(csvItem);
                    str >> dummy;
                    if (std::isnan(dummy) || std::isinf(dummy)) {
                        dummy = 0.0f;
                    }
                    out_arr[0][c] = dummy;
                    // fileErrors << "String: " << tempx << endl;
                    // fileErrors << out_arr_x[c] << endl;
                    // fileErrors << "Rec: " << rec << endl;
                    // fileErrors << "posTrajType: " << posTrajType << endl;
                    // fileErrors << "trial : " << trial << endl;                   
                    c++;
                    str.clear();
                }
            }
        }

        c = 0;
        for (int l = 0; getline(yCoord, line) && l <= trial; l++) {
            if (l == trial) {
                istringstream myline(line);
                while(getline(myline, csvItem, ',')) {
                    std::istringstream str(csvItem);
                    str >> dummy;
                    if (std::isnan(dummy) || std::isinf(dummy)) {
                        dummy = 0.0f;
                    }
                    out_arr[1][c] = dummy;
                    c++;
                    str.clear();
                }
            }
        }
    }
    else {
        fileErrors << "Wrong pos/traj code" << endl;
    }

    xCoord.close();
    yCoord.close();
    
    fileErrors.close();
    return out_arr;
}


void getRealTestFRs (string dataDir, int numNeurons, int numNeuronsAll, int numPaths, float** out_arr) {

    string inputVar = "meanLR";

    string dataFile;
    vector<vector<float> > data;

    // string ch1;
    // stringstream strs1;
    // std::string end(".csv");

    int numBins;
    int numBinsIn = 140;
    int numBinsOut = 197;
    int startBinInd;
    int pathInd;

    for (int path = 0; path < numPaths; path ++) {

        if (path < 4) {
            numBins = numBinsIn;
            startBinInd = path * numBinsIn;
        } else {
            numBins = numBinsOut;
            startBinInd = numBinsIn * 4 + numBinsOut * (path-4);
        } 
        pathInd = path+1;

        // load data from one path
        dataFile = getFileName(dataDir, inputVar, pathInd);
        data = loadData(dataFile, numNeuronsAll, numBins); 

        for (int neur = 0; neur < numNeurons; neur ++) {
            for (int bin = 0; bin < numBins; bin ++) {
                if (isnan(data[neur][bin]) || isinf(data[neur][bin])) {
                    data[neur][bin] = 0.0f;
                }
                out_arr[neur][startBinInd+bin] = data[neur][bin];
            }
        }
    }
}

vector<float> getAngVel (string dataDir, int rec, int posTrajType, bool testing, int trial, int bins) {
    vector<float> out_arr(bins, 0.0f);

    std::ofstream fileErrors;
    fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);
    std::string line, csvItem;

    std::ifstream angVel;

    std::string fileAngVel(dataDir+"angVel/angVel");
    std::string stim("Stim_");

    std::string combType;
    std::string trialType;
    stringstream strs2;
    strs2 << (posTrajType+1);
    combType = "Path_"+strs2.str()+"_";

    if (testing == 0) {
        trialType = "Tr";
    }
    else {
        trialType = "Te";
    }

    std::string end(".csv");

    stringstream strs;
    std::string ch;

    strs << (rec+1);
    ch = strs.str();

    std::string tempFileAngVel = (fileAngVel+trialType+stim+combType+ch+end);

    angVel.open(tempFileAngVel.c_str());
    if (angVel.fail()){
        fileErrors << tempFileAngVel << endl;
        fileErrors << trialType << " Ang Vel " << combType << " " << (rec+1) << " failed" << std::endl;
    }

    if (posTrajType >= 0 && posTrajType < 6) { 
        int c = 0;
        float dummy;
        for (int l = 0; getline(angVel, line) && l <= trial; l++) {
            if (l == trial) {
                istringstream myline(line);
                while(getline(myline, csvItem, ',')) {
                    std::istringstream str(csvItem);
                    str >> dummy;
                    if (std::isnan(dummy) || std::isinf(dummy)) {
                        dummy = 0.0f;
                    }
                    out_arr[c] = dummy;
                    c++;
                    str.clear();
                }
            }
        }
    }
    else {
        fileErrors << "Wrong pos/traj code" << endl;
    }

    angVel.close();
    fileErrors.close();
    return out_arr;
}

vector<float> getLinVel (string dataDir, int rec, int posTrajType, bool testing, int trial, int bins) {
    vector<float> out_arr(bins, 0.0f);

    std::ofstream fileErrors;
    fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);
    std::string line, csvItem;

    std::ifstream linVel;

    std::string fileLinVel(dataDir+"linVel/linVel");
    std::string stim("Stim_");

    std::string combType;
    std::string trialType;
    stringstream strs2;
    strs2 << (posTrajType+1);
    combType = "Path_"+strs2.str()+"_";

    if (testing == 0) {
        trialType = "Tr";
    }
    else {
        trialType = "Te";
    }

    std::string end(".csv");

    stringstream strs;
    std::string ch;

    strs << (rec+1);
    ch = strs.str();

    std::string tempFileLinVel = (fileLinVel+trialType+stim+combType+ch+end);

    linVel.open(tempFileLinVel.c_str());
    if (linVel.fail()){
        fileErrors << tempFileLinVel << endl;
       fileErrors << trialType << " Lin Vel " << combType << " " << (rec+1) << " failed" << std::endl;
    }

    if (posTrajType >= 0 && posTrajType < 6) {
        int c = 0;
        float dummy;
        for (int l = 0; getline(linVel, line) && l <= trial; l++) {
            if (l == trial) {
                istringstream myline(line);
                while(getline(myline, csvItem, ',')) {
                    std::istringstream str(csvItem);
                    str >> dummy;
                    if (std::isnan(dummy) || std::isinf(dummy)) {
                        dummy = 0.0f;
                    }
                    out_arr[c] = dummy;
                    c++;
                    str.clear();
                }
            }
        }
    }
    else {
        fileErrors << "Wrong pos/traj code" << endl;
    }
    linVel.close();
    fileErrors.close();
    return out_arr;
}

vector<float> getHeadDir (string dataDir, int rec, int posTrajType, bool testing, int trial, int bins) {
    
    vector<float> out_arr(bins, 0.0f);

    std::ofstream fileErrors;
    fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);
    std::string line, csvItem;
       
    std::ifstream headDir;

    std::string fileHeadDir(dataDir+"headDir/headDir");
    std::string stim("Stim_");

    std::string combType;
    std::string trialType;
    stringstream strs2;
    strs2 << (posTrajType+1);
    combType = "Path_"+strs2.str()+"_";
    
    if (testing == 0) {
        trialType = "Tr";
    }
    else {
        trialType = "Te";
    }

    std::string end(".csv");

    stringstream strs;
    std::string ch;

    strs << (rec+1);
    ch = strs.str();

    std::string tempFileHeadDir = (fileHeadDir+trialType+stim+combType+ch+end);

    headDir.open(tempFileHeadDir.c_str());
    if (headDir.fail()){
        fileErrors << tempFileHeadDir << endl;
        fileErrors << trialType << " Head Dir " << combType << " " << (rec+1) << " failed" << std::endl;
    }

    if (posTrajType >= 0 && posTrajType < 6) {
        int c = 0;
        float dummy;
        for (int l = 0; getline(headDir, line) && l <= trial; l++) {
            if (l == trial) {
                istringstream myline(line);
                while(getline(myline, csvItem, ',')) {
                    std::istringstream str(csvItem);
                    str >> dummy;
                    if (std::isnan(dummy) || std::isinf(dummy)) {
                        dummy = 0.0f;
                    }
                    out_arr[c] = dummy;
                    c++;
                    str.clear();
                }
            }
        }
    }
    else {
        fileErrors << "Wrong pos/traj code" << endl;
    }

    headDir.close();
    fileErrors.close();
    return out_arr;
}

vector<float> getOcc (string dataDir, int rec, int posTrajType, bool testing, int trial, int bins) {
    vector<float> out_arr(bins, 0.0f);

    std::ofstream fileErrors;
    fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);
    std::string line, csvItem;

    std::ifstream occ;

    std::string fileOcc(dataDir+"occ/occ");
    std::string stim("Stim_");

    std::string combType;
    std::string trialType;
    stringstream strs2;
    strs2 << (posTrajType+1);
    combType = "Path_"+strs2.str()+"_";

    if (testing == 0) {
        trialType = "Tr";
    }
    else {
        trialType = "Te";
    }

    std::string end(".csv");

    stringstream strs;
    std::string ch;

    strs << (rec+1);
    ch = strs.str();

    std::string tempFileOcc = (fileOcc+trialType+stim+combType+ch+end);

    occ.open(tempFileOcc.c_str());
    if (occ.fail()){
        fileErrors << tempFileOcc << endl;
        fileErrors << trialType << " Occ " << combType << " " << (rec+1) << " failed" << std::endl;
    }

    if (posTrajType >= 0 && posTrajType < 6) { 
        int c = 0;
        float dummy;
        for (int l = 0; getline(occ, line) && l <= trial; l++) {
            if (l == trial) {
                istringstream myline(line);
                while(getline(myline, csvItem, ',')) {
                    std::istringstream str(csvItem);
                    str >> dummy;
                    if (std::isnan(dummy) || std::isinf(dummy)) {
                        dummy = 0.0f;
                    }
                    out_arr[c] = dummy;
                    c++;
                    str.clear();
                }
            }
        }
    }
    else {
        fileErrors << "Wrong pos/traj code" << endl;
    }

    occ.close();
    fileErrors.close();
    return out_arr;
}


// void getCoordPrefs (float * out_arr_x, float * out_arr_y, int size) {
//     std::ofstream fileErrors;
//     fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);

//     std::ifstream xpref;
//     xpref.open(dataDir+"PreferredXCoords.csv");
//     std::ifstream ypref;
//     ypref.open(dataDir+"PreferredYCoords.csv");

//     if (xpref.fail()) {
//         fileErrors << "Couldn't open preferred x coords " << endl;
//     }
//     if (ypref.fail()) {
//         fileErrors << "Couldn't open preferred y coords " << endl;
//     }

//     for (int i = 0; i < size; i++) {
//         xpref >> out_arr_x[i];
//         ypref >> out_arr_y[i]; 
//         xpref.get();
//         ypref.get();
//     }

//     xpref.close();
//     ypref.close();
// }

// void getProgPrefs (float * out_arr, int size) {
//     std::ofstream fileErrors;
//     fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);
    
//     std::ifstream progpref;
//     progpref.open(dataDir+"PreferredProg.csv");

//     if (progpref.fail()) {
//         fileErrors << "Couldn't open preferred progression " << endl;
//     }

//     for (int i = 0; i < size; i++) {
//         progpref >> out_arr[i];
//         progpref.get();
//     }

//     progpref.close();
// }

// void getLVPrefs (float * out_arr, int size) {
//     std::ofstream fileErrors;
//     fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);
    
//     std::ifstream linvelpref;
//     linvelpref.open(dataDir+"PreferredLV.csv");

//     if (linvelpref.fail()) {
//         fileErrors << "Couldn't open preferred lin vel " << endl;
//     }

//     for (int i = 0; i < size; i++) {
//         linvelpref >> out_arr[i];
//         linvelpref.get();
//     }

//     linvelpref.close();
// }

// void getAVPrefs (float * out_arr, int size) {
//     std::ofstream fileErrors;
//     fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);

//     std::ifstream angvelpref;
//     angvelpref.open(dataDir+"PreferredAV.csv");

//     if (angvelpref.fail()) {
//         fileErrors << "Couldn't open preferred ang vel " << endl;
//     }

//     for (int i = 0; i < size; i++) {
//         angvelpref >> out_arr[i];
//         angvelpref.get();
//     }

//     angvelpref.close();
// }

vector<vector<float> > loadData(string file, int numRow, int numCol) {
    vector<vector<float> > out_data(numRow, vector<float>(numCol, 0.0f));
    // out_data = new float*[numRow];
    // for (int i = 0; i < numRow; i++) {
    //     out_data[i] = new float[numCol];
    // }

    std::ofstream fileErrors;
    fileErrors.open("./data_file_errors.txt", std::ofstream::out | std::ofstream::app);

    ifstream ip;
    ip.open(file.c_str());
    
    if (ip.fail()) {
        fileErrors << "Couldn't open " << file << endl;
    }

    string line, field;
    if (!ip.eof()) {
        for (int i = 0; i < numRow; i++) {
            getline(ip, line);
            istringstream s(line);

            for (int j = 0; j < numCol; j++) {
                getline(s, field, ',');
                istringstream str(field);
                str >> out_data[i][j];
            }
        }
    }
    ip.close();
    return out_data;
}

string getFileName(string dataDir, string inputVar, int path) {
    string inputFold = (dataDir + inputVar + "/");
    string pathString = "_Path_";
    string name_suffix = ".csv";
    string filename;

    stringstream path_name_id_ss;
    path_name_id_ss << path;
    string path_name_id = path_name_id_ss.str();

    filename = (inputFold + inputVar + pathString + path_name_id + name_suffix);
    return filename;
}
