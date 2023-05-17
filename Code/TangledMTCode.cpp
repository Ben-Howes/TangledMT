#include <iostream>
#include <cmath>
#include <random>
#include<fstream>
#include<algorithm>
#include <fstream>
#include <experimental/filesystem>
#include <sys/stat.h> // mkdir
#include <dirent.h>
#include <thread>
#include <chrono>

using namespace std;

////////////////////////////
// User-Defined Variables //
////////////////////////////

// These variables can be easily changed to alter the output of the model

string dir = "../Results/TNM_Output/";   // Directory for output of model
int fragmentation = 0;                      // If 0 creates new community in non-fragmented landscape, if 1 loads in community and loads in your own landscape from command line argument 4.
                                            // but this can be changed in the command line argument, as argument 3.

int defSeed = 1;                            // This is the default seed that will be used if one is not provided in the command line argument (recommend using command line

const int cellRows = 11;                    //Sets the number of cells in the rows of the landscape (Note: must match landscape file used in code, if using file).
const int cellCols = 11;                    // Sets the number of cells in the columns of the landscape (Note: must match landscape file used in code).
const int numCells = cellRows * cellCols;   // Sets the number of cells in the landscape (this needs to be created outside of this code, likely in R, with appropriate distances etc).
int landscapeArray[cellCols][cellRows];     // Landscape array for filling with default values, or reading in from landscape created in .txt format.
int landscapeCoords[numCells][2];           // Coordinates (x,y) of each cell in the landscape for use in distance calculations.
double distArray[numCells][numCells];       // Distance from each cell to every other cell, for use in dispersal and interactions.
int Rfr = 10;                               // Set carrying capacity which will be the same for all cells.

const int numSpec = 100;                    // Number of species in the model (each species will have interacitons and mass associated with it).
const int numGenPre = 10000;                 // Number of generations to run the model pre fragmentation. Each generation is broken into time steps.
const int numGenPost = 10000;                //2000 Number of generations to run the model after fragmentation. Each generation is broken into time steps.
const int initPop = 100;                    // Number of individuals to put into each cell at the start of the model.

const float probDeath = 0.15;               // Probability of individual dying if chosen.
double probImm = 0.001;                      // Probability of an individual immigrating into a cell (individual of a random species suddenly occurring in a cell).
double probImmFrag = 0.001;                  // Probability of an individual immigrating into a cell after fragmentation (individual of a random species suddenly occurring in a cell).
float probDisp = 0.001;                       // Probability of an individual dispersing from one cell to another cell. This is a baseline, and will increase from this with increasing density.
const float probInt = 0.5;                  // Fraction of non-zero interactions. Aka approximate fraction of interactions that will occur.
double dispDist = 1;                         // Store dispersal distance

int weightInt = 5;                         // Weighting of importance of interactions. Higher value puts more importance on interactions in calculating probability of offspring.

///////////////////////
// Define Variables //
//////////////////////

// These variables should not be changed unless the model itself is being edited

int totalGen = numGenPre + numGenPost;      // Total generations the model will run for.
static double J[numSpec][numSpec];                 // J-matrix which includes interactions between all species, with number of rows and cols equal to number of species. 
static double M[numSpec];                        // Stores body mass of species
static double disp[numSpec];                       // Stores dispersal ability for each species.
int totalPop = 0;                           // Stores the total population across all cells in the model at a given generation.
int cellPop[numCells];                      // Stores the total population in each cell at a given generation.
vector <int> totalPopSpec[2];                      // Stores the total population of each species at a given generation.
int cellPopSpec[numCells][numSpec];             // Stores the total population of each species in each cell at a given generation.
int totalRich = 0;                          // Stores the total species richness of the model at a given generation.
int cellRich[numCells];                     // Stores the species richness of a specific cell at a given generation.
int cellList[numCells][2];                     // Lists numbers of cell (e,g with 6 cells it would read, 0,1,2,3,4,5), and whether they are non-forest or forest (0 = non-forest, 1 = forest).
vector <int> forestCellList;                    // Lists cells that are forested.
int cellOrder[numCells];                        // Array to store the randomly selected order of cells for dynamics.
double distMatrix[numCells][numCells];          // Stores distances between all cells in the landscape.
double probDispDen = 0;                         // Store density dependent probability of diserpsal
int seed;
int immNum = 0;                                 // Counts number of immigrations that occur
int dispNum = 0;                                // Counts number of dispersals that occur

int tmax = 50;                              // Number of times each dynamic could happen per cell per generation (e.g probDeath could happen 50 times per cell per gen).

static const double two_pi  = 2.0*3.141592653;

//////////////////////////
// Initialise functions //
//////////////////////////

void createJMatrix(double (&J)[numSpec][numSpec], double probInt, mt19937& eng);
void createM(double (&M)[numSpec], mt19937& eng);
void createDisp(double (&disp)[numSpec]);
double uniform(mt19937& eng);
double gaussian(mt19937& eng);
void initialisePop(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], 
    int &totalRich, int numSpec, int numCells, mt19937& eng);
int chooseInRange(int a, int b, mt19937& eng);
void immigration(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], int &totalRich ,
double prob, int cell, int numSpec, int &immNum, mt19937& eng);
void kill(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], int &totalRich,
double prob, int cell, int numSpec, mt19937& eng);
void reproduction(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells], 
int (&cellList)[numCells][2], double (&J)[numSpec][numSpec], int cell, int numSpec, int Rfr, mt19937& eng);
double calculateInteractions(double (&J)[numSpec][numSpec], int (&cellPopSpec)[numCells][numSpec], int (&cellPop)[numCells],
    int cell, int numSpec, int ind);
void shuffle(int arr[], int arrElements, mt19937& eng);
int randomInd(int (&cellPopSpec)[numCells][numSpec], int (&cellPop)[numCells], int numSpec, int cell, mt19937& eng);
void cellCoords(int (&landscapeCoords)[numCells][2], int cols, int rows, int numCells);
void getDistances(double (&distArray)[numCells][numCells], int landscapeCoords[][2], int cellRows, int cellCols);
void dispersal(int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], int &totalRich,
double (&distArray)[numCells][numCells], double (&disp)[numSpec], double prob, int cell, int numSpec, int &dispNum, mt19937& eng); 
vector<int> findValidCells(double (&distArray)[numCells][numCells], double distance, int cell);
double dispersalProb(int cellPop[numCells],int Rfr, int cell, double probDeath, double probDisp);
void fillCellList(int landscapeArray[][cellRows], int cellList[][2], int cellCols, int cellRows);
void getForestCellList(int cellList[][2], vector <int> &forestCellList);

void store2ColFiles(ofstream &stream, int firstCol, int secondCol);
void store3ColFiles(ofstream &stream, int firstCol, int secondCol, int thirdCol[]);
void store3ColVec(ofstream &stream, int firstCol, vector <int> secondCol, vector <int> thirdCol);
void store4ColFiles(ofstream &stream, int firstCol, int secondCol, int thirdCol, int fourthCol[][numSpec]);
void store4ColVec(ofstream &stream, int firstCol, vector <int> secondCol, int thirdCol, vector <vector<int>> fourthCol);
void storeVec(vector <int> vec, int vecElements, string fileName, string outpath);
void store2DVec(vector <vector<int>> vec, int c, int r, string fileName, string outpath);
void storeNum(int num, string fileName, string outpath);
void storeParam(string fileName, string outpath);
void readNum(int &num, string fileName, string inpath);

///////////////////////////////
//          Templates        //
///////////////////////////////

template <typename type>
void printArray(type arr[], int arrElements) {
    for (int i = 0; i < arrElements; i++){
        cout << arr[i] << " " << endl;
    }
    
}

template <typename type, int rows>
void print2DArray(type arr[][rows], int c, int r) {

    for (int i = 0; i < c; i++){
        for (int j = 0; j < r; j++){
            cout << arr[i][j] << " ";
        }
    cout << endl;
    }
}

template <typename type>
void storeArray(type arr[], int arrElements, string fileName, string outpath) {

    ofstream array;
    array.open(outpath + fileName);

    for (int i = 0; i < arrElements; i++){
        array << arr[i] << endl;
    }

    array.close();
    
}

template <typename type, int second>
void store2DArray(type arr[][second], int c, int r, string fileName, string outpath) {

    ofstream array;
    array.open(outpath + fileName);

    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            array << arr[i][j] << " ";
        }
        array << endl;
    }

    array.close();
    
}

template <typename type>
void readArray(type arr[], int arrElements, string fileName, string inpath) {

    ifstream array(inpath + fileName);

    for (int i = 0; i < arrElements; i++){
        array >> arr[i];
    }

}

template <typename type, int rows>
void read2DArray(type arr[][rows], int c, int r, string fileName, string inpath) {

    ifstream array(inpath + fileName);

    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            array >> arr[i][j];
        }
    }
}

template <                                                  //timer template FC
    class result_t   = std::chrono::seconds,
    class clock_t    = std::chrono::steady_clock,
    class duration_t = std::chrono::seconds
>
auto since(std::chrono::time_point<clock_t, duration_t> const& start)
{
    return std::chrono::duration_cast<result_t>(clock_t::now() - start);
}

///////////////////////////////
// Main Tangled Nature Model //
///////////////////////////////


int main(int argc, char *argv[]) {
    auto start = std::chrono::steady_clock::now(); // start timer FC

    string fraginpath, fragName;
    int newSeed;
       
    /////////////
    // Check command Line
    /////////////

    // We only want to run under two conditions
    // 1) With just the seed defined, which means we are making a new community
    // 2) With the old community seed defined, a new seed, fragmentation (set to 1), and a landscape file to fragment to
    // Under all other circumstances we want to exit the program (other than no arguments, when we just make a new community with seed 1 as an example)
    
    if(argc == 2) {seed = atof(argv[1]);} else {seed = defSeed;}
    
    if(argc == 7) {seed = atof(argv[1]); dispDist = atof(argv[2]); probImm = atof(argv[3]);
        probDisp = atof(argv[6]);} 
    
    if(argc == 10) {
        seed = atof(argv[1]);
        newSeed = atof(argv[2]);
        dispDist = atof(argv[3]); 
        probImmFrag = atof(argv[4]);
        probDisp = atof(argv[5]);
        fragmentation = atof(argv[6]);
        fragName = argv[7];
        fraginpath = "../../Data/Fragments/";
    }

    if(argc != 2 && argc != 7 && argc != 10) {
        cout << "Unexpected number of command line arguments" << endl;
        cout << "Either input a seed for new community" << endl;
        cout << "Or input the old seed, new seed, fragmentation (1), and a fragmented landscape file path" << endl;
        cout << "Make sure there are input for dispersal distance, interaction distance and Intrascesific competition" << endl;
        exit(0);
    }

    // Set working directories based on the seed taken from the command line argument
    string outpath = dir + "Seed_" + to_string(seed);       // Folder that outputs are saved to.
    string respath = outpath + "/Results";                  // Folder for output of results.
    string landpath = outpath + "/Landscape";                // Folder for output of landscape files.
    string finpath = outpath + "/Final_Output";              // Folder for final generation outputs for use in future models.
    string fragpath = outpath + "/Fragmentation";               // Folder to hold future fragmentation outcomes
    string finfragpath = outpath + "/Final_Frag_Output";              // Folder for final generation outputs for use in future models.

    // Make new directory at output location named after the seed

    DIR* checkDir = opendir(outpath.c_str());
    if (checkDir && fragmentation == 0) {
        closedir(checkDir);
        cout << "Directory for seed " << seed << " already exists, exiting....." << endl;
        exit(0);
    } else if(checkDir && fragmentation == 1) {
        closedir(checkDir);
        cout << "Community exists, fragmenting...." << endl;
    } else if (!checkDir && fragmentation == 0){
        mkdir(outpath.c_str(),07777);
        cout << "Directory does not exist, creating... at " << outpath << endl;
    } else if(!checkDir && fragmentation == 1) {
        cout << "Directory doesnt exist, so no community to load....exiting" << endl;
        exit(0);
    }

    // Add more directories. 
    // One to store the final outputs for use in future models, another for landscape files, and another for the full results.
    if(fragmentation == 0) {
        mkdir(respath.c_str(), 07777);
        mkdir(finpath.c_str(), 07777);
        mkdir(landpath.c_str(), 07777);
    }
    // If fragmentation is occuring make a main fragmentation folder if it doesn't already exist
    if(fragmentation == 1) mkdir(fragpath.c_str(), 07777);
    // Then make a folder within that for the specific seed used for this fragmentation
    string fragseedpath = fragpath + "/Seed_" + to_string(newSeed);
    if(fragmentation == 1) mkdir(fragseedpath.c_str(), 07777);
    // Then add a folder within the fragmentaiton folder for this specific landscape structure, which will be named after the landscape
    string fragoutpath = fragseedpath + "/" + fragName;
    if(fragmentation == 1) mkdir(fragoutpath.c_str(), 07777);

    // Fill cellList with number of each cell (0,1,2,3 etc)
    // second column in cellList currently blank as it could change depending on the landscape (1 = forest, 0 = non-forest)
    for (int i = 0; i < numCells; i++) {cellList[i][0] = i; cellOrder[i] = i;}

    if(fragmentation == 1) {seed = newSeed;}

    mt19937 eng(seed); // Seed random numbers

    // Check if we are loading a previous community to undergo fragmentation
    // or if we need to create a new community
    if(fragmentation == 1) {

        read2DArray<double, numSpec>(J, numSpec, numSpec, "/JMatrix.txt", outpath);
        readArray<double>(M, numSpec, "/M.txt", outpath);
        readArray<double>(disp, numSpec, "/disp.txt", outpath);

        // Read in predefined landscape with 1s for forest cells and 0 for non-forest cells
        read2DArray<int, cellRows>(landscapeArray, cellCols, cellRows, fragName, fraginpath);
        fillCellList(landscapeArray, cellList, cellCols, cellRows);
        store2DArray<int, 2>(cellList, 2, numCells, "/cellList.txt", fragoutpath);
        store2DArray<int, cellRows>(landscapeArray, cellCols, cellRows, "/landscape.txt", fragoutpath);

        // Get list of cells which are forest cells
        getForestCellList(cellList, forestCellList);

        // Get coordinate for each cell in the landscape for distance calculation
        cellCoords(landscapeCoords, cellCols, cellRows, numCells);
        // Calculate distances between all cells
        getDistances(distArray, landscapeCoords, cellRows, cellCols);

        // store2DArray<double, numCells>(distArray, numCells, numCells, "/distArray_frag.txt", fragpath);

        readNum(totalPop, "/final_totalPop.txt", finpath);
        readNum(totalRich, "/final_totalRich.txt", finpath);
        readArray<int>(cellPop, numCells, "/final_cellPop.txt", finpath);
        readArray<int>(cellRich, numCells, "/final_cellRich.txt", finpath);
        // readArray<int>(totalPopSpec, numSpec, "/final_totalPopSpec.txt", finpath);
        read2DArray<int, numSpec>(cellPopSpec, numSpec, numCells, "/final_cellPopSpec.txt", finpath);
        
       
    } else {

        // Create and store species traits
        createJMatrix(J, probInt, eng);
        store2DArray<double, numSpec>(J, numSpec, numSpec, "/JMatrix.txt", outpath);
        createM(M, eng);
        storeArray<double>(M, numSpec, "/M.txt", outpath);
        createDisp(disp);
        storeArray<double>(disp, numSpec, "/disp.txt", outpath);
        storeParam("/Parameters.txt", outpath);

        // Make default landscape based on rows, columns and numCells provided (aka all cells forest)
        for (int i = 0; i < cellRows; i++) {fill_n(landscapeArray[i], cellCols, 1);}
        for (int i = 0; i < numCells; i++) {cellList[i][1] = 1;}
        store2DArray<int, cellRows>(landscapeArray, cellCols, cellRows, "/landscape.txt", landpath);
        store2DArray<int, 2>(cellList, 2, numCells, "/cellList.txt", landpath);

        // Get list of cells which are forest cells
        getForestCellList(cellList, forestCellList);

        // Get coordinates for each cell for use in calculating distance
        cellCoords(landscapeCoords, cellCols, cellRows, numCells);

        // Calculate distances between all cells
        getDistances(distArray, landscapeCoords, cellRows, cellCols);
        store2DArray<double, numCells>(distArray, numCells, numCells, "/distArray.txt", landpath);

        // Set all cells to have no individuals in them
        // fill_n(totalPopSpec, numSpec, 0);
        fill_n(cellRich, numCells, 0);
        fill_n(cellPop, numCells, 0);
        for (int i = 0; i < numCells; i++) {fill_n(cellPopSpec[i], numSpec, 0);}

        // Initialise population and store the initial population to file
        initialisePop(totalPopSpec, cellPopSpec, totalPop, cellPop, cellRich, totalRich, numSpec, numCells, eng);
        store2DArray<int, numSpec>(cellPopSpec, numSpec, numCells, "/initial_Pop.txt", outpath);

    }

    // Open files to output to
    if (fragmentation == 1) {
        respath = fragoutpath + "/Results";
        finfragpath = fragoutpath + "/Final_Frag_Output";              // Folder for final generation outputs for use in future models.
        storeParam("/FragParameters.txt", fragoutpath);

        mkdir(respath.c_str(), 07777);     
        mkdir(finfragpath.c_str(), 07777);

    }

    ofstream s_totalPop;
    s_totalPop.open(respath + "/totalPop.txt");
    ofstream s_totalRich;
    s_totalRich.open(respath + "/totalRich.txt");
    ofstream s_cellRich;
    s_cellRich.open(respath + "/cellRich.txt");
    ofstream s_cellPop;
    s_cellPop.open(respath + "/cellPop.txt");
    ofstream s_totalPopSpec;
    s_totalPopSpec.open(respath + "/totalPopSpec.txt");
    ofstream s_cellPopSpec;
    s_cellPopSpec.open(respath + "/cellPopSpec.txt");
    
    // Is the number of generations for before fragmentation or after fragmentation?
    int numGen;
    if(fragmentation == 0) {numGen = numGenPre;} else {numGen = numGenPost;}

    // Start model dynamics
    for (int i = 0; i < numGen; i++) {
        // Each dynamic happens tmax times per generation per cell (tmax = 50)
        for (int t = 0; t < tmax; t++) {
            shuffle(cellOrder, numCells, eng);
            //Loop through each cell of the landscape
            for (int j = 0; j < numCells; j++) {
                int cell = cellOrder[j];
                kill(totalPopSpec, cellPopSpec, totalPop, cellPop,  cellRich, totalRich, probDeath, cell, numSpec, eng);
                reproduction(totalPopSpec, cellPopSpec, totalPop, cellPop, cellList, J, cell, numSpec, Rfr, eng);
                probDispDen = dispersalProb(cellPop, Rfr, cell, probDeath, probDisp);
                dispersal(cellPopSpec, totalPop, cellPop, cellRich, totalRich, distArray, disp, probDispDen, cell, numSpec, dispNum, eng);
                if(fragmentation == 0) {
                    immigration(totalPopSpec, cellPopSpec, totalPop, cellPop, cellRich, totalRich, probImm, cell, numSpec, immNum, eng);
                } else {
                    immigration(totalPopSpec, cellPopSpec, totalPop, cellPop, cellRich, totalRich, probImmFrag, j, numSpec, immNum, eng);
                }
            }
        }

        // Live output to console
        if(i%50 == 0) { cout << "Gen: " << i << "/" << numGen << " | Total Pop: " << totalPop << " | Total Richness: " << totalRich << "\n";}

        if ((i+1) < 50 || (i+1)%50 == 0) {

            // After each generation store all outputs
            store2ColFiles(s_totalPop, i, totalPop);
            store2ColFiles(s_totalRich, i, totalRich);
            store3ColFiles(s_cellRich, i, numCells, cellRich);
            store3ColFiles(s_cellPop, i, numCells, cellPop);
            // store3ColFiles(s_totalPopSpec, i, numSpec, totalPopSpec);
            store4ColFiles(s_cellPopSpec, i, numCells, numSpec, cellPopSpec);
        }
    }
    // Store final outputs of burn ins
    if(fragmentation == 0) {
        storeNum(totalPop, "/final_totalPop.txt", finpath);
        storeNum(totalRich, "/final_totalRich.txt", finpath);
        storeArray<int>(cellPop, numCells, "/final_cellPop.txt", finpath);
        storeArray<int>(cellRich, numCells, "/final_cellRich.txt", finpath);
        // storeArray<int>(totalPopSpec, numSpec, "/final_totalPopSpec.txt", finpath);
        store2DArray<int, numSpec>(cellPopSpec, numSpec, numCells, "/final_cellPopSpec.txt", finpath);
        storeNum(immNum, "/total_Immigrations.txt", finpath);
        storeNum(dispNum, "/total_Dispersals.txt", finpath);
        // storeVec(totalPopSpec, totalPopSpec.size(), "/totalPopSpec.txt", finpath);
    }
    // Store final outputs for fragment runs
    if(fragmentation == 1) {
        storeNum(totalPop, "/final_frag_totalPop.txt", finfragpath);
        storeNum(totalRich, "/final_frag_totalRich.txt", finfragpath);
        storeArray<int>(cellPop, numCells, "/final_frag_cellPop.txt", finfragpath);
        storeArray<int>(cellRich, numCells, "/final_frag_cellRich.txt", finfragpath);
        // storeArray<int>(totalPopSpec, numSpec, "/final_frag_totalPopSpec.txt", finfragpath);
        store2DArray<int, numSpec>(cellPopSpec, numSpec, numCells, "/final_frag_cellPopSpec.txt", finfragpath);
        storeNum(immNum, "/total_Immigrations.txt", finfragpath);
        storeNum(dispNum, "/total_Dispersals.txt", finfragpath);
    }


    // Close all streams to files
    s_totalPop.close(); s_totalRich.close(); s_cellPop.close(); s_cellRich.close(); s_totalPopSpec.close(); s_cellPopSpec.close();

    // end timer FC
    cout << "Elapsed(s)=" << since(start).count() << endl; 
    return 0;

}

///////////////////////
//     Functions    //
//////////////////////

void fillCellList(int landscapeArray[][cellRows], int cellList[][2], int cellCols, int cellRows) {

    int i = 0;

    while (i < numCells) {
        for (int j = 0; j < cellRows; j++) {
            for (int k = 0; k < cellCols; k++) {
                cellList[i][1] = landscapeArray[j][k];
                i++;
            } 
        }

    }
}

void getForestCellList(int cellList[][2], vector <int> &forestCellList) {

    for (int i = 0; i < numCells; i++) {
        if(cellList[i][1] == 1) {
            forestCellList.push_back(cellList[i][0]);
        }
    }
}

void readNum(int &num, string fileName, string inpath) {

    ifstream file (inpath + fileName);

    file >> num;

}

void storeNum(int num, string fileName, string outpath) {

    ofstream file;
    file.open(outpath + fileName);

    file << num << "\n";

    file.close();

}

void storeParam(string fileName, string outpath) {

    ofstream file;
    file.open(outpath + fileName);

    file << "Dispersal_Distance " << dispDist << "\n";
    file << "Immigration_RateB " << probImm << "\n";
    file << "Immigration_RateF " << probImmFrag << "\n";
    file << "Dispersal_Probability " << probDisp << "\n";

    file.close();

}

void store4ColFiles(ofstream &stream, int firstCol, int secondCol, int thirdCol, int fourthCol[][numSpec]) {

    double H;
    for (int j = 0; j < secondCol; j++) {
        for (int k = 0; k < thirdCol; k++){
            if(fourthCol[j][k] > 0) {
                H = calculateInteractions(J, cellPopSpec, cellPop, j, numSpec, k);
                stream << firstCol+1 << " " << j+1 << " " << k+1 << " " << fourthCol[j][k] << " " << H << "\n";
            }
        }
    }
}

void store3ColVec(ofstream &stream, int firstCol, vector <int> secondCol, vector <int> thirdCol) {

    int i = 0;
    for (int j: secondCol) {
        stream << firstCol+1 << " " << j+1 << " " << thirdCol[i] << "\n";
        i++;
    }
}

void store4ColVec(ofstream &stream, int firstCol, vector <int> secondCol, int thirdCol, vector <vector<int>> fourthCol) {

    int i = 0; double H;
    for (int j: secondCol) {
        for (int k = 0; k < thirdCol; k++) {
            if(fourthCol[i][k] > 0) {
                H = calculateInteractions(J, cellPopSpec, cellPop, j, numSpec, k);
                stream << firstCol+1 << " " << j+1 << " " << k+1 << " " << fourthCol[i][k] << " " << H << "\n";
            }
        }
        i++;
    }
}


void store3ColFiles(ofstream &stream, int firstCol, int secondCol, int thirdCol[]) {

    for (int j = 0; j < secondCol; j++) {
        stream << firstCol+1 << " " << j+1 << " " << thirdCol[j] << "\n";
    }
}

void store2ColFiles(ofstream &stream, int firstCol, int secondCol) {

    stream << firstCol+1 << " " << secondCol << "\n";

}

void storeVec(vector <int> vec, int vecElements, string fileName, string outpath) {

    ofstream out;
    out.open(outpath + fileName);

    for (int i = 0; i < vecElements; i++){
        out << vec[i] << endl;
    }

    out.close();
    
}

void store2DVec(vector <vector<int>> vec, int c, int r, string fileName, string outpath) {

    ofstream out;
    out.open(outpath + fileName);

    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            out << vec[i][j] << " ";
        }
        out << endl;
    }

    out.close();
    
}

double dispersalProb(int cellPop[numCells],int Rfr, int cell, double probDeath, double probDisp) {

    double Nequ, pDispEff, cutOff, grad;
    pDispEff = probDisp;

    cutOff = Rfr*(4+log10((1-probDeath)/probDeath));        // Cut off for density dependent linear, max is 47.5
    // Nequ = Rfr*(3+log((1-probDeath)/probDeath));			//calculate density dependent Pmig if N/Nequ>.75

    // Make density dependence linear with a cut-off
    
    grad = 1/cutOff;
    if(cellPop[cell] > cutOff) {pDispEff = probDisp;} else {pDispEff = (grad*cellPop[cell])*probDisp;}
    
    // if (cellPop[cell]/Nequ > 1.5) {
    //     pDispEff = 1/(1+(1-probDisp)*pow((1+probDisp),3)/pow(probDisp,4)*pow((probDisp/(1+probDisp)),(4*(double)cellPop[cell]/Nequ)));
    // }

    return pDispEff;

}

vector <int> findValidCells(double (&distArray)[numCells][numCells], double distance, int cell) {

    vector<int>validCells;

    // Go through row of distArray for the cell and identify all cells within distance
    for (int i = 0; i < numCells; i++) {
        if(distArray[cell][i] <= distance) {
            validCells.push_back(i);
        }
    }
    return validCells;
}

void dispersal(int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], int &totalRich,
double (&distArray)[numCells][numCells], double (&disp)[numSpec], double prob, int cell, int numSpec, int &dispNum, mt19937& eng) {

    if(cellPop[cell] > 0) { // Check cell isn't empty

        int chosenInd, dispCell;
        double dispAbility;
        vector<int> validCells;

        if (uniform(eng) <= prob) {
            chosenInd = randomInd(cellPopSpec, cellPop, numSpec, cell, eng);
            dispAbility = disp[chosenInd];

            if(dispAbility >= 1) {  // If dispAbility less than 1 then the individual can't disperse
                validCells = findValidCells(distArray, dispAbility, cell);
                int arr[validCells.size()]; 
                copy(validCells.begin(), validCells.end(), arr); // Copy into array for use in shuffle function
                shuffle(arr, validCells.size(), eng);                 // Shuffle valid cells into random order
                if(arr[0] != cell) {
                    dispCell = arr[0];                             // Choose first cell in shuffled order which isn't the original cell
                } else {
                    dispCell = arr[1];
                }

                // Move individual from cell to chosenCell
                if(cellPopSpec[cell][chosenInd] == 1) {cellRich[cell]--;}
                if(cellPopSpec[dispCell][chosenInd] == 0) {cellRich[dispCell]++;}
                cellPopSpec[cell][chosenInd]--;
                cellPopSpec[dispCell][chosenInd]++;
                cellPop[cell]--;
                cellPop[dispCell]++;

                // Count dispersal
                dispNum++;                                    
            }
        }

    } 
}

void getDistances(double (&distArray)[numCells][numCells], int landscapeCoords[][2], int cellRows, int cellCols) {

double x_diff, y_diff;

    for (int i = 0; i < numCells; i++) {
        for (int j = 0; j < numCells; j++) {
            x_diff = abs(landscapeCoords[i][0] - landscapeCoords[j][0]);
            if(x_diff > cellRows/2) {x_diff = cellRows - x_diff;}
            y_diff = abs(landscapeCoords[i][1] - landscapeCoords[j][1]);
            if(y_diff > cellCols/2) {y_diff = cellCols - y_diff;}
            distArray[i][j] = sqrt((x_diff*x_diff) + (y_diff*y_diff));
        }
    }
    
}

void cellCoords(int (&landscapeCoords)[numCells][2], int cols, int rows, int numCells) {

    int i = 0;

    while(i < numCells) {
        for (int j = 0; j < rows; j++) {
            for (int k = 0; k < cols; k++) {
                landscapeCoords[i][0] = j;
                landscapeCoords[i][1] = k;
                i++;
            }
        }
    }
}

double calculateInteractions(double (&J)[numSpec][numSpec], int (&cellPopSpec)[numCells][numSpec], int (&cellPop)[numCells], int cell, int numSpec, int ind) {

    double H;
    H = 0; 

    for (int i = 0; i < numSpec; i++) {H += J[ind][i]*cellPopSpec[cell][i];}

    return H;

}

void reproduction(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells], 
int (&cellList)[numCells][2], double (&J)[numSpec][numSpec], int cell, int numSpec, int Rfr, mt19937& eng) {

    if(cellPop[cell] > 0) {

        int chosenInd;
        double H, pOff;

        chosenInd = randomInd(cellPopSpec, cellPop, numSpec, cell, eng);
        H = calculateInteractions(J, cellPopSpec, cellPop, cell, numSpec, chosenInd);
        H = (H*weightInt/cellPop[cell]) - (cellPop[cell]/Rfr); // 10 chosen as arbitrary carrying capactiy
        pOff = exp(H) / (1 + exp(H));

        if (uniform(eng) <= pOff) {
            for (int i = 0; i < totalPopSpec[0].size(); i++) {if(totalPopSpec[0][i] == chosenInd) {totalPopSpec[1][i] += 1; break;}}
            cellPopSpec[cell][chosenInd]++;
            totalPop++;
            cellPop[cell]++;
        }
    }
}

void kill(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], int &totalRich,
double prob, int cell, int numSpec, mt19937& eng) {
    
    if(cellPop[cell] > 0){ 

    int chosenInd;

        if (uniform(eng) <= prob) {
            chosenInd = randomInd(cellPopSpec, cellPop, numSpec, cell, eng);
            for (int i = 0; i < totalPopSpec[0].size(); i++){
                if(totalPopSpec[0][i] == chosenInd) {
                    if(totalPopSpec[1][i] == 1) {
                            totalPopSpec[0].erase(totalPopSpec[0].begin() + i);
                            totalPopSpec[1].erase(totalPopSpec[1].begin() + i);
                            totalRich--; 
                        } else {
                        totalPopSpec[1][i] -= 1;
                        break;
                    }
                }
            }
            
            if(cellPopSpec[cell][chosenInd] == 1) {cellRich[cell]--;}
            if(cellPopSpec[cell][chosenInd] > 0) {cellPopSpec[cell][chosenInd]--;}
            if(totalPop > 0) {totalPop--;}
            if(cellPop[cell] > 0) {cellPop[cell]--;}
        }
    }
}

// Choose a random individual (proportional to how many of that species there are)
int randomInd(int (&cellPopSpec)[numCells][numSpec], int (&cellPop)[numCells], int numSpec, int cell, mt19937& eng) {
    
    double sum, threshold;
    sum = 0;
    // Create threshold value by multiplying the total number of individuals in the cell by a number between 0-1
    threshold = uniform(eng)*cellPop[cell];

    // Get an individual of a random species (proportional to how many of that species there are) by summing the number of individuals of each species in the cell
    // and then choosing the species whose population get the sum to the threshold value
    for (int i = 0; i < numSpec; i++) {
        sum += cellPopSpec[cell][i];
        if (sum > threshold) {return i;}
    }
    cout << "Threshold for randomInd not hit, this is likely causing huge errors in results" << "\n";
    return 0;
}

void shuffle(int arr[], int arrElements, mt19937& eng) {
    
    int randomNum;

    for (int i = 0; i < arrElements; i++){
        randomNum = chooseInRange(0, arrElements-1, eng);
        swap(arr[i], arr[randomNum]);
    }
}

void immigration(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], int &totalRich ,
double prob, int cell, int numSpec, int &immNum, mt19937& eng) {
    
    int chosenSpec;
    bool exists = false;

    if(uniform(eng) <= prob) {  
        chosenSpec = chooseInRange(0, numSpec-1, eng);
        for (int i = 0; i < totalPopSpec[0].size(); i++) {
            if(totalPopSpec[0][i] == chosenSpec) {
                totalPopSpec[1][i] += 1;
                exists = true;
            }
        }
        if(exists = false) {
            totalRich++;
            totalPopSpec[0].push_back(chosenSpec);
            totalPopSpec[1].push_back(1);
        }
        
        if(cellPopSpec[cell][chosenSpec] == 0) {cellRich[cell]++;}
        cellPopSpec[cell][chosenSpec]++;
        totalPop++;
        cellPop[cell]++;
        immNum++;
    }
}

int chooseInRange(int a, int b, mt19937& eng) {
    uniform_int_distribution<> choose(a, b); // define the range [a,b], extremes included.
    return choose(eng);
}

void initialisePop(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], 
    int &totalRich, int numSpec, int numCells, mt19937& eng) {

    int chosenSpec;
    bool exists = false;

    for (int i = 0; i < numCells; i++){
        for (int j = 0; j < initPop; j++) {
            exists = false;
            chosenSpec = chooseInRange(0, numSpec-1, eng); 
            for (int i = 0; i < totalPopSpec[0].size(); i++) {
                if(totalPopSpec[0][i] == chosenSpec) {
                    totalPopSpec[1][i] += 1;
                    exists = true;}
            }
            if(exists == false) {
                totalPopSpec[0].push_back(chosenSpec);
                totalPopSpec[1].push_back(1);
                totalRich++;
            }

            if(cellPopSpec[i][chosenSpec] == 0) {cellRich[i]++;}
            cellPopSpec[i][chosenSpec]++;
            totalPop++;
            cellPop[i]++;
        }
    }
}

double uniform(mt19937& eng) {
    uniform_real_distribution<> rand(0, 1); // define the range [a,b], extremes included.
    return rand(eng);
}

double gaussian(mt19937& eng) {
    double u1, u2, normalNum;

    u1 = uniform(eng); // Create two random numbers to normally distribute
    u2 = uniform(eng);

    normalNum = 0.25*(sqrt(-2.0 * log(u1)) * cos(two_pi * u2)); // Create normally distributed vales with mean 0 and sd 0.25
    return normalNum;
}

void createJMatrix(double (&J)[numSpec][numSpec], double probInt, mt19937& eng) {
    for (int i = 0; i < numSpec; i++){
        for (int j = 0; j <= i; j++) {
            if(i == j ) {
                J[i][j] = 0;
            } else if (uniform(eng) <= probInt){
                J[i][j] = 0;
                J[j][i] = 0;
            }
            else {
                J[i][j] = gaussian(eng);
                J[j][i] = gaussian(eng);
            } 
        }
    } 
}

void createM(double (&M)[numSpec], mt19937& eng) {
    for (int i = 0; i < numSpec; i++) {
        M[i] = uniform(eng);
    }
}

void createDisp(double (&disp)[numSpec]) {
    for (int i = 0; i < numSpec; i++) {
        disp[i] = dispDist;
    }
}