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

const int cellRows = 2;                    //Sets the number of cells in the rows of the landscape (Note: must match landscape file used in code, if using file).
const int cellCols = 2;                    // Sets the number of cells in the columns of the landscape (Note: must match landscape file used in code).
int Rfr = 10;                               // Set carrying capacity which will be the same for all cells.

const int L = 8;                           // Length of binary identifiers to use in the model (genome sequences)
const int numSpec = 256;                    // Number of species in the model, the number of species must equal 2^L .
const int t = 20000;                         // Number of time steps in the model
const int initPop = 50;                    // Number of individuals to put into each cell at the start of the model.

const float probDeath = 0.15;               // Probability of individual dying if chosen.
double probImm = 0.001;                      // Probability of an individual immigrating into a cell (individual of a random species suddenly occurring in a cell).
double probImmFrag = 0.001;                  // Probability of an individual immigrating into a cell after fragmentation (individual of a random species suddenly occurring in a cell).
float probDisp = 0.001;                       // Probability of an individual dispersing from one cell to another cell. This is a baseline, and will increase from this with increasing density.
double dispDist = 1;                         // Store dispersal distance
double probMut = 0.0001;                      // Probability of a number in the genome sequence switching from 0 -> 1, or 1 -> 0

// Metabolic theory variables
double ppProb = 0.25;                       // Sets proportion of species that are primary producers
double T = 20;                             // Set temperature in kelvin (273.15 kelvin = 0 celsius)
double k = 8.6173*(10^-5);                  // Boltzmann constant

///////////////////////
// Define Variables //
//////////////////////

// These variables should not be changed unless the model itself is being edited

const int numCells = cellRows * cellCols;   // Sets the number of cells in the landscape (this needs to be created outside of this code, likely in R, with appropriate distances etc).
int landscapeArray[cellCols][cellRows];     // Landscape array for filling with default values, or reading in from landscape created in .txt format.
int landscapeCoords[numCells][2];           // Coordinates (x,y) of each cell in the landscape for use in distance calculations.
double distArray[numCells][numCells];       // Distance from each cell to every other cell, for use in dispersal and interactions.

static double J[numSpec][numSpec];                 // J-matrix which includes interactions between all species, with number of rows and cols equal to number of species. 
static double traits[numSpec][2];                        // Stores base traits of species, mass, primary producer etc
static double disp[numSpec];                       // Stores dispersal ability for each species.
int totalPop = 0;                           // Stores the total population across all cells in the model at a given generation.
vector <int> cellPop[2];                    // Stores the total population in each cell at a given generation.
vector <int> totalPopSpec[2];                      // Stores the total population of each species at a given generation.
vector <int> cellPopSpec[numCells][4];             // Stores which individuals are in which cells (1st array = species, 2nd array = mass, 3rd array = primary producer, 4th = dispersal)
vector<vector<double>> searchRates[numCells];                // Stores search rate of individual i with individual j
vector<vector<double>> attackProbs[numCells];                // Stores attack probability of individual i with individual j
vector<vector<double>> handlingTimes[numCells];                // Stores handling time of individual i with individual j
vector<vector<double>> profitabilities[numCells];                // Stores profitability of individual i with individual j
int totalRich = 0;                          // Stores the total species richness of the model at a given generation.
int cellList[numCells][2];                     // Lists numbers of cell (e,g with 6 cells it would read, 0,1,2,3,4,5), and whether they are non-forest or forest (0 = non-forest, 1 = forest).
vector <int> forestCellList;                    // Lists cells that are forested.
int cellOrder[numCells];                        // Array to store the randomly selected order of cells for dynamics.
double distMatrix[numCells][numCells];          // Stores distances between all cells in the landscape.
double probDispDen = 0;                         // Store density dependent probability of diserpsal
int seed;
int immNum = 0;                                 // Counts number of immigrations that occur
int dispNum = 0;                                // Counts number of dispersals that occur

static const double two_pi  = 2.0*3.141592653;

// Metabolic theory variables
static double T0 = 273.15;                        // 0 celsius in Kelvin - used to add to the temperature we set in C to turn into Kelvin

//////////////////////////
// Initialise functions //
//////////////////////////

void createJMatrix(double (&J)[numSpec][numSpec], double probInt, mt19937& eng);
void createTraits(double (&traits)[numSpec][2], mt19937& eng, double ppProb);
void createDisp(double (&disp)[numSpec]);
double uniform(mt19937& eng);
double gaussian(mt19937& eng);
void initialisePop(vector <int> (&cellPopSpec)[numCells][4], double (&traits)[numSpec][2], int numSpec, int numCells, mt19937& eng);
int chooseInRange(int a, int b, mt19937& eng);
void immigration(vector <int> (&cellPopSpec)[numCells][4], double prob, int cell, int numSpec, int &immNum, mt19937& eng);
void kill(vector <int> (&cellPopSpec)[numCells][4], double prob, int cell, int numSpec, mt19937& eng);
void reproduction(vector <int> (&cellPopSpec)[numCells][4], 
int (&cellList)[numCells][2], double (&traits)[numSpec][2], int cell, int numSpec, int Rfr, mt19937& eng);
double calculateInteractions(vector <int> (&cellPopSpec)[numCells][4], double (&traits)[numSpec][2], int cell, int numSpec, int ind);
void shuffle(int arr[], int arrElements, mt19937& eng);
int randomIndex(vector <int> (&cellPopSpec)[numCells][4], int pop, int numSpec, int cell, mt19937& eng);
void cellCoords(int (&landscapeCoords)[numCells][2], int cols, int rows, int numCells);
void getDistances(double (&distArray)[numCells][numCells], int landscapeCoords[][2], int cellRows, int cellCols);
void dispersal(vector <int> (&cellPopSpec)[numCells][4], double (&distArray)[numCells][numCells], double (&disp)[numSpec], double prob, int cell, int numSpec, int &dispNum, mt19937& eng); 
vector<int> findValidCells(double (&distArray)[numCells][numCells], double distance, int cell);
double dispersalProb(vector <int> (&cellPopSpec)[numCells][4], int Rfr, int cell, double probDeath, double probDisp);
void fillCellList(int landscapeArray[][cellRows], int cellList[][2], int cellCols, int cellRows);
void getForestCellList(int cellList[][2], vector <int> &forestCellList);
void removeSpecies(vector <int> (&cellPopSpec)[numCells][4], int cell, int chosenInd);
 void addSpecies(vector <int> (&cellPopSpec)[numCells][4], int cell, int chosenInd, double (&traits)[numSpec][2]);

void store2ColFiles(ofstream &stream, int firstCol, int secondCol);
void storeVec(ofstream &stream, vector <int> vec[], int gen, int cols);
void storeVecEnd(vector <int> vec[], int cols, string fileName, string outpath);
void storeNum(int num, string fileName, string outpath);
void storeParam(string fileName, string outpath);
void storeCellPopSpec(ofstream &stream, vector <int> (&cellPopSpec)[numCells][4], int gen, double (&traits)[numSpec][2]);

void calculateTotalPopSpec(vector <int> (&cellPopSpec)[numCells][4], vector <int> (&totalPopSpec)[2]);
void calculateCellPop(vector <int> (&cellPopSpec)[numCells][4], vector <int> (&cellPop)[2]);
int calcCellPop(vector <int> (&cellPopSpec)[numCells][4], int cell);
int getPop(int ind, int cell, vector <int> (&cellPopSpec)[numCells][4]);

// Mutations
int mutation(vector <int> (&cellPopSpec)[numCells][4], double prob, int chosenInd, mt19937& eng);
int ConvertToDec(int (&arr)[L]);
void ConvertToBinary(int n, int j, int (&vec)[L]);
void Bin_recursive(int n, int j, int (&vec)[L]);

// Metabolic Theory Functions
double searchRate(int Si, int Sj, double E, double (&traits)[numSpec][2]);
double attackProb(int Si, int Sj, double (&traits)[numSpec][2]);
double handlingTime(int Si, int Sj, double E, double (&traits)[numSpec][2]);
double consumptionRate(int Si, int Sj, double E, double (&traits)[numSpec][2], int Nj);
double getCellMass(int cell, vector <int> cellPopSpec[numCells][4]);
double arrhenius(double E);
void storeConsumptionRate(ofstream &stream, vector <int> (&cellPopSpec)[numCells][4], int gen, double (&traits)[numSpec][2]);
double profitability(int Si, int Sj, double (&traits)[numSpec][2]);
void optimalForaging(int Si, int cell, double (&traits)[numSpec][2], vector <int> (&cellPopSpec)[numCells][4], vector <int> (&optimalSet));
void findOptimalSet(int Si, int cell, vector <double> (&sortedProfits)[2], double (&traits)[numSpec][2], vector <int> (&optimalSet));
void sortedCols(vector <double> (&profits)[2], vector <double> (&sortedProfits)[2]);

///////////////////////////////
//          Templates        //
///////////////////////////////

template <typename type>
void printArray(type arr[], int arrElements) {
    for (int i = 0; i < arrElements; i++){
        std::cout << arr[i] << " " << endl;
    }
    
}

template <typename type, int rows>
void print2DArray(type arr[][rows], int c, int r) {

    for (int i = 0; i < c; i++){
        for (int j = 0; j < r; j++){
            std::cout << arr[i][j] << " ";
        }
    std::cout << endl;
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
        std::cout << "Unexpected number of command line arguments" << endl;
        std::cout << "Either input a seed for new community" << endl;
        std::cout << "Or input the old seed, new seed, fragmentation (1), and a fragmented landscape file path" << endl;
        std::cout << "Make sure there are input for dispersal distance, interaction distance and Intrascesific competition" << endl;
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
        std::cout << "Directory for seed " << seed << " already exists, exiting....." << endl;
        exit(0);
    } else if(checkDir && fragmentation == 1) {
        closedir(checkDir);
        std::cout << "Community exists, fragmenting...." << endl;
    } else if (!checkDir && fragmentation == 0){
        mkdir(outpath.c_str(),07777);
        std::cout << "Directory does not exist, creating... at " << outpath << endl;
    } else if(!checkDir && fragmentation == 1) {
        std::cout << "Directory doesnt exist, so no community to load....exiting" << endl;
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
        read2DArray<double, 2>(traits, 2, numSpec, "/traits.txt", outpath);
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

        store2DArray<double, numCells>(distArray, numCells, numCells, "/distArray_frag.txt", fragpath);
        // read2DArray<int, numSpec>(cellPopSpec, numSpec, numCells, "/final_cellPopSpec.txt", finpath);
        
       
    } else {

        // Create and store species traits
        // createJMatrix(J, probInt, eng);
        // store2DArray<double, numSpec>(J, numSpec, numSpec, "/JMatrix.txt", outpath);
        createTraits(traits, eng, ppProb);
        store2DArray<double, 2>(traits, 2, numSpec, "/traits.txt", outpath);
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

        // Initialise population
        initialisePop(cellPopSpec, traits, numSpec, numCells, eng);
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
    ofstream s_cellPop;
    s_cellPop.open(respath + "/cellPop.txt");
    ofstream s_totalPopSpec;
    s_totalPopSpec.open(respath + "/totalPopSpec.txt");
    ofstream s_cellPopSpec;
    s_cellPopSpec.open(respath + "/cellPopSpec.txt");
    ofstream s_consumptionRate;
    s_consumptionRate.open(respath + "/consumptionRate.txt");
    
    
    // Start model dynamics
    for (int i = 0; i < t; i++) {
        shuffle(cellOrder, numCells, eng);
        //Loop through each cell of the landscape
        for (int j = 0; j < numCells; j++) {
            int cell = cellOrder[j];
            kill(cellPopSpec, probDeath, cell, numSpec, eng);
            // reproduction(cellPopSpec, cellList, traits, cell, numSpec, Rfr, eng);
            // // probDispDen = dispersalProb(cellPopSpec, Rfr, cell, probDeath, probDisp);
            dispersal(cellPopSpec, distArray, disp, probDisp, cell, numSpec, dispNum, eng); // Now set to a constant probability
            if(fragmentation == 0) {
                immigration(cellPopSpec, probImm, cell, numSpec, immNum, eng);
            } else {
                immigration(cellPopSpec, probImmFrag, j, numSpec, immNum, eng);
            }
        }

        // Calculate richness and abundance metrics
        // Dependent on how often you want to save them
        if((i+1)%1000 == 0) {
            calculateTotalPopSpec(cellPopSpec, totalPopSpec);
            totalPop = 0;
            for (int i = 0; i < totalPopSpec[0].size(); i++){totalPop += totalPopSpec[1][i];}
            totalRich = totalPopSpec[0].size();
            calculateCellPop(cellPopSpec, cellPop);

            //After each generation store all outputs
            store2ColFiles(s_totalPop, i, totalPop);
            store2ColFiles(s_totalRich, i, totalRich);
            storeVec(s_totalPopSpec, totalPopSpec, i, 2);
            storeVec(s_cellPop, cellPop, i, 2);
            storeCellPopSpec(s_cellPopSpec, cellPopSpec, i, traits);
            // storeConsumptionRate(s_consumptionRate,  cellPopSpec, i, traits);

            //Live output to console
            std::cout << "Time Step: " << i + 1 << "/" << t << " | Total Pop: " << totalPop << " | Total Richness: " << totalRich << "\n";

        }
    }
    // // Store final outputs of burn ins
    if(fragmentation == 0) {
        storeNum(immNum, "/total_Immigrations.txt", finpath);
        storeNum(dispNum, "/total_Dispersals.txt", finpath);
        // storeVecEnd(cellPopSpec, 3, "/final_cellPopSpec.txt", finpath);
    }
    // Store final outputs for fragment runs
    if(fragmentation == 1) {
    //   save cellpopspec
        storeNum(immNum, "/total_Immigrations.txt", finfragpath);
        storeNum(dispNum, "/total_Dispersals.txt", finfragpath);
        //  storeVecEnd(cellPopSpec, 3, "/final_cellPopSpec.txt", finfragpath);
    }
    // Close all streams to files
    s_totalPop.close(); s_totalRich.close(); s_cellPop.close(); s_totalPopSpec.close(); s_cellPopSpec.close(), s_consumptionRate.close();

    // end timer FC
    std::cout << "Elapsed(s)=" << since(start).count() << endl; 
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

void store2ColFiles(ofstream &stream, int firstCol, int secondCol) {

    stream << firstCol+1 << " " << secondCol << "\n";

}

void storeVec(ofstream &stream, vector <int> vec[], int gen, int cols) {

    for (int i = 0; i < vec[0].size(); i++) {
        if(cols == 1) {
            stream << gen+1 << " " << vec[0][i] << "\n";
        } else if( cols == 2) {
            stream << gen+1 << " " << vec[0][i] + 1 << " " << vec[1][i] << "\n";
        } else if (cols == 3) {
            stream << gen+1 << " " << vec[0][i] + 1 << " " << vec[1][i] + 1 << " " << vec[2][i] << "\n";
        }
    }
}

void storeVecEnd(vector <int> vec[], int cols, string fileName, string outpath) {

    ofstream out;
    out.open(outpath + fileName);

    for (int i = 0; i < vec[0].size(); i++) {
        if(cols == 1) { 
            out << vec[0][i] << endl;
        } else if (cols == 2) {
            out << vec[0][i] << " " << vec[1][i] << endl;
        } else if (cols == 3) {
            out << vec[0][i] << " " << vec[1][i] << " " << vec[2][i] << endl;
        }
    }

    out.close();
}

double dispersalProb(vector <int> (&cellPopSpec)[numCells][4],int Rfr, int cell, double probDeath, double probDisp) {

    double Nequ, pDispEff, cutOff, grad;
    int pop = cellPopSpec[cell][0].size();
    pDispEff = probDisp;

    cutOff = Rfr*(3+log10((1-probDeath)/probDeath));        // Cut off for density dependent linear, max is 47.5

    // Make density dependence linear with a cut-off
    
    grad = 1/cutOff;
    if(pop > cutOff) {pDispEff = probDisp;} else {pDispEff = (grad*pop)*probDisp;}

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

void dispersal(vector <int> (&cellPopSpec)[numCells][4], double (&distArray)[numCells][numCells], 
double (&disp)[numSpec], double prob, int cell, int numSpec, int &dispNum, mt19937& eng) {

    int pop = cellPopSpec[cell][0].size();

    if(pop > 0) { // Check cell isn't empty
        int chosenInd, chosenIndex, dispCell;
        double dispAbility; // Store dispersal ability of chosen species
        vector<int> validCells; // store cells that individual could disperse to

        if (uniform(eng) <= prob) {
            // Get index of random individual
            chosenIndex = randomIndex(cellPopSpec, pop, numSpec, cell, eng);
            dispAbility = cellPopSpec[cell][3][chosenIndex];
            chosenInd = cellPopSpec[cell][0][chosenIndex];

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

                // Remove individual from cell it was in
                removeSpecies(cellPopSpec, cell, chosenIndex);
                // Add individual to dispCell
                addSpecies(cellPopSpec, dispCell, chosenInd, traits);
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

double calculateInteractions(vector <int> (&cellPopSpec)[numCells][4], double (&traits)[numSpec][2], int cell, int numSpec, int ind) {

    double H = 0;
    double CE = 0.5; // Conversation efficiency
    double N; // Store abundance of focal species
    vector<int> optimalSet; // Store optimal set of resources
    int Sj; // Store identities of non-focal species
    int Nj; // Store abundance of non-focal species in the cell

    // Only run the below if there is more than one species in the cell
    // as species can't interact (feed) on themselves
    if(cellPopSpec[cell][0].size() > 1) {

        // Find the set of species which are in the optimal set for focal species Si (ind)
        optimalForaging(ind, cell, traits, cellPopSpec, optimalSet);

        // Iterate over set finding consumption rate and summing for
        // total consumption rate for focal species Si
        for (int i = 0; i < optimalSet.size(); i++) {
            int Sj = optimalSet[i];
            // Per capita rate of consumption kg per ind per second x conversion efficiency
            H += CE*consumptionRate(ind, Sj, 0, traits, getPop(Sj, cell, cellPopSpec));
        }

        // Then need to find optimal set of resources for all other species in the cell
        // if our focal species is in that set then the consumption rate of consumer Sj on focal species Si
        // will be calculated and used as the loss term
        for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
            // Don't calculate optimal set for our focal species again
            if(ind != cellPopSpec[cell][0][i]) {
                Sj = cellPopSpec[cell][0][i];
                Nj = cellPopSpec[cell][1][i];
                optimalForaging(Sj, cell, traits, cellPopSpec, optimalSet);
                for (int j = 0; j < optimalSet.size(); j++) {
                    // Only include interactions with our focal species
                    if(optimalSet[j] == ind) {
                        H -= Nj*consumptionRate(Sj, ind, 0, traits, 1); 
                    }              
                }
            }
        }
    }

    return H;

}

void reproduction(vector <int> (&cellPopSpec)[numCells][4], 
int (&cellList)[numCells][2], double (&traits)[numSpec][2], int cell, int numSpec, int Rfr, mt19937& eng) {

    int pop = cellPopSpec[cell][0].size();
    double cellMass; double B0 = 4.15*std::pow(10, -8);

    if(pop > 0) {

        int chosenInd, chosenIndex;
        double H, pOff;
        
        cellMass = getCellMass(cell, cellPopSpec); // Get the total mass of all individuals in the cell
        chosenIndex = randomIndex(cellPopSpec, pop, numSpec, cell, eng);
        chosenInd = cellPopSpec[cell][0][chosenInd];
        H = calculateInteractions(cellPopSpec, traits, cell, numSpec, chosenIndex) - (cellMass/Rfr) - (B0*std::pow(cellPopSpec[cell][2][chosenIndex], 0.75)); // Rfr 10 chosen as arbitrary carrying capactiy
        pOff = exp(H) / (1 + exp(H));

        if (uniform(eng) <= pOff) {
            int mutSpec = mutation(cellPopSpec, probMut, chosenInd, eng);
            addSpecies(cellPopSpec, cell, mutSpec, traits);
        }
    }
}

void kill(vector <int> (&cellPopSpec)[numCells][4], double prob, int cell, int numSpec, mt19937& eng) {
    
    int pop = cellPopSpec[cell][0].size();

    if(pop > 0) { 

    int chosenInd;

        if (uniform(eng) <= prob) {
            chosenInd = randomIndex(cellPopSpec, pop, numSpec, cell, eng);
            removeSpecies(cellPopSpec, cell, chosenInd);
        }
    }
}

// Choose a random individual (proportional to how many of that species there are)
// NOTE: THIS NOW RETURNS THE INDEX OF THE INDIVIDUAL IN CELLPOPSPEC
int randomIndex(vector <int> (&cellPopSpec)[numCells][4], int pop, int numSpec, int cell, mt19937& eng) {
    
    double sum, threshold;
    sum = 0;
    // Create threshold value by multiplying the total number of individuals in the cell by a number between 0-1
    threshold = uniform(eng)*pop;

    // Get a random individual by going the threshold amount through the cellpopspec vector for that cell
    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        sum += 1;
        if (sum > threshold) {return i;}
    }
    
    std::cout << "Threshold for randomInd not hit, this is likely causing huge errors in results" << "\n";
    return 0;
}

void shuffle(int arr[], int arrElements, mt19937& eng) {
    
    int randomNum;

    for (int i = 0; i < arrElements; i++){
        randomNum = chooseInRange(0, arrElements-1, eng);
        swap(arr[i], arr[randomNum]);
    }
}

void immigration(vector <int> (&cellPopSpec)[numCells][4], double prob, int cell, int numSpec, int &immNum, mt19937& eng) {
    
    int chosenSpec;
    bool exists = false;

    if(uniform(eng) <= prob) {  
        chosenSpec = chooseInRange(0, numSpec-1, eng);
        addSpecies(cellPopSpec, cell, chosenSpec, traits);
        immNum++;
    }
}

int chooseInRange(int a, int b, mt19937& eng) {
    uniform_int_distribution<> choose(a, b); // define the range [a,b], extremes included.
    return choose(eng);
}

void initialisePop(vector <int> (&cellPopSpec)[numCells][4], double (&traits)[numSpec][2], int numSpec, int numCells, mt19937& eng) {

    int chosenSpec;

    for (int i = 0; i < numCells; i++){
        for (int j = 0; j < initPop; j++) {
            chosenSpec = chooseInRange(0, numSpec-1, eng); 

            addSpecies(cellPopSpec, i, chosenSpec, traits);
            // updateSeachRate
            // updateAttackProb
            // updateHandling
            // updateProfitability
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
                double oppoInt = gaussian(eng);
                if(J[i][j] > 0 & oppoInt < 0) {
                    J[j][i] = oppoInt;
                } else if (J[i][j] > 0 & oppoInt > 0) {
                    J[j][i] = -oppoInt;
                } else if (J[i][j] < 0 & oppoInt < 0) {
                    J[j][i] = -oppoInt;
                } else if (J[i][j] < 0 & oppoInt > 0) {
                    J[j][i] = oppoInt;
                }
            }

        } 
    }
} 

void createTraits(double (&traits)[numSpec][2], mt19937& eng, double ppProb) {
    std::lognormal_distribution<double> distribution(-3.0, 3.0);
    for (int i = 0; i < numSpec; i++) {
        traits[i][0] = distribution(eng);
        if(ppProb < uniform(eng)) {traits[i][1] = 1;} else {traits[i][1] = 0;};
    }
}

void createDisp(double (&disp)[numSpec]) {
    for (int i = 0; i < numSpec; i++) {
        disp[i] = dispDist;
    }
}

void removeSpecies(vector <int> (&cellPopSpec)[numCells][4], int cell, int chosenInd) {

    for (int j = 0; j < 3; j++) {cellPopSpec[cell][j].erase(cellPopSpec[cell][j].begin() + chosenInd);}

}

// Function to add a species which may not already exist in the cell - immigration, mutation, initialisation
void addSpecies(vector <int> (&cellPopSpec)[numCells][4], int cell, int chosenInd, double (&traits)[numSpec][2]) {

    cellPopSpec[cell][0].push_back(chosenInd);
    cellPopSpec[cell][1].push_back(traits[0][chosenInd]);
    cellPopSpec[cell][2].push_back(traits[1][chosenInd]);
    cellPopSpec[cell][3].push_back(1); // Placeholder for dispersal ability

}

// Function for when a species already exists and is being replicated (maybe with intraspecific variation)
void baby(vector <int> (&cellPopSpec)[numCells][4], int cell, int chosenInd, double (&traits)[numSpec][2]) {

    int spec = cellPopSpec[cell][0][chosenInd];

    cellPopSpec[cell][0].push_back(spec);
    cellPopSpec[cell][1].push_back(traits[0][spec]);
    cellPopSpec[cell][2].push_back(traits[1][spec]);
    cellPopSpec[cell][3].push_back(1); // Placeholder for dispersal ability

}

void calculateTotalPopSpec(vector <int> (&cellPopSpec)[numCells][4], vector <int> (&totalPopSpec)[2]) {

    totalPopSpec[0].clear(); totalPopSpec[1].clear();

    for (int i = 0; i < numCells; i++) {
        for (int j = 0; j < cellPopSpec[i][0].size(); j++) {
            bool exists = false;
            for (int k = 0; k < totalPopSpec[0].size(); k++) {
                if(cellPopSpec[i][0][j] == totalPopSpec[0][k]) {
                    totalPopSpec[1][k] += 1;
                    exists = true;
                }
            }
            if(exists == false) {
                totalPopSpec[0].push_back(cellPopSpec[i][0][j]);
                totalPopSpec[1].push_back(1);
            }
        }
    }
}

int calcCellPop(vector <int> (&cellPopSpec)[numCells][4], int cell) {

    int pop = cellPopSpec[cell][0].size();

    return pop;
}

void calculateCellPop(vector <int> (&cellPopSpec)[numCells][4], vector <int> (&cellPop)[2]) {

    cellPop[0].clear(); cellPop[1].clear();

    for (int i = 0; i < numCells; i++) {
        int pop = cellPopSpec[i][0].size();
        if(pop > 0) {
            cellPop[0].push_back(i);
            cellPop[1].push_back(pop);
        }
    }
}

void storeCellPopSpec(ofstream &stream, vector <int> (&cellPopSpec)[numCells][4], int gen, double (&traits)[numSpec][2]) {

    for (int i = 0; i < numCells; i++) {
        for (int j = 0; j < cellPopSpec[i][0].size(); j++) {
            stream << gen+1 << " " << i+1 << " " << cellPopSpec[i][0][j] + 1 << " " << cellPopSpec[i][1][j] << " " << 
            traits[cellPopSpec[i][0][j]][0] << "\n";
        }
    }
}

// void storeCellPopSpec(ofstream &stream, vector <int> (&cellPopSpec)[numCells][4], int gen, double (&traits)[numSpec][2]) {

//     for (int i = 0; i < numCells; i++) {
//         for (int j = 0; j < cellPopSpec[i][0].size(); j++) {
//             stream << gen+1 << " " << i+1 << " " << cellPopSpec[i][0][j] + 1 << " " << cellPopSpec[i][1][j] << " " << 
//             traits[cellPopSpec[i][0][j]][0] << " " << calculateInteractions(cellPopSpec, traits, i, numSpec, cellPopSpec[i][0][j]) << " " 
//             << (getCellMass(i)/Rfr) << " " << std::pow(traits[cellPopSpec[i][0][j]][0], 0.75) << " " << 
//             exp(calculateInteractions(cellPopSpec, traits, i, numSpec, cellPopSpec[i][0][j]) - 
//             (getCellMass(i)/Rfr) - std::pow(traits[cellPopSpec[i][0][j]][0], 0.75))/(1 + exp(calculateInteractions(cellPopSpec, traits, i, numSpec, cellPopSpec[i][0][j]) - 
//             (getCellMass(i)/Rfr) - std::pow(traits[cellPopSpec[i][0][j]][0], 0.75))) << "\n";
//         }
//     }
// }

int getPop(int ind, int cell, vector <int> (&cellPopSpec)[numCells][4]) {

    int pop;

    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        if(ind == cellPopSpec[cell][0][i]) {
            pop = cellPopSpec[cell][1][i];
            break;
        }
    }
    
    return pop;

}

///////////////////////
// Mutation functions
///////////////////////

void Bin_recursive(int n, int j, int (&vec)[L]) { // It must be called setting j=L.
    j--;
    if (n / 2 != 0) {
        Bin_recursive(n/2, j, vec);
    }
    vec[j] = n % 2;
}


void ConvertToBinary(int n, int j, int (&vec)[L]) {
    int i;
    if(n >= std::pow(2,L)){
        std::cout << "You have set the number of species greater than 2^" << L << endl;
        exit (EXIT_FAILURE);
    }
    for (i=0; i<L; i++) {
        vec[i]=0;
    }
    Bin_recursive(n, j, vec);
}


int ConvertToDec(int (&arr)[L]) {
    int i, dec=0;
    for (i=0; i<L;i++) {
        if(arr[L-1-i] == 1){
            dec += std::pow(2,i);
        }
    }
    return dec;
}

int mutation(vector <int> (&cellPopSpec)[numCells][4], double prob, int chosenInd, mt19937& eng) {
    
    int b1[L]; // To store the binary sequence
    int mutSpec; // Integer identifier of mutated species that will be created

    ConvertToBinary(chosenInd, L, b1);

    for (int i = 0; i < L; i++) {
        if (uniform(eng) <= prob) {
            if (b1[i] == 1) {
                b1[i] = 0;
            } else {
                b1[i] = 1;
            }
        }
    }

    mutSpec = ConvertToDec(b1); // Convert our mutated species to the binary representation

    return mutSpec;

}

///////////////////////
// Metabolic Theory Functions
///////////////////////

double searchRate(int Si, int Sj, double E, double (&traits)[numSpec][2]) {
    
    double V0 = 0.33; double d0 = 1.62; 
    double a; double Mi = traits[Si][0];
    double Mj = traits[Sj][0];

    a = 2*V0*d0*(std::pow(Mi, 0.63))*(std::pow(Mj, 0.21))*arrhenius(E);

    return a;
    
}

double attackProb(int Si, int Sj, double (&traits)[numSpec][2]) {

    double Rp = 0.1;
    double Mi = traits[Si][0];
    double Mj = traits[Sj][0];

    double A = (1/(1 + 0.25*(std::pow(exp(1), -std::pow(Mi, 0.33)))))*std::pow(1/(1 + std::pow(log10(Rp*(Mi/Mj)), 2)),0.2);

    return A;

}

double handlingTime(int Si, int Sj, double E, double (&traits)[numSpec][2]) {

    double h0 = 1;
    double Rp = 0.1;
    double Mi = traits[Si][0];
    double Mj = traits[Sj][0];

    double h = h0*std::pow(Mi, -0.75)*(1-exp(-(std::pow((Mj/Mi) - Rp, 2))/2))*arrhenius(E);

    return h;

}

double consumptionRate(int Si, int Sj, double E, double (&traits)[numSpec][2], int Nj) {

    double c = (searchRate(Si, Sj, E, traits)*attackProb(Si, Sj, traits)*Nj*traits[Sj][0])/
    (1 + (searchRate(Si, Sj, E, traits)*attackProb(Si, Sj, traits)*handlingTime(Si, Sj, E, traits)*Nj));

    return c;

}

double arrhenius(double E) {

    return std::pow(exp(1), -E/(k*(T+T0)));

}

double getCellMass(int cell, vector <int> cellPopSpec[numCells][4]) {

    double cellMass = 0;

    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        cellMass += cellPopSpec[cell][1][i]; 
    }

    return cellMass;

}

void storeConsumptionRate(ofstream &stream, vector <int> (&cellPopSpec)[numCells][4], int gen, double (&traits)[numSpec][2]) {

    vector<int> optimalSet;

    for (int i = 0; i < numCells; i++) {
        for (int j = 0; j < cellPopSpec[i][0].size(); j++) {         
            // Find the set of species which are in the optimal set for focal species Si (ind)
            optimalForaging(cellPopSpec[i][0][j], i, traits, cellPopSpec, optimalSet);
            // Iterate over set and save consumption rate
            for (int k = 0; k < optimalSet.size(); k++) {
                int Sj = optimalSet[k];
                stream << gen+1 << " " << i+1 << " " << cellPopSpec[i][0][j] + 1 << " " << traits[cellPopSpec[i][0][j]][0] << " " << cellPopSpec[i][1][j]
                << " " << (Sj + 1) << " " << traits[Sj][0] << " " << getPop(Sj, i, cellPopSpec) << " " << 
                searchRate(cellPopSpec[i][0][j], Sj, 0, traits) << " " << attackProb(cellPopSpec[i][0][j], Sj, traits) << 
                " " << handlingTime(cellPopSpec[i][0][j], Sj, 0, traits) << " " << 
                consumptionRate(cellPopSpec[i][0][j], Sj, 0, traits, getPop(Sj, i, cellPopSpec)) << "\n";;
            }
        }
    }

}

double profitability(int Si, int Sj, double (&traits)[numSpec][2]) {
    
    double Mj = traits[Sj][0];
    double pij;

    pij = (attackProb(Si, Sj, traits)*Mj)/handlingTime(Si, Sj, 0, traits);

    return pij;

}

void sortedCols(vector <double> (&profits)[2], vector <double> (&sortedProfits)[2]) {

    const int size = profits[0].size();

    // Creating a vector of indices to track the sorting order
    vector<int> indices(size);
    for (int i = 0; i < size; i++) {
        indices[i] = i;
    }

    // Sorting the indices based on the values in profits[1]
    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
        return profits[1][a] > profits[1][b];
    });

    // Rearranging the columns based on the sorted indices
    for (int i = 0; i < size; i++) {
        sortedProfits[0][i] = profits[0][indices[i]];
        sortedProfits[1][i] = profits[1][indices[i]];
    }

}

void findOptimalSet(int Si, int cell, vector <double> (&sortedProfits)[2], double (&traits)[numSpec][2], vector <int> (&optimalSet)) {

    int Sj, pop;
    double OF = 0, OFnew = 0;
    double upper = 0; double lower = 0;
    double search, attack;
    optimalSet.clear(); // Clear current optimal set so we don't add on

    // We run through each set of species comparing against the previous set
    // e.g set 1 is just Sj = 1, then set 2 is Sj = {1,2}, set 3 is Sj = {1,2,3} etc
    // To aid with computation time we simply add the additions from the new set to the raw values from the old set
    // e.g for set 2 we add upper and lower for Si = 2, to the upper and lower results for Si = 1.
    for (int i = 0; i < sortedProfits[0].size(); i++) {
        OF = OFnew;

        Sj = sortedProfits[0][i];
        pop = getPop(Sj, cell, cellPopSpec);
        search = searchRate(Si, Sj, 0, traits);
        attack = attackProb(Si, Sj, traits);

        upper += search*attack*pop*traits[Sj][0];
        lower += search*attack*handlingTime(Si, Sj, 0, traits)*pop;

        // Calculate the consumption rate with this set of resource species
        OFnew = upper/(1 + lower);

        // If the OF (optimal foraging) calculated prior to the addition of this species
        // is greater than the OF with this species, stop running, we have found the optimal set
        // Save this as the number of rows down the sortedProfit vector we've gone - 1, as the most recent
        // species added doesn't help
        // And trim sortedProfits to just this profitable set of species

        if(OFnew < OF || i == (sortedProfits[0].size() - 1)) {
            int optimal = i;
            for (int j = 0; j <= optimal; j++) {
                optimalSet.push_back(sortedProfits[0][j]); // Store species ID
                }
            break;
        }
    }    
}

void optimalForaging(int Si, int cell, double (&traits)[numSpec][2], vector <int> (&cellPopSpec)[numCells][4], vector <int> (&optimalSet)) {

    vector <double> profits[2];
    vector<double> sortedProfits[2];
    int Sj;

    // Calculate profitability of consumer Si feeding on
    // each resource Sj
    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        Sj = cellPopSpec[cell][0][i]; 
        // Don't calculate profitability for a species with itself
        if(Si != Sj) {
            profits[0].push_back(Sj);
            profits[1].push_back(profitability(Si, Sj, traits));
        }
    }

    // Now we know how big sortedProfits needs to be
    // so we can preassign it
    sortedProfits[0].resize(profits[0].size());
    sortedProfits[1].resize(profits[1].size());

    // Sort in descending order of profitability
    // Store output in sortedProfits
    sortedCols(profits, sortedProfits);

    // Next we need to optimise the optimal foraging equation
    findOptimalSet(Si, cell, sortedProfits, traits, optimalSet); 

}

void updateSeachRate(vector<vector<double>> searchRates[numCells], int cell, int chosenInd) {

    // searchRates

}