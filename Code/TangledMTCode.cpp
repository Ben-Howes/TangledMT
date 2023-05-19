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

const int numSpec = 10000;                    // Number of species in the model (each species will have interacitons and mass associated with it).
const int numGenPre = 1000;                 // Number of generations to run the model pre fragmentation. Each generation is broken into time steps.
const int numGenPost = 10000;                //2000 Number of generations to run the model after fragmentation. Each generation is broken into time steps.
const int initPop = 50;                    // Number of individuals to put into each cell at the start of the model.

const float probDeath = 0.15;               // Probability of individual dying if chosen.
double probImm = 0.00001;                      // Probability of an individual immigrating into a cell (individual of a random species suddenly occurring in a cell).
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
vector <int> cellPop[2];                    // Stores the total population in each cell at a given generation.
vector <int> totalPopSpec[2];                      // Stores the total population of each species at a given generation.
vector <int> cellPopSpec[3];                       // Stores the total population of each species in each cell at a given generation.
int totalRich = 0;                          // Stores the total species richness of the model at a given generation.
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
void initialisePop(vector <int> (&cellPopSpec)[3], int numSpec, int numCells, mt19937& eng);
int chooseInRange(int a, int b, mt19937& eng);
void immigration(vector <int> (&cellPopSpec)[3], double prob, int cell, int numSpec, int &immNum, mt19937& eng);
void kill(vector <int> (&cellPopSpec)[3], double prob, int cell, int numSpec, mt19937& eng);
void reproduction(vector <int> (&cellPopSpec)[3], int (&cellList)[numCells][2], double (&J)[numSpec][numSpec], int cell, int numSpec, int Rfr, mt19937& eng);
double calculateInteractions(double (&J)[numSpec][numSpec], vector <int> (&cellPopSpec)[3], int cell, int numSpec, int ind);
void shuffle(int arr[], int arrElements, mt19937& eng);
int randomInd(vector <int> (&cellPopSpec)[3], int pop, int numSpec, int cell, mt19937& eng);
void cellCoords(int (&landscapeCoords)[numCells][2], int cols, int rows, int numCells);
void getDistances(double (&distArray)[numCells][numCells], int landscapeCoords[][2], int cellRows, int cellCols);
void dispersal(vector <int> (&cellPopSpec)[3], double (&distArray)[numCells][numCells], double (&disp)[numSpec], double prob, int cell, int numSpec, int &dispNum, mt19937& eng); 
vector<int> findValidCells(double (&distArray)[numCells][numCells], double distance, int cell);
double dispersalProb(vector <int> (&cellPopSpec)[3], int Rfr, int cell, double probDeath, double probDisp);
void fillCellList(int landscapeArray[][cellRows], int cellList[][2], int cellCols, int cellRows);
void getForestCellList(int cellList[][2], vector <int> &forestCellList);
void removeSpecies(vector <int> (&cellPopSpec)[3], int cell, int chosenInd);
void addSpecies(vector <int> (&cellPopSpec)[3], int cell, int chosenInd);

void store2ColFiles(ofstream &stream, int firstCol, int secondCol);
void storeVec(ofstream &stream, vector <int> vec[], int gen, int cols);
void storeVecEnd(vector <int> vec[], int cols, string fileName, string outpath);
void storeNum(int num, string fileName, string outpath);
void storeParam(string fileName, string outpath);

void calculateTotalPopSpec(vector <int> (&cellPopSpec)[3], vector <int> (&totalPopSpec)[2]);
void calculateCellPop(vector <int> (&cellPopSpec)[3], vector <int> (&cellPop)[2]);
int getPop(vector <int> (&cellPopSpec)[3], int cell);

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
        // read2DArray<int, numSpec>(cellPopSpec, numSpec, numCells, "/final_cellPopSpec.txt", finpath);
        
       
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

        // Initialise population and store the initial population to file
        initialisePop(cellPopSpec, numSpec, numCells, eng);
        // store2DArray<int, numSpec>(cellPopSpec, numSpec, numCells, "/initial_Pop.txt", outpath);

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
                kill(cellPopSpec, probDeath, cell, numSpec, eng);
                reproduction(cellPopSpec, cellList, J, cell, numSpec, Rfr, eng);
                probDispDen = dispersalProb(cellPopSpec, Rfr, cell, probDeath, probDisp);
                dispersal(cellPopSpec, distArray, disp, probDispDen, cell, numSpec, dispNum, eng);
                if(fragmentation == 0) {
                    immigration(cellPopSpec, probImm, cell, numSpec, immNum, eng);
                } else {
                    immigration(cellPopSpec, probImmFrag, j, numSpec, immNum, eng);
                }
            }
        }

        // Calculate richness and abundance metrics
        // Dependent on how often you want to save them
        if((i+1)%10 == 0) {
        totalPop = 0;
        for (int i = 0; i < cellPopSpec[0].size(); i++) {totalPop += cellPopSpec[2][i];}
        calculateTotalPopSpec(cellPopSpec, totalPopSpec);
        totalRich = totalPopSpec[0].size();
        calculateCellPop(cellPopSpec, cellPop);

        // After each generation store all outputs
        store2ColFiles(s_totalPop, i, totalPop);
        store2ColFiles(s_totalRich, i, totalRich);
        storeVec(s_cellPopSpec, cellPopSpec, i, 3);
        storeVec(s_totalPopSpec, totalPopSpec, i, 2);
        storeVec(s_cellPop, cellPop, i, 2);

        // Live output to console
        cout << "Gen: " << i + 1 << "/" << numGen << " | Total Pop: " << totalPop << " | Total Richness: " << totalRich << "\n";

        }
    }
    // Store final outputs of burn ins
    if(fragmentation == 0) {
        storeNum(immNum, "/total_Immigrations.txt", finpath);
        storeNum(dispNum, "/total_Dispersals.txt", finpath);
        storeVecEnd(cellPopSpec, 3, "/final_cellPopSpec.txt", finpath);
    }
    // Store final outputs for fragment runs
    if(fragmentation == 1) {
        // save cellpopspec
        storeNum(immNum, "/total_Immigrations.txt", finfragpath);
        storeNum(dispNum, "/total_Dispersals.txt", finfragpath);
        storeVecEnd(cellPopSpec, 3, "/final_cellPopSpec.txt", finfragpath);
    }

    // Close all streams to files
    s_totalPop.close(); s_totalRich.close(); s_cellPop.close(); s_totalPopSpec.close(); s_cellPopSpec.close();

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

void storeCellPopSpec(ofstream &stream, int gen, vector <int> cellPopSpec[3]) {
    for (int i = 0; i < cellPopSpec[0].size(); i++) {
        stream << gen+1 << " " << cellPopSpec[0][i] + 1 << " " << cellPopSpec[1][i] + 1 << " " << cellPopSpec[2][i] << "\n";
    }
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

double dispersalProb(vector <int> (&cellPopSpec)[3],int Rfr, int cell, double probDeath, double probDisp) {

    double Nequ, pDispEff, cutOff, grad;
    int pop = getPop(cellPopSpec, cell);
    pDispEff = probDisp;

    cutOff = Rfr*(3+log10((1-probDeath)/probDeath));        // Cut off for density dependent linear, max is 47.5
    // Nequ = Rfr*(3+log((1-probDeath)/probDeath));			//calculate density dependent Pmig if N/Nequ>.75

    // Make density dependence linear with a cut-off
    
    grad = 1/cutOff;
    if(pop > cutOff) {pDispEff = probDisp;} else {pDispEff = (grad*pop)*probDisp;}
    
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

void dispersal(vector <int> (&cellPopSpec)[3], double (&distArray)[numCells][numCells], 
double (&disp)[numSpec], double prob, int cell, int numSpec, int &dispNum, mt19937& eng) {

    int pop = getPop(cellPopSpec, cell);

    if(pop > 0) { // Check cell isn't empty

        int chosenInd, dispCell;
        double dispAbility; // Store dispersal ability of chosen species
        vector<int> validCells; // store cells that individual could disperse to

        if (uniform(eng) <= prob) {
            chosenInd = randomInd(cellPopSpec, pop, numSpec, cell, eng);
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

                // Remove individual from cell it was in
                removeSpecies(cellPopSpec, cell, chosenInd);
                // Add individual to dispCell
                addSpecies(cellPopSpec, dispCell, chosenInd);

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

double calculateInteractions(double (&J)[numSpec][numSpec], vector <int> (&cellPopSpec)[3], int cell, int numSpec, int ind) {

    double H;
    H = 0; 

    for (int i = 0; i < cellPopSpec[0].size(); i++) {
        // Check that they are in the same cell
        if(cellPopSpec[0][i] == cell) {
            // Find J matrix interaction between ind and the other species in the cell
            // and multiply by the abundance of the other species
            H += J[ind][cellPopSpec[1][i]]*cellPopSpec[2][i];
        }
    }

    return H;

}

void reproduction(vector <int> (&cellPopSpec)[3], 
int (&cellList)[numCells][2], double (&J)[numSpec][numSpec], int cell, int numSpec, int Rfr, mt19937& eng) {

    int pop = getPop(cellPopSpec, cell);

    if(pop > 0) {

        int chosenInd;
        double H, pOff;

        chosenInd = randomInd(cellPopSpec, pop, numSpec, cell, eng);
        H = calculateInteractions(J, cellPopSpec, cell, numSpec, chosenInd);
        H = (H*weightInt/pop) - (pop/Rfr); // 10 chosen as arbitrary carrying capactiy
        pOff = exp(H) / (1 + exp(H));

        if (uniform(eng) <= pOff) {
            addSpecies(cellPopSpec, cell, chosenInd);
        }
    }
}

void kill(vector <int> (&cellPopSpec)[3], double prob, int cell, int numSpec, mt19937& eng) {
    
    int pop = getPop(cellPopSpec, cell);

    if(pop > 0){ 

    int chosenInd;

        if (uniform(eng) <= prob) {
            chosenInd = randomInd(cellPopSpec, pop, numSpec, cell, eng);
            removeSpecies(cellPopSpec, cell, chosenInd);
        }
    }
}

// Choose a random individual (proportional to how many of that species there are)
int randomInd(vector <int> (&cellPopSpec)[3], int pop, int numSpec, int cell, mt19937& eng) {
    
    double sum, threshold;
    sum = 0;
    // Create threshold value by multiplying the total number of individuals in the cell by a number between 0-1
    threshold = uniform(eng)*pop;

    // Get an individual of a random species (proportional to how many of that species there are) by summing the number of individuals of each species in the cell
    // and then choosing the species whose population get the sum to the threshold value
    for (int i = 0; i < cellPopSpec[0].size(); i++) {
        if(cellPopSpec[0][i] == cell) {
            sum += cellPopSpec[2][i];
            if (sum > threshold) {return cellPopSpec[1][i];}
        }
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

void immigration(vector <int> (&cellPopSpec)[3], double prob, int cell, int numSpec, int &immNum, mt19937& eng) {
    
    int chosenSpec;
    bool exists = false;

    if(uniform(eng) <= prob) {  
        chosenSpec = chooseInRange(0, numSpec-1, eng);
        addSpecies(cellPopSpec, cell, chosenSpec);
        immNum++;
    }
}

int chooseInRange(int a, int b, mt19937& eng) {
    uniform_int_distribution<> choose(a, b); // define the range [a,b], extremes included.
    return choose(eng);
}

void initialisePop(vector <int> (&cellPopSpec)[3], int numSpec, int numCells, mt19937& eng) {

    int chosenSpec;
    bool exists = false;

    for (int i = 0; i < numCells; i++){
        for (int j = 0; j < initPop; j++) {
            chosenSpec = chooseInRange(0, numSpec-1, eng); 
            addSpecies(cellPopSpec, i, chosenSpec);
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

void removeSpecies(vector <int> (&cellPopSpec)[3], int cell, int chosenInd) {
    for (int i = 0; i < cellPopSpec[0].size(); i++) {
        if(cellPopSpec[0][i] == cell & cellPopSpec[1][i] == chosenInd) {
            if(cellPopSpec[2][i] > 1) {
                cellPopSpec[2][i] -= 1;
                break;
            } else if (cellPopSpec[2][i] == 1){
                for (int j = 0; j < 3; j++) {cellPopSpec[j].erase(cellPopSpec[j].begin() + i);}
                break;
            }
        }
    }
}

 void addSpecies(vector <int> (&cellPopSpec)[3], int cell, int chosenInd) {

    bool exists = false;

    for (int i = 0; i < cellPopSpec[0].size(); i++) {
        if(cellPopSpec[0][i] == cell & cellPopSpec[1][i] == chosenInd) {
            cellPopSpec[2][i] +=1;
            exists = true;
            break;
        }
    }
    if(exists == false) {  
        cellPopSpec[0].push_back(cell);
        cellPopSpec[1].push_back(chosenInd);
        cellPopSpec[2].push_back(1);
    }
}

void calculateTotalPopSpec(vector <int> (&cellPopSpec)[3], vector <int> (&totalPopSpec)[2]) {

    totalPopSpec[0].clear(); totalPopSpec[1].clear();

    for (int i = 0; i < cellPopSpec[0].size(); i++) {
        bool exists = false;
        for (int j = 0; j < totalPopSpec[0].size(); j++) {
            if(cellPopSpec[1][i] == totalPopSpec[0][j]) {
                totalPopSpec[1][j] += cellPopSpec[2][i];
                exists = true;
                break;
            } 
        }
        if(exists == false) {
            totalPopSpec[0].push_back(cellPopSpec[1][i]);
            totalPopSpec[1].push_back(cellPopSpec[2][i]);
        }
    }
}

int getPop(vector <int> (&cellPopSpec)[3], int cell) {

    int pop = 0;

    for (int i = 0; i < cellPopSpec[0].size(); i++) {
        if(cellPopSpec[0][i] == cell) {
            pop += cellPopSpec[2][i];
        }   
    }
    return pop;
}

void calculateCellPop(vector <int> (&cellPopSpec)[3], vector <int> (&cellPop)[2]) {

    cellPop[0].clear(); cellPop[1].clear();

    for (int i = 0; i < cellPopSpec[0].size(); i++) {
        bool exists = false;
        for (int j = 0; j < cellPop[0].size(); j++) {
            if(cellPopSpec[0][i] == cellPop[0][j]) {
                cellPop[1][j] += cellPopSpec[2][i];
                exists = true;
                break;
            } 
        }
        if(exists == false) {
            cellPop[0].push_back(cellPopSpec[0][i]);
            cellPop[1].push_back(cellPopSpec[2][i]);
        }
    }

}