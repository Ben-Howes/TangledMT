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
const int numGenPre = 100;                 // Number of generations to run the model pre fragmentation. Each generation is broken into time steps.
const int numGenPost = 100;                //2000 Number of generations to run the model after fragmentation. Each generation is broken into time steps.
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

void initialisePop(vector <int>(&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], 
    int &totalRich, int numSpec, int numCells, mt19937& eng);
int chooseInRange(int a, int b, mt19937& eng);
void kill(vector <int> (&totalPopSpec)[2], int (&cellPopSpec)[numCells][numSpec], int &totalPop, int (&cellPop)[numCells],  int (&cellRich)[numCells], int &totalRich,
double prob, int cell, int numSpec, mt19937& eng);
int randomInd(int (&cellPopSpec)[numCells][numSpec], int (&cellPop)[numCells], int numSpec, int cell, mt19937& eng);
double uniform(mt19937& eng);

int main(int argc, char *argv[]) {

  mt19937 eng(seed); // Seed random numbers

for (int i = 0; i < 2; i++) {totalPopSpec[i].resize(0);}


// Initialise population and store the initial population to file
    initialisePop(totalPopSpec, cellPopSpec, totalPop, cellPop, cellRich, totalRich, numSpec, numCells, eng);
    for (int cell = 0; cell < 121; cell++){
        for (int i = 0; i < 720; i++) {
            kill(totalPopSpec, cellPopSpec, totalPop, cellPop,  cellRich, totalRich, probDeath, cell, numSpec, eng);
        }
    }
    


    ofstream out;
    out.open("totalPopSpec.txt");

    for (int i = 0; i < totalPopSpec[0].size(); i++) {
        for (int j = 0; j < 2; j++) {
            out << totalPopSpec[j][i] << " ";
        }
        out << endl;
    }

        out.close();
    
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
            }

            if(cellPopSpec[i][chosenSpec] == 0) {cellRich[i]++;}
            cellPopSpec[i][chosenSpec]++;
            totalPop++;
            cellPop[i]++;
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

double uniform(mt19937& eng) {
    uniform_real_distribution<> rand(0, 1); // define the range [a,b], extremes included.
    return rand(eng);
}