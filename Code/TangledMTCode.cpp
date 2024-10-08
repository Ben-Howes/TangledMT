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
#include <omp.h>

using namespace std;

////////////////////////////
// User-Defined Variables //
////////////////////////////

// These variables can be easily changed to alter the output of the model

string dir = "../Results/TNM_Output/";   // Directory for output of model
int habitatLoss = 0;                      // If 0 the landscape stays the same, if 1 then the landscape looses habitat
int climateChange = 0;                      // if 0 the temperature stays the same, if 1 then the temperature changes

int defSeed = 1;                            // This is the default seed that will be used if one is not provided in the command line argument (recommend using command line

const int cellRows = 1;                    //Sets the number of cells in the rows of the landscape (Note: must match landscape file used in code, if using file).
const int cellCols = 1;                    // Sets the number of cells in the columns of the landscape (Note: must match landscape file used in code).

const int L = 10;                           // Length of binary identifiers to use in the model (genome sequences)
const int numSpec = 1024;                    // Number of species in the model, the number of species must equal 2^L .
const int t = 5000000;                         // Number of time steps in the model
int lossT = t/2;                            // Time step at which habitat loss occurs if habitatLoss == 1
int climateT = t/2;                         // Time step at which climate change occurs if climateChange == 1
const int initPop = numSpec;                // Number of individuals to put into each cell at the start of the model.

const float probDeath = 0.2;               // Probability of individual dying if chosen.
double probImm = 0.00001;                      // Probability of an individual immigrating into a cell (individual of a random species suddenly occurring in a cell).
float probDisp = 0.001;                       // Probability of an individual dispersing from one cell to another cell. 
double probMut = 0;                      // Probability of a number in the genome sequence switching from 0 -> 1, or 1 -> 0
double propLoss = 0.5;                  // Proportion of cells that are lost if habitatLoss == 1 (0 - 1)


// Metabolic theory variables
double ppProp = 0.2;                       // Sets proportion of species that are primary producers
double T = 20;                             // Set temperature in kelvin (273.15 kelvin = 0 celsius)
double changeT = 5;                        // Amount to change temperature by if climateChange == 1
double E = 0.65;                           // Activation energy
double k = 8.6173*pow(10, -5);                  // Boltzmann constant
double r0 = 100;                              // Multiplier for gain in mass of primary producers
double K0 = 1000;                             // Weighting of carrying capacity of each primary producing species Ki. Increased K0 increases primary producer abundance linearly.
double I0 = 0.5;                         // Constant affecting the influence of interference (intraspecific competition), higher I0 = higher intraspecific competition
double G0;                               // Store normalising constant for generation time which will be equal to 1/(mass of smallest species in pool)^0.25
double D0;                               // Store normalising constant for dispersal distance, which will be based on the mass of the smallest spcies in the pool
double alpha = 0.4;                        // Sets slope of pOff, and therefore timescale, you want alpha to be large enough that species never hit their maximum pOff
double H0 = log((1/(1*pow(10, -10)))-1)/alpha;      // Offset value so pR = a negligible value (1x10^-10) when H = 0

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
int totalPop = 0;                           // Stores the total population across all cells in the model at a given generation.
vector <int> cellPop[2];                    // Stores the total population in each cell at a given generation.
vector <int> totalPopSpec[2];                      // Stores the total population of each species at a given generation.
vector <double> cellPopSpec[numCells][3];                 // Stores 0 = species, 1 = abundance of species in cell, 2 = biomass of species in cell
vector <double> cellPopInd[numCells][3];             // Stores which individuals are in which cells (0 = species, 1 = mass, 2 = primary producer
int totalRich = 0;                          // Stores the total species richness of the model at a given generation.
int cellList[numCells][2];                     // Lists numbers of cell (e,g with 6 cells it would read, 0,1,2,3,4,5), and whether they are non-forest or forest (0 = non-forest, 1 = forest).
vector <int> forestCellList;                    // Lists cells that are forested.
int cellOrder[numCells];                        // Array to store the randomly selected order of cells for dynamics.
double distMatrix[numCells][numCells];          // Stores distances between all cells in the landscape.
int seed;
int immNum = 0;                                 // Counts number of immigrations that occur
int dispNum = 0;                                // Counts number of dispersals that occur
double maxDisp = 0;                             // Stores maximum dispersal distance        

double pMax = 0;                                   // Count the number of times the pOff of a species is within 1% of its maximum

static const double two_pi  = 2.0*3.141592653;

// Metabolic theory variables
static double T0 = 273.15;                        // Kevlin when Celsius = 0, to use as offset with T

//////////////////////////
// Initialise functions //
//////////////////////////

void createTraits(double (&traits)[numSpec][2], mt19937& eng, double ppProp);
double uniform(mt19937& eng);
double gaussian(mt19937& eng);
void initialisePop(vector <double> (&cellPopInd)[numCells][3], double (&traits)[numSpec][2], int numSpec, int numCells, mt19937& eng);
int chooseInRange(int a, int b, mt19937& eng);
void immigration(vector <double> (&cellPopInd)[numCells][3], double prob, int cell, int numSpec, int &immNum, mt19937& eng);
void kill(vector <double> (&cellPopInd)[numCells][3], double prob, int cell, int numSpec, mt19937& eng);
void reproduction(vector <double> (&cellPopInd)[numCells][3], 
int (&cellList)[numCells][2], double (&traits)[numSpec][2], int cell, int numSpec, int gen, mt19937& eng);
double calculateInteractions(vector <double> (&cellPopInd)[numCells][3], double (&traits)[numSpec][2], int cell, int numSpec, int ind,
    vector <double> (&cellPopSpec)[numCells][3], int gen);
void shuffle(int arr[], int arrElements, mt19937& eng);
int randomIndex(vector <double> (&cellPopInd)[numCells][3], int pop, int numSpec, int cell, mt19937& eng);
void cellCoords(int (&landscapeCoords)[numCells][2], int cols, int rows, int numCells);
void getDistances(double (&distArray)[numCells][numCells], int landscapeCoords[][2], int cellRows, int cellCols);
void dispersal(vector <double> (&cellPopInd)[numCells][3], double (&distArray)[numCells][numCells], double prob, int cell, int numSpec, int &dispNum, mt19937& eng); 
vector <int> findValidCells(double (&distArray)[numCells][numCells], int (&cellList)[numCells][2], double distance, int cell);
void fillCellList(int landscapeArray[][cellRows], int cellList[][2], int cellCols, int cellRows);
void getForestCellList(int cellList[][2], vector <int> &forestCellList);
void removeInd(vector <double> (&cellPopInd)[numCells][3], int cell, int chosenIndex, vector <double> (&cellPopSpec)[numCells][3]);
void addInd(vector <double> (&cellPopInd)[numCells][3], int cell, int chosenSpec, double (&traits)[numSpec][2], vector <double> (&cellPopSpec)[numCells][3]);
void loseHabitat(int (&landscapeArray)[cellCols][cellRows], int (&cellList)[numCells][2], int cellRows, int cellCols, int numCells, mt19937& eng, double propLoss);

void store2ColFiles(ofstream &stream, int firstCol, int secondCol);
void storeVec(ofstream &stream, vector <int> vec[], int gen, int cols);
void storeVecEnd(vector <int> vec[], int cols, string fileName, string outpath);
void storeNum(int num, string fileName, string outpath);
void storeParam(string fileName, string outpath);
void storecellPopInd(ofstream &stream, vector <double> (&cellPopInd)[numCells][3], int gen);
void storeCellPopSpec(ofstream &stream, vector <double> (&cellPopSpec)[numCells][3], int gen, double (&traits)[numSpec][2]);

void calculateTotalPopSpec(vector <double> (&cellPopInd)[numCells][3], vector <int> (&totalPopSpec)[2]);
void calculateCellPop(vector <double> (&cellPopInd)[numCells][3], vector <int> (&cellPop)[2]);
void calculateCellPopSpec(vector <double> (&cellPopInd)[numCells][3], vector <double> (&cellPopSpec)[numCells][3]);

// Mutations
int mutation(vector <double> (&cellPopInd)[numCells][3], double prob, int chosenSpec, mt19937& eng);
int ConvertToDec(int (&arr)[L]);
void ConvertToBinary(int n, int j, int (&vec)[L]);
void Bin_recursive(int n, int j, int (&vec)[L]);

// Metabolic Theory Functions
double searchRate(double Mi, double Mj, double E, double (&traits)[numSpec][2]);
double attackProb(double Mi, double Mj, double (&traits)[numSpec][2]);
double handlingTime(double Mi, double Mj, double E, double (&traits)[numSpec][2]);
double consumptionRate(double Mi, double Mj, double E, double (&traits)[numSpec][2], int Nj);
double getCellMass(int cell, vector <double> (&cellPopInd)[numCells][3]);
double getSpeciesCellMass(int cell, int chosenSpec, vector <double> (&cellPopInd)[numCells][3]);
double arrhenius(double E);
void storeConsumptionRate(ofstream &stream, vector <double> (&cellPopInd)[numCells][3], int gen, double (&traits)[numSpec][2]);
void storeReproduction(ofstream &stream, vector <double> (&cellPopInd)[numCells][3], int gen, double (&traits)[numSpec][2]);
double profitability(int Si, int Sj, double (&traits)[numSpec][2]);
double dispersalDist(double Mi, double E, mt19937& eng);

///////////////////////////////
//          Templates        //
///////////////////////////////

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
       
    /////////////
    // Check command Line
    /////////////
    
    if(argc == 2) {seed = atof(argv[1]);} else {seed = defSeed;}
    
    if(argc == 5) {
        seed = atof(argv[1]);
        r0 = atof(argv[2]);
        K0 = atof(argv[3]);
        I0 = atof(argv[4]);
    }

    if(argc == 7) {
        seed = atof(argv[1]);
        r0 = atof(argv[2]);
        K0 = atof(argv[3]);
        I0 = atof(argv[4]);
        habitatLoss = atof(argv[5]);
        propLoss = atof(argv[6]);
    }

    if(argc != 2 && argc != 5 && argc != 7) {
        std::cout << "Incorrect number of arguments entered" << endl;
        std::cout << "Either input a seed, or a seed and arguments for primary producer gain multiplier, carrying capacity, and competition" << endl;
        std::cout << "Or input a seed, arguments for primary producer gain multiplier, carrying capacity, competition, habitat loss (0 or 1), and proportion of cells lost" << endl;
        exit(0);
    }

    // Set working directories based on the seed taken from the command line argument
    string outpath = dir + "Seed_" + to_string(seed);       // Folder that outputs are saved to.
    string respath = outpath + "/Results";                  // Folder for output of results.
    string landpath = outpath + "/Landscape";                // Folder for output of landscape files.
    string finpath = outpath + "/Final_Output";              // Folder for final generation outputs for use in future models.

    // Make new directory at output location named after the seed

    DIR* checkDir = opendir(outpath.c_str());
    if (checkDir) {
        // closedir(checkDir);
        std::cout << "Directory for seed " << seed << " already exists, overwriting....." << endl;
        // exit(0);
    } else if (!checkDir){
        mkdir(outpath.c_str(),07777);
        std::cout << "Directory does not exist, creating... at " << outpath << endl;
    }

    // Add more directories. 
    // One to store the final outputs for use in future models, another for landscape files, and another for the full results.
    mkdir(respath.c_str(), 07777);
    mkdir(finpath.c_str(), 07777);
    mkdir(landpath.c_str(), 07777);

    // Fill cellList with number of each cell (0,1,2,3 etc)
    // second column in cellList currently blank as it could change depending on the landscape (1 = forest, 0 = non-forest)
    for (int i = 0; i < numCells; i++) {cellList[i][0] = i; cellOrder[i] = i;}

    mt19937 eng(seed); // Seed random numbers

    // Create and store species traits
    // createJMatrix(J, probInt, eng);
    // store2DArray<double, numSpec>(J, numSpec, numSpec, "/JMatrix.txt", outpath);
    createTraits(traits, eng, ppProp);
    store2DArray<double, 2>(traits, 2, numSpec, "/traits.txt", outpath);

    // Calculate G0 based on minimum mass of species in the pool
    // and on the maximum temperature the landscape could be at
    double minMi = traits[0][0];
    for (int i = 0; i < numSpec; i++){if(traits[i][0] < minMi){minMi = traits[i][0];}}

    if(climateChange == 1) {
        T = T + changeT;
        G0 = 1/(pow(minMi, 0.25)*arrhenius(-E));
        T = T - changeT;
    } else {
        G0 = 1/(pow(minMi, 0.25)*arrhenius(-E));
    }

    // Calculate D0 based on minimum mass of species in the pool
    // Calculate the value of D0 needed to make the smallest species have a 1% chance
    // of dispersing at least 1 cell distance
    // Do this based on maximum temperature of the landscape (if climateChange occurs)
    if(climateChange == 1) {
        T = T + changeT;
        D0 = 0.01/((pow(exp(1), (-1/1.5))*pow(minMi, 0.63)*arrhenius(E)));
        T = T - changeT;
    } else {
        D0 = 0.01/((pow(exp(1), (-1/1.5))*pow(minMi, 0.63)*arrhenius(E)));
    }

    // Store parameters used in model
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
    initialisePop(cellPopInd, traits, numSpec, numCells, eng);

    ofstream s_totalPop;
    s_totalPop.open(respath + "/totalPop.txt");
    ofstream s_totalRich;
    s_totalRich.open(respath + "/totalRich.txt");
    ofstream s_cellPop;
    s_cellPop.open(respath + "/cellPop.txt");
    ofstream s_totalPopSpec;
    s_totalPopSpec.open(respath + "/totalPopSpec.txt");
    ofstream s_cellPopInd;
    s_cellPopInd.open(respath + "/cellPopInd.txt");
    ofstream s_cellPopSpec;
    s_cellPopSpec.open(respath + "/cellPopSpec.txt");
    ofstream s_consumptionRate;
    s_consumptionRate.open(respath + "/consumptionRate.txt");
    ofstream s_reproduction;
    s_reproduction.open(respath + "/reproduction.txt");
    ofstream s_pMax;
    s_pMax.open(respath + "/pMax.txt");

    // Start model dynamics //
    ///////////////////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < t; i++) {
        
        // Check if this is the timestep when habitat loss occurs
        if (habitatLoss == 1 && i == lossT) {
            loseHabitat(landscapeArray, cellList, cellRows, cellCols, numCells, eng, propLoss);
            // Save the new cellList with cells showing as (0 = no habitat, 1 = habitat)
            store2DArray<int, 2>(cellList, 2, numCells, "/lossCellList.txt", landpath);
            // Save teh new landscape array so it's easy to visualise the loss to the landscape
            store2DArray<int, cellRows>(landscapeArray, cellCols, cellRows, "/lossLandscape.txt", landpath);
        }

        // Check if this is the timestep when climate change occurs
        if (climateChange == 1 && i == climateT) {
            // Update T to equal T + changeT
            T += changeT;
        }
        

        shuffle(cellOrder, numCells, eng);
        //Loop through each cell of the landscape
        // Going over each cell in parallel
        // Note that we only need to worry about dispersal as all other things are cell
        // specific anyway
        for (int j = 0; j < numCells; j++) {
            int cell = cellOrder[j];
            // Only run the model dynamics if the cell has habitat (cellList[cell][1] == 1)
            if(cellList[cell][1] == 1) {
                kill(cellPopInd, probDeath, cell, numSpec, eng);
                reproduction(cellPopInd, cellList, traits, cell, numSpec, i, eng);
                dispersal(cellPopInd, distArray, probDisp, cell, numSpec, dispNum, eng); // Now set to a constant probability
                immigration(cellPopInd, probImm, cell, numSpec, immNum, eng);
            }
        }

        // Calculate richness and abundance metrics
        // And save total population, richness, number of each species in each cell etc
        // Do this t/100 times, so 100 saves total
        if((i+1)%(t/500) == 0) {
            calculateTotalPopSpec(cellPopInd, totalPopSpec);
            totalPop = 0;
            for (int j = 0; j < totalPopSpec[0].size(); j++){totalPop += totalPopSpec[1][j];}
            totalRich = totalPopSpec[0].size();
            calculateCellPop(cellPopInd, cellPop);

            //After each generation store all outputs
            store2ColFiles(s_totalPop, i, totalPop);
            store2ColFiles(s_totalRich, i, totalRich);
            // Store number of times p_R was within 1% of pMax, and then reset to 0
            store2ColFiles(s_pMax, i, pMax); pMax = 0;
            storeVec(s_totalPopSpec, totalPopSpec, i, 2);
            storeVec(s_cellPop, cellPop, i, 2);
            // storecellPopInd(s_cellPopInd, cellPopInd, i);
            storeCellPopSpec(s_cellPopSpec, cellPopSpec, i, traits);

            // Store consumption rate and probability of reproduction
            // Do this less frequently (t/10) times which saves it 10 times total
            // THis is because they create huge files
            if((i+1)%(t/20) == 0) {
                storeConsumptionRate(s_consumptionRate,  cellPopInd, i, traits);
                storeReproduction(s_reproduction, cellPopInd, i, traits);
            }

            //Live output to console
            std::cout << "Time Step: " << i + 1 << "/" << t << " | Total Pop: " << totalPop << " | Total Richness: " << totalRich << "\n";

        }
    }
    // // Store final output files
    storeNum(immNum, "/total_Immigrations.txt", finpath);
    storeNum(dispNum, "/total_Dispersals.txt", finpath);
    storeNum(maxDisp, "/maximum_Dispersal.txt", finpath);

    // Close all streams to files
    s_totalPop.close(); s_totalRich.close(); s_cellPop.close(); s_totalPopSpec.close(); s_cellPopInd.close(); s_cellPopSpec.close();
    s_consumptionRate.close(); s_reproduction.close(); s_pMax.close();

    // end timer
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

    file << "Species_Pool" << " S " << numSpec << "\n";
    file << "PP_Gain_Multiplier" << " r0 " << r0 << "\n";
    file << "Carrying_Capactiy" << " K0 " << K0 << "\n";
    file << "Interference" << " I0 " << I0 << "\n";
    file << "Generation_Constant" << " G0 " << G0 << "\n";
    file << "Dispersal_Constant" << " D0 " << D0 << "\n";
    file << "pOff_Slope" << " alpha " << alpha << "\n";
    file << "pOff_Offset" << " H0 " << H0 << "\n";
    file << "Temperature" << " T " << T << "\n";
    file << "Climate_Change" << " climateChange " << climateChange << "\n";
    if(climateChange == 1) {
        file << "Climate_Change_Time" << " climateT " << climateT << "\n";
        file << "Climate_Change_Amount" << " changeT " << changeT << "\n";
    }
    file << "habitatLoss" << " loss " << habitatLoss << "\n";
    if(habitatLoss == 1) {
        file << "Habitat_Loss_Time" << " lossT " << lossT << "\n";
        file << "Habitat_Loss_Prop" << " propLoss " << propLoss << "\n";
    }

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

vector <int> findValidCells(double (&distArray)[numCells][numCells], int (&cellList)[numCells][2], double distance, int cell) {

    vector<int>validCells;

    // Go through row of distArray for the cell and identify all cells within distance
    for (int i = 0; i < numCells; i++) {
        // Check if cell is within distance and has habitat (cellList[i][1] == 1)
        if(distArray[cell][i] <= distance && cellList[i][1] == 1) {
            validCells.push_back(i);
        }
    }

    return validCells;
}

void dispersal(vector <double> (&cellPopInd)[numCells][3], double (&distArray)[numCells][numCells], 
double prob, int cell, int numSpec, int &dispNum, mt19937& eng) {

    int pop = cellPopInd[cell][0].size();

    if(pop > 0) { // Check cell isn't empty
        int chosenSpec, chosenIndex, dispCell;
        double disp; // Store distance individual is dispersing
        vector<int> validCells; // store cells that individual could disperse to

        if (uniform(eng) <= prob) {
            // Get index of random individual
            chosenIndex = randomIndex(cellPopInd, pop, numSpec, cell, eng);
            chosenSpec = cellPopInd[cell][0][chosenIndex];

            // Calculate distance individual is dispersing
            disp = dispersalDist(cellPopInd[cell][1][chosenIndex], E, eng);

            // Check if this is the largest dispersal distance
            // for storing
            if(disp > maxDisp) {
                maxDisp = disp;
            }
            

            if(disp >= 1) {  // If disp less than 1 then the individual can't disperse
                validCells = findValidCells(distArray, cellList, disp, cell);
                int arr[validCells.size()]; 
                copy(validCells.begin(), validCells.end(), arr); // Copy into array for use in shuffle function
                shuffle(arr, validCells.size(), eng);                 // Shuffle valid cells into random order
                if(arr[0] != cell) {
                    dispCell = arr[0];                             // Choose first cell in shuffled order which isn't the original cell
                } else {
                    dispCell = arr[1];
                }

                // Remove individual from cell it was in
                removeInd(cellPopInd, cell, chosenIndex, cellPopSpec);
                // Add individual to dispCell
                addInd(cellPopInd, dispCell, chosenSpec, traits, cellPopSpec);
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

double calculateInteractions(vector <double> (&cellPopInd)[numCells][3], double (&traits)[numSpec][2], int cell, int numSpec, int ind,
    vector <double> (&cellPopSpec)[numCells][3], int gen) {

    double H = 0; // Store energy state after interactions
    double CE = 0.5; // Conversion efficiency
    int Si = cellPopInd[cell][0][ind]; // Store identity of the focal species
    double Mi = cellPopInd[cell][1][ind]; // Mass of chosen individual
    double Mj; // Store mass of consumer/other species, which will be calculated from the average mass of all individuals of that species
    // in the cell (cellPopSpec[cell][3][spec]/cellPopSpec[cell][2][spec])
    int Ni; // Store number of individuals of same species as focal individual in cell
    int Sj; // Store identities of non-focal species
    int Nj; // Store abundance of non-focal species in the cell
    double Ii = 0; // Interference term for focal individual i

    double Xi; // Biomass of population of focal species i in cell C (primary producers onlu)
    double Ki; // Carrying capacity of population of focal individual i (primary prodiucers onlu)
    
    for (int k = 0; k < cellPopSpec[cell][0].size(); k++) {
        if(Si == cellPopSpec[cell][0][k]) {
            Ni = cellPopSpec[cell][1][k]; 
            break;
        }
    }
    
    // Check if focal species is a primary producer or not (1 = primary producer)
    if(cellPopInd[cell][2][ind] == 0) {
        // Loop over consumption rate for our focal species on 
        // resources (all species) in the same cell
        #pragma omp parallel for private(Sj, Nj, Mj) reduction(+:H)
        for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
            Sj = cellPopSpec[cell][0][i];
            Nj = cellPopSpec[cell][1][i];
            Mj = cellPopSpec[cell][2][i]/cellPopSpec[cell][1][i]; // Divide total biomass by number of individuals
            // If the focal individual is the same species as the chosen individual
            // then don't calculate interactions (no cannibalism)
            if(Si != Sj) {
                H += CE*Nj*consumptionRate(Mi, Mj, E, traits, Nj);
            }
        }
        // Calculate interference of focal individual i with conspecifics (only applicable for non-primary producers)
        Ii = Ni*searchRate(Mi, Mi, E, traits);
        H -= I0*Ii*Mi;       
    } else if (cellPopInd[cell][2][ind] == 1){
        Xi = getSpeciesCellMass(cell, Si, cellPopInd); // Mass of individuals in the cell of the same species
        Ki = K0*pow(Mi, 0.25); // pow(individualsMass, -0.75) * individualsMass, same as individualsMass^0.25
        // Calculate energy intake of ourprimary producer as base energy rate (r0)*mass*density function including species specific carrying capacity
        H += r0*Mi*(1-(Xi/Ki));
    }

    // Loop over consumption of focal species
    #pragma omp parallel for private(Sj, Nj, Mj) reduction(+:H)
    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        Sj = cellPopSpec[cell][0][i];
        // Check if primary producer as primary producers can't consume other species
        // Also check if it's it's the same species as the focal species
        // as species are not cannibalistic
        if(traits[Sj][1] == 0 && Si != Sj) {
            Nj = cellPopSpec[cell][1][i];
            // Calculate average mass of Mj in the cell
            // If we want a truly individual based model then we should
            // calculate the consumption rate of each individual feeding on our focal species (see addConsumptionRate branch)
            Mj = cellPopSpec[cell][2][i]/cellPopSpec[cell][1][i];
            H -= Nj*consumptionRate(Mj, Mi, E, traits, Ni);
        }
    }

    return H;

}

void reproduction(vector <double> (&cellPopInd)[numCells][3], 
int (&cellList)[numCells][2], double (&traits)[numSpec][2], int cell, int numSpec, int gen, mt19937& eng) {

    int pop = cellPopInd[cell][0].size();
    double cellMass; double B0 = 1;

    if(pop > 0) {
        int chosenSpec, chosenIndex;
        double H, pOff, G;
        
        chosenIndex = randomIndex(cellPopInd, pop, numSpec, cell, eng);
        chosenSpec = cellPopInd[cell][0][chosenIndex];
        H = calculateInteractions(cellPopInd, traits, cell, numSpec, chosenIndex, cellPopSpec, gen) - 
            (B0*pow(cellPopInd[cell][1][chosenIndex], 0.75));
        // Divide H (energy state) by the mass of the invidiaul to make the energy relative to the mass of the individual
        G = G0*pow(cellPopInd[cell][1][chosenIndex], 0.25)*arrhenius(-E);
        H = H/cellPopInd[cell][1][chosenIndex];
        pOff = (1/G)*(1/(1 + exp(-alpha*(H - H0))));

        // If pOff is within 1% of the maximum pOff (1/G), then warn user at the end
        if(pOff > (1/G) - (0.01*(1/G))) {
            pMax += 1;
        }

        if (uniform(eng) <= pOff) {
            int mutSpec = mutation(cellPopInd, probMut, chosenSpec, eng);
            addInd(cellPopInd, cell, mutSpec, traits, cellPopSpec);
        }
    }
}

void kill(vector <double> (&cellPopInd)[numCells][3], double prob, int cell, int numSpec, mt19937& eng) {
    
    int pop = cellPopInd[cell][0].size();

    if(pop > 0) { 

    int chosenIndex;
    double G;

    chosenIndex = randomIndex(cellPopInd, pop, numSpec, cell, eng);
    G = G0*pow(cellPopInd[cell][1][chosenIndex], 0.25)*arrhenius(-E);

    prob = prob/G;

        if (uniform(eng) <= prob) {
            removeInd(cellPopInd, cell, chosenIndex, cellPopSpec);
        }
    }
}

// Choose a random individual (proportional to how many of that species there are)
// NOTE: THIS NOW RETURNS THE INDEX OF THE INDIVIDUAL IN cellPopInd
int randomIndex(vector <double> (&cellPopInd)[numCells][3], int pop, int numSpec, int cell, mt19937& eng) {
    
    int sum, threshold;
    sum = 0;

    // Create threshold value by multiplying the total number of individuals in the cell by a number between 0-1
    threshold = uniform(eng)*pop;
    
    // Return threshold value which is the index of cellPopInd
    // This works as pop is the size of cellPopInd, but cellPopInds
    // addresses got to pop-1, so when you cast threshold to an int
    // it can be at most pop-1, since the numbers after the decimal are removed
    return threshold;
}

void shuffle(int arr[], int arrElements, mt19937& eng) {
    
    int randomNum;

    for (int i = 0; i < arrElements; i++){
        randomNum = chooseInRange(0, arrElements-1, eng);
        swap(arr[i], arr[randomNum]);
    }
}

void immigration(vector <double> (&cellPopInd)[numCells][3], double prob, int cell, int numSpec, int &immNum, mt19937& eng) {
    
    int chosenSpec;

    if(uniform(eng) <= prob) {  
        chosenSpec = chooseInRange(0, numSpec-1, eng);
        addInd(cellPopInd, cell, chosenSpec, traits, cellPopSpec);
        immNum++;
    }
}

int chooseInRange(int a, int b, mt19937& eng) {
    uniform_int_distribution<> choose(a, b); // define the range [a,b], extremes included.
    return choose(eng);
}

void initialisePop(vector <double> (&cellPopInd)[numCells][3], double (&traits)[numSpec][2], int numSpec, int numCells, mt19937& eng) {

    int chosenSpec;

    for (int i = 0; i < numCells; i++){
        for (int j = 0; j < initPop; j++) {
            chosenSpec = chooseInRange(0, numSpec-1, eng); 

            addInd(cellPopInd, i, chosenSpec, traits, cellPopSpec);
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

void createTraits(double (&traits)[numSpec][2], mt19937& eng, double ppProp) {
    std::lognormal_distribution<double> distribution(0, 2.0);
    for (int i = 0; i < numSpec; i++) {
        traits[i][0] = distribution(eng);
        if(uniform(eng) < ppProp) {traits[i][1] = 1;} else {traits[i][1] = 0;};
    }
}

void removeInd(vector <double> (&cellPopInd)[numCells][3], int cell, int chosenIndex, vector <double> (&cellPopSpec)[numCells][3]) {

    int Si = cellPopInd[cell][0][chosenIndex]; // Species of chosen individual

    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        if(Si == cellPopSpec[cell][0][i] && cellPopSpec[cell][1][i] > 1) {
            // Remove individul from cellpopspec (-1 from population)
            cellPopSpec[cell][1][i] -= 1;
            // Subtract weight of individual from total biomass of species
            cellPopSpec[cell][2][i] -= cellPopInd[cell][1][chosenIndex];
            break;
        } else if (Si == cellPopSpec[cell][0][i] && cellPopSpec[cell][1][i] == 1) {
            cellPopSpec[cell][0].erase(cellPopSpec[cell][0].begin() + i);
            cellPopSpec[cell][1].erase(cellPopSpec[cell][1].begin() + i);
            cellPopSpec[cell][2].erase(cellPopSpec[cell][2].begin() + i);
        }
    }

    for (int j = 0; j < 3; j++) {cellPopInd[cell][j].erase(cellPopInd[cell][j].begin() + chosenIndex);}

}

// Function to add a species which may not already exist in the cell - immigration, mutation, initialisation
void addInd(vector <double> (&cellPopInd)[numCells][3], int cell, int chosenSpec, double (&traits)[numSpec][2], vector <double> (&cellPopSpec)[numCells][3]) {

    cellPopInd[cell][0].push_back(chosenSpec);
    cellPopInd[cell][1].push_back(traits[chosenSpec][0]);
    cellPopInd[cell][2].push_back(traits[chosenSpec][1]);

    bool exists = false;
    // Parallelise below loop
    // Slower format when cellPopSpec is small (equal at about 1000 species)
    // But quicker when it's large
    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        if(chosenSpec == cellPopSpec[cell][0][i]) {
            exists = true;
            cellPopSpec[cell][1][i] += 1;
            cellPopSpec[cell][2][i] += traits[chosenSpec][0];
            break;
        }
    }

    if(exists == false) {
        cellPopSpec[cell][0].push_back(chosenSpec);
        cellPopSpec[cell][1].push_back(1);
        cellPopSpec[cell][2].push_back(traits[chosenSpec][0]);
    }
    

}

// Function for when a species already exists and is being replicated (maybe with intraspecific variation)
void baby(vector <double> (&cellPopInd)[numCells][3], int cell, int chosenIndex, double (&traits)[numSpec][2]) {

    int chosenSpec = cellPopInd[cell][0][chosenIndex];

    cellPopInd[cell][0].push_back(chosenSpec);
    cellPopInd[cell][1].push_back(traits[chosenSpec][0]);
    cellPopInd[cell][2].push_back(traits[chosenSpec][1]);
}

void calculateTotalPopSpec(vector <double> (&cellPopInd)[numCells][3], vector <int> (&totalPopSpec)[2]) {

    totalPopSpec[0].clear(); totalPopSpec[1].clear();
    
    for (int i = 0; i < numCells; i++) {
        for (int j = 0; j < cellPopSpec[i][0].size(); j++) {
            bool exists = false;
            // Parallelise inner loop
            // As outer loop has push back which I don't know how to deal with
            // when parallelising at the moment
            for (int k = 0; k < totalPopSpec[0].size(); k++) {
                if(cellPopSpec[i][0][j] == totalPopSpec[0][k]) {
                    totalPopSpec[1][k] += cellPopSpec[i][1][j];
                    exists = true;
                }
            }
            if(exists == false) {
                totalPopSpec[0].push_back(cellPopSpec[i][0][j]);
                totalPopSpec[1].push_back(cellPopSpec[i][1][j]);
            }
        }
    }
}

void calculateCellPop(vector <double> (&cellPopInd)[numCells][3], vector <int> (&cellPop)[2]) {

    cellPop[0].clear(); cellPop[1].clear();

    int pop = 0;

    for (int i = 0; i < numCells; i++) {
        pop = cellPopInd[i][0].size();
        if(pop > 0) {
            cellPop[0].push_back(i);
            cellPop[1].push_back(pop);
        }
    }
}

void storecellPopInd(ofstream &stream, vector <double> (&cellPopInd)[numCells][3], int gen) {

    for (int i = 0; i < numCells; i++) {
        for (int j = 0; j < cellPopInd[i][0].size(); j++) {
            stream << gen+1 << " " << i+1 << " " << cellPopInd[i][0][j] + 1 << " " << cellPopInd[i][1][j] << " " << 
            cellPopInd[i][2][j] << "\n";
        }
    }
}

void storeCellPopSpec(ofstream &stream, vector <double> (&cellPopSpec)[numCells][3], int gen, double (&traits)[numSpec][2]) {

    for (int i = 0; i < numCells; i++) {
        double cellMass = getCellMass(i, cellPopInd);
        for (int j = 0; j < cellPopSpec[i][0].size(); j++) {
            int Si = cellPopSpec[i][0][j];
            stream << gen+1 << " " << i+1 << " " << cellPopSpec[i][0][j] + 1 << " " << cellPopSpec[i][1][j] << " " << 
            cellPopSpec[i][2][j] << " " << traits[Si][0] << " " << traits[Si][1] << "\n";
        }
    }
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

int mutation(vector <double> (&cellPopInd)[numCells][3], double prob, int chosenSpec, mt19937& eng) {
    
    int b1[L]; // To store the binary sequence
    int mutSpec; // Integer identifier of mutated species that will be created

    ConvertToBinary(chosenSpec, L, b1);

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

void loseHabitat(int (&landscapeArray)[cellCols][cellRows], int (&cellList)[numCells][2], int cellRows, int cellCols, int numCells, mt19937& eng, double propLoss) {

    int numCellsLoss = numCells*propLoss; // Number of cells to be lost

    int lost = 0; // Counts number of cells lost

    // Update cellList to indicate which cells have lost habitat
    while(lost < numCellsLoss) {
        int cell = chooseInRange(0, numCells-1, eng); // Choose a random cell to start the loss from
        // Change second column of cellList from 1 to 0 to indicate they no longer have habitat
        if(cellList[cell][1] == 1) {
            cellList[cell][1] = 0;
            lost++; // Only add one to lost if the cell was covered in habitat before
        }
    }

    // Update landscapeArray to indicate which cells have lost habitat
    for (int i = 0; i < numCells; i++) {
        if(cellList[i][1] == 0) {
            int col = cellList[i][0] / cellCols;
            int row = cellList[i][0] % cellRows;
            landscapeArray[col][row] = 0;
        }
    }
    
    

    // Then remove all individuals from cells that are no longer covered in habitat
    for (int i = 0; i < numCells; i++) {
        if(cellList[i][1] == 0) {
            for (int j = 0; j < 3; j++) {
                cellPopInd[i][j].clear();
                cellPopSpec[i][j].clear();
            }
            
        }
    }
    


}

///////////////////////
// Metabolic Theory Functions
///////////////////////

double searchRate(double Mi, double Mj, double E, double (&traits)[numSpec][2]) {
    
    double V0 = 0.33; double d0 = 1.62; 
    double a;

    a = 2*V0*d0*(pow(Mi, 0.63))*arrhenius(E);

    return a;
    
}

double attackProb(double Mi, double Mj, double (&traits)[numSpec][2]) {

    double Rp = 0.1;

    double A = (1/(1 + 0.25*(pow(exp(1), -pow(Mi, 0.33)))))*pow(1/(1 + pow(log10(Rp*(Mi/Mj)), 2)), 5);

    return A;

}

double handlingTime(double Mi, double Mj, double E, double (&traits)[numSpec][2]) {

    double h0 = 1;
    double Rp = 0.1;

    double h = h0*std::pow(Mi, -0.75)*(1-exp(-(std::pow((Mj/Mi) - Rp, 2))/2))*arrhenius(E);

    return h;

}

double consumptionRate(double Mi, double Mj, double E, double (&traits)[numSpec][2], int Nj) {

    double c = (searchRate(Mi, Mj, E, traits)*attackProb(Mi, Mj, traits)*Mj)/
    (1 + (searchRate(Mi, Mj, E, traits)*attackProb(Mi, Mj, traits)*handlingTime(Mi, Mj, -E, traits)*Nj));

    return c;

}

double arrhenius(double E) {

    return pow(exp(1), (E*((T + T0) - (T0+20)))/(k*(T+T0)*(T0+20)));

}

double getCellMass(int cell, vector <double> (&cellPopInd)[numCells][3]) {

    double cellMass = 0;

    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        cellMass += cellPopSpec[cell][2][i]; 
    }

    return cellMass;

}

double getSpeciesCellMass(int cell, int chosenSpec, vector <double> (&cellPopInd)[numCells][3]) {
    
    double speciesCellMass = 0;

    // Parallelize this loop
    // Slower when cellPopSpec is small
    // Should be quicker as cellPopSpec gets larger though
    for (int i = 0; i < cellPopSpec[cell][0].size(); i++) {
        if(cellPopSpec[cell][0][i] == chosenSpec) {
            speciesCellMass = cellPopSpec[cell][2][i];
            break;
        }
    }
    
    return speciesCellMass;

}

void storeConsumptionRate(ofstream &stream, vector <double> (&cellPopInd)[numCells][3], int gen, double (&traits)[numSpec][2]) {

    // This is currently calculated at the species level as it was very slow at
    // the individual level, even though we may want to store at individual level
    // once we have intraspecific variation

    for (int i = 0; i < numCells; i++) {
        for (int j = 0; j < cellPopSpec[i][0].size(); j++) {
            int Si = cellPopSpec[i][0][j];         
            // Check if the species is a primary producer (no consumption)
            if(traits[Si][1] == 0) {
                for (int k = 0; k < cellPopSpec[i][0].size(); k++) {
                    int Sj = cellPopSpec[i][0][k];
                    if(Si != Sj) {
                        stream << gen+1 << " " << i+1 << " " << Si + 1 << " " << traits[Si][0] << " " << 
                        cellPopSpec[i][1][j] << " " << (Sj + 1) << " " << traits[Sj][0] << 
                        " " << cellPopSpec[i][1][k] << " " <<
                        consumptionRate(Si, Sj, E, traits, cellPopSpec[i][1][k]) << "\n";
                    }
                }
            } else {
                for (int k = 0; k < cellPopSpec[i][0].size(); k++) {
                    int Sj = cellPopSpec[i][0][k];
                    if(Si != Sj) {
                        stream << gen+1 << " " << i+1 << " " << Si + 1 << " " << traits[Si][0] << " " << 
                        cellPopSpec[i][1][j] << " " << (Sj + 1) << " " << traits[Sj][0] << 
                        " " << cellPopSpec[i][1][k] << " " << 0 << "\n";
                    }
                }
            }
        }
    }

}

// THIS NEEDS TO BE HEAVILY MODIFIED WHEN WE CONVERT TO AN INTRASPECIFIC VARIATION MODEL
void storeReproduction(ofstream &stream, vector <double> (&cellPopInd)[numCells][3], int gen, double (&traits)[numSpec][2]) {

    // This is currently calculated at the species level as it was very slow at
    // the individual level, even though we may want to store at individual level
    // once we have intraspecific variation

    for (int i = 0; i < numCells; i++) {
        for (int j = 0; j < cellPopSpec[i][0].size(); j++) {

        double cellMass; double B0 = 1;
        // double B0 = 4.15*pow(10, -8);

        int chosenSpec, chosenIndex;
        double H, pOff, G;
        
        cellMass = getCellMass(i, cellPopInd); // Get the total mass of all individuals in the cell
        chosenSpec = cellPopSpec[i][0][j];
        // Now we need to cheat and find the index of an individual of the species in cellPopInd
        // later I need to overhaul this approach
        for (int k = 0; k < cellPopInd[i][0].size(); k++) {if(cellPopInd[i][0][k] == chosenSpec){chosenIndex = k;}}
        H = calculateInteractions(cellPopInd, traits, i, numSpec, chosenIndex, cellPopSpec, gen) - 
            (B0*std::pow(traits[chosenSpec][0], 0.75));
        // Divide H (energy state) by the mass of the invidiaul to make the energy relative to the mass of the individual
        G = G0*pow(traits[chosenSpec][0], 0.25)*arrhenius(-E);
        H = H/traits[chosenSpec][0];
        pOff = (1/G)*(1/(1 + exp(-alpha*(H - H0))));

        stream << gen + 1 << " " << i+1 << " " << chosenSpec+1 << " " << traits[chosenSpec][0] << " " << H << " " << 1/G << " " << 
        pOff << " " << probDeath/(G0*pow(traits[chosenSpec][0], 0.25)*arrhenius(-E)) << "\n";

        }
    }
}

double dispersalDist(double Mi, double E, mt19937& eng) {

    std::exponential_distribution<double> expDist(1.5);

    double disp = D0*expDist(eng)*pow(Mi, 0.63)*arrhenius(E);

    return disp;

}