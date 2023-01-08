/////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- BFSK MODEL AND "PROBABILITY OF Bit ERROR" vs "RATIO OF SIGNAL TO NOISE ENERGY (dB)"
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>
#include <fstream>

// #define one_million 1000000

using namespace std;



// Function of generate vector form of signal, this function generate the upper value of the vector.
// Input is the source bit stream and energy of bit.
// Output is the vector containing all upper values based on the source bits. 
vector<double>TransmittedSymbol_columnVector1(vector<double> sourceBits, double& energyOfbit, const int one_million)
{
    vector <double> transmittedSymbol;

     for(int i=0; i<one_million; i++){
        if(sourceBits[i]== 0){
            transmittedSymbol.insert(transmittedSymbol.end(), sqrt(energyOfbit));
        }
        else{
            transmittedSymbol.insert(transmittedSymbol.end(), 0);
        }

    }

    return transmittedSymbol;
}



// Function of generate vector form of signal, this function generate the lower value of the vector.
// Input is the source bit stream and energy of bit.
// Output is the vector containing all lower values based on the source bits. 
vector<double>TransmittedSymbol_columnVector2(vector<double> sourceBits, double& energyOfbit, const int one_million)
{
    vector <double> transmittedSymbol;

     for(int i=0; i<one_million; i++){
        if(sourceBits[i]== 0){
            transmittedSymbol.insert(transmittedSymbol.end(), 0);
        }
        else{
            transmittedSymbol.insert(transmittedSymbol.end(), sqrt(energyOfbit));
        }

    }

    return transmittedSymbol;
}



// Function for generating random noise based on gaussian distribution N(mean, variance).
// Input mean and standard deviation.
// Output is the vector that contain gaussian noise as an element.
vector<double>GnoiseVector(double& mean, double& stddev, const int one_million)
{
    vector<double> noise;
    
    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    std::normal_distribution<double> dist(mean, stddev);

    // Add Gaussian noise
    for (int i =0; i<=one_million; i++) {
        noise.insert(noise.end(),dist(generator));
    }

    return noise;
}



// Function for take decision after receiving the vector
// Input is the 2D receiver vector
// Output is the decision {0, 1} 
vector <double>DecisionBlock(vector<double> receive1, vector<double> receive2)
{
    vector<double> decision;

    for(int j =0; j<receive1.size(); j++){

        if(receive1[j]< receive2[j]){
            decision.insert(decision.end(), 1);
        }else{
            decision.insert(decision.end(), 0);
        }
    }

    return decision;
}


