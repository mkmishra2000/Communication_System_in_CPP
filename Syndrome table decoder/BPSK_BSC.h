/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organization:- IIITDM KANCHEEPURAM
// Topic:- Header file to model binary symmetric channel with BPSK modulation scheme.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <random>
#include <chrono>

using namespace std;

// Modulator block 
////////////////////////////////////////

// Function for mapping bits to symbol energy. 
// Input is a binary bit vector. Here 0---> -(sqrt(Energy)) and 1---> (sqrt(Energy))
// Output is a vector that contains transmitted symbols. Output type <double>
vector<double> bit_maps_to_symbol_of_energy_E(vector<int> sourceBits, double energyOfSymbol)
{
    vector<double> transmittedSymbol;

    for(int i=0; i<sourceBits.size(); i++){
        if(sourceBits[i]== 0){
            transmittedSymbol.insert(transmittedSymbol.end(), -sqrt(energyOfSymbol));
        }
        else{
            transmittedSymbol.insert(transmittedSymbol.end(), sqrt(energyOfSymbol));
        }

    }

    return transmittedSymbol;
}

// Noise generation (AWG noise)
///////////////////////////////////////////////////////////////////

// Function for generating random noise based on gaussian distribution N(mean, variance).
// Input mean and standard deviation.
// Output is the vector that contain gaussian noise as an element.
vector<double> GnoiseVector(double mean, double stddev, const int one_million)
{
    std::vector<double> data;
    
    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    std::normal_distribution<double> dist(mean, stddev);

    // Add Gaussian noise
    for (int i =0; i<one_million; i++) {
        data.insert(data.end(),dist(generator));
    }

    return data;
}

// AWGN channel model 
/////////////////////////////////////////////////////////////////////

// Function for modeling additive channel. Here gaussian noise adds to the transmitted bit.
// Inputs are the transmitted bits and gaussian noise with mean 0 and defined variance.
// Output is the receive bits.
vector<double> receiveBits(vector<double> transBit, vector<double> gnoise)
{
    vector<double> recievebits;

    for(int j =0; j<transBit.size(); j++){
        recievebits.insert(recievebits.end(), transBit[j]+gnoise[j]); //Additive noise.
    }

    return recievebits;

}

// Demodulator block 
////////////////////////////////////////////////////////////////

// Function for deciding the bit value from the received bits 
// Input is the received bits.
// Output is the decoded bits.
// Decision rule :- if receiveBit >0 then 1 otherwise 0 (simple Binary detection)
vector<double> decisionBlock(vector<double> receiveBits)
{
    vector<double> decodedBits;

    for(int i =0; i<receiveBits.size(); i++){
        if(receiveBits[i]>0){
            decodedBits.insert(decodedBits.end(), 1);
        }
        else{
            decodedBits.insert(decodedBits.end(), 0);
        }
    }

    return decodedBits;
}