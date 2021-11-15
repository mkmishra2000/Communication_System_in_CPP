// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- BPSK MODEL AND "PROBABILITY OF SYMBOL ERROR" vs "RATIO OF SIGNAL TO NOISE ENERGY (dB)"
// Starting date:- 29-SEP-2021
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>

#define one_million 1000000 

using namespace std;


// Function for generating binary bits at source side. Each bit is equiprobable.
// Input is nothing
// Output is a vector that contains the binary bits of length one_million*1.
vector<double> sourceVector()
{
    vector<double> sourceBits;

    // Use current time as seed for random generator
    srand(time(0));
 
    for(int i = 0; i<one_million; i++){
        sourceBits.insert(sourceBits.end(), rand()%2);
    }

    return sourceBits;
}


// Function for mapping bits to symbol. 
// Input is a binary bit vector. Here 0---> -(sqrt(Energy)) and 1---> (sqrt(Energy))
// Output is a vector that contains transmitted symbols.
vector<double> bit_maps_to_symbol_of_energy_E(vector<double> sourceBits, double energyOfSymbol)
{
    vector<double> transmittedSymbol;

    for(int i=0; i<one_million; i++){
        if(sourceBits[i]== 0){
            transmittedSymbol.insert(transmittedSymbol.end(), -sqrt(energyOfSymbol));
        }
        else{
            transmittedSymbol.insert(transmittedSymbol.end(), sqrt(energyOfSymbol));
        }

    }

    return transmittedSymbol;
}


// Function for generating random noise based on gaussian distribution N(mean, variance).
// Input mean and standard deviation.
// Output is the vector that contain gaussian noise as an element.
vector<double> GnoiseVector(double mean, double stddev)
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



// Function for modeling additive channel. Here gaussian noise adds to the transmitted bit.
// Inputs are the transmitted bit and gaussian noise with mean 0 and variance 1.
// Output is the receive bits.
vector<double> receiveBits(vector<double> transBit, vector<double> gnoise)
{
    vector<double> recievebits;

    for(int j =0; j<transBit.size(); j++){
        recievebits.insert(recievebits.end(), transBit[j]+gnoise[j]);
    }

    return recievebits;

}


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

// Function to count number of errors in the received bits.
// Inputs are the sourcebits and decodedbits
// OUtput is the number of error in received bits.
// error: if sourcebit  != receivebit
int errorCalculation (vector<double> sourceBits, vector<double> decodedBits)
{
    int countError =0;
    for(int i =0; i<sourceBits.size();i++){
        if(sourceBits[i]!= decodedBits[i]){
            countError++;
        }
    }

    return countError;
}


int main(){

    // source defination
    vector<double> sourceBits;

    // Mapping of bits to symbols;
    vector<double> transmittedSymbol;

    // Noise definition
    vector<double> gnoise;

    // Receive bits
    vector<double> recevBIts;

    // Decision block
    vector<double> decodedBits;

    sourceBits = sourceVector();
    // std::copy(begin(sourceBits), end(sourceBits),  std::ostream_iterator<double>(std::cout, " "));

    double Energy_of_symbol =1; 

    transmittedSymbol = bit_maps_to_symbol_of_energy_E(sourceBits, Energy_of_symbol);

    // std::copy(begin(transmittedSymbol), end(transmittedSymbol),  std::ostream_iterator<double>(std::cout, " "));

    gnoise = GnoiseVector(0.0, 1);

    recevBIts = receiveBits(transmittedSymbol, gnoise);

    decodedBits = decisionBlock(recevBIts);

    int countErrors = errorCalculation(sourceBits, decodedBits);
    cout<< endl;

    cout<< countErrors << endl;
    cout<< endl;
    // std::copy(begin(recevBIts), end(recevBIts),  std::ostream_iterator<double>(std::cout, " "));



    return 0;
}