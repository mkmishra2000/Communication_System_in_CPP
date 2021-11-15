// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- BPSK MODEL AND "PROBABILITY OF SYMBOL ERROR" vs "RATIO OF SIGNAL TO NOISE ENERGY (dB)"
// Starting date:- 29-SEP-2021
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>
#include <fstream>

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
double errorCalculation (vector<double> sourceBits, vector<double> decodedBits)
{
    double countError =0;
    for(int i =0; i<sourceBits.size();i++){
        if(sourceBits[i]!= decodedBits[i]){
            countError++;
        }
    }

    return countError;
}


// Function to store the data in the file (.dat)
// Input is the SNR per bit in dB and calculated probability of error
// Output is the nothing but in processing it is creating a file and writing data into it.
void datafile(vector<double> xindB, vector<double> Prob_error)
{
    ofstream outfile;

    outfile.open("outputfile/Filedata_E_changing.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<xindB.size(); i++){
        outfile<< xindB[i] << " "<<"\t"<<" "<< Prob_error[i]<< endl;
    }

    outfile.close();
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


    vector<double> x_indB;
    for(float i =0.5; i<=12; i=i+0.07)
    {
        x_indB.insert(x_indB.end(), i);
    }

    // N_o corresponds to the variance of noise
    vector <double> N_o; // can be any positive real number. 
    // N_o= 3, 4, 5, 6, 7, 8, 9, 10 my experiment.

    double p, stdnoise;

    double energyOfSymbol = 9;

    vector<double> Prob_error;

    for(int i =0; i<x_indB.size(); i++){
        p = x_indB[i]/10;
        N_o.insert(N_o.end(), energyOfSymbol/pow(10, p));
    }


    for(int step =0; step<N_o.size(); step++){
        sourceBits = sourceVector();
   
        transmittedSymbol = bit_maps_to_symbol_of_energy_E(sourceBits, energyOfSymbol);

        stdnoise = sqrt(N_o[step]/2); // std of noise.
        gnoise = GnoiseVector(0.0, stdnoise);

        recevBIts = receiveBits(transmittedSymbol, gnoise);

        decodedBits = decisionBlock(recevBIts);

        double countErrors = errorCalculation(sourceBits, decodedBits);
        cout<<" Energy of symbol "<< energyOfSymbol<<endl;
        cout<<" Std of noise "<< stdnoise<<endl;
        cout<< endl;

        cout<< "Error count "<< countErrors << endl;
        cout<< endl;

        double pe = countErrors/one_million;

        Prob_error.insert(Prob_error.end(), pe);

    }

    datafile(x_indB, Prob_error);

    return 0;
}