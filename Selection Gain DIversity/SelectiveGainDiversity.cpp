/////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- Performance of selective gain combiner using BPSK in Rayleigh fadding channel
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>
#include <fstream>
#include "BPSK.h"
#include "MyMatrixOperation.h"

#define one_million 10

using namespace std;



// Function for printing the vector on the console output.
void PrintVectorDouble(vector<double> vectr)
{
    std::copy(begin(vectr), end(vectr),  std::ostream_iterator<double>(std::cout, "   "));
    cout<<endl;
}


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



// Function for Rayleigh fadding cofficients
// Inputs are two vectors one is the gaussion noise as real part and second as Gaussian noise as imazinary part 
// Output is a vector that contain rayliegh noise coff, sqrt(real_part^2 + imz_part^2)
vector<double> RayleighFaddingCoff(vector<double> realGaussian, 
                                   vector<double> ImziGaussian)
{
    vector<double> rayleighNoise;
    double temp;

    for(int times=0; times<realGaussian.size(); times++){

        temp = sqrt(pow(realGaussian[times], 2)+pow(ImziGaussian[times], 2));
        rayleighNoise.insert(rayleighNoise.end(), temp);
    }

    return rayleighNoise;
}



vector<double> MaximumSNRIndex(vector<vector<double>> SNRMat)
{
    vector<double> MaxIndexes;
    double index;

    for(int i =0; i<SNRMat[0].size(); i++){
        index =0;

        for(int j =0; j<SNRMat.size(); j++){

            if(SNRMat[index][i]<SNRMat[j][i]){
                index = j;
            }
        }

        MaxIndexes.insert(MaxIndexes.end(), index);
    }

    return MaxIndexes;
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


int main(){

    // source defination
    vector<double> sourceBits;
    sourceBits = sourceVector();

    // Mapping of bits to symbols;
    vector<double> transmittedSymbol;

    // Noise definition
    vector<double> gnoise;

    vector<double> realGaussian;
    vector<double> imziGaussian;

    vector<double> RayleighNoise;
    vector<double> MaxiSNRdetection;
    vector<double> decodedBits;

    double energy =64;
    double sigmaSquare = 0.5;
    double stddevRayleigh = sqrt(sigmaSquare);
    double N_o =4;
    double p, stdnoise;
    stdnoise = sqrt(N_o);
    double counterror, P_error;

    transmittedSymbol = bit_maps_to_symbol_of_energy_E(sourceBits, energy, one_million);

    int L =2;

    vector<vector<double>> MultipleSignals;
    vector<vector<double>> RaleighMat;
    vector<vector<double>> SNRMat;
    vector<vector<double>> GnoiseMat;
    vector<vector<double>> NewMat;
    vector<vector<double>> ReceiveMat;

    for(int i =0; i<L; i++){
        MultipleSignals.insert(MultipleSignals.end(), transmittedSymbol);
    }

    for(int i =0; i<L; i++){
        realGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);
        imziGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);

        RayleighNoise = RayleighFaddingCoff(realGaussian, imziGaussian);

        RaleighMat.insert(RaleighMat.end(), RayleighNoise);

        realGaussian.clear();
        imziGaussian.clear();
        RayleighNoise.clear();
    }

    for(int i=0; i<L; i++){
        gnoise = GnoiseVector(0.0, stdnoise, one_million);
        // PrintVectorDouble(gnoise);
        GnoiseMat.insert(GnoiseMat.end(), gnoise);
        gnoise.clear();
    }

    NewMat = ElementWiseMultiplication(MultipleSignals, RaleighMat);
    // sourceBits = sourceVector();

    PrintVectorDouble(sourceBits);

    // cout<<"Transmission"<<endl;
    // PrintMat(MultipleSignals);

    // cout<<"Rayleigh noise"<<endl;
    // PrintMat(RaleighMat);

    // cout<<"Rayleigh operation"<<endl;
    // PrintMat(NewMat);

    // cout<<"Gnoise matrix"<<endl;
    // PrintMat(GnoiseMat);

    // cout<<"Received signals"<<endl;
    ReceiveMat = ElementWiseAddition(NewMat, GnoiseMat);
    // PrintMat(NewMat2);

    // cout<<"SNR coff"<<endl;
    SNRMat = ElementWiseMultiplication(RaleighMat, RaleighMat);
    // PrintMat(SNRMat);

    MaxiSNRdetection = MaximumSNRIndex(SNRMat);
    // cout<<"Maximum SNR index"<<endl;
    // PrintVectorDouble(MaxiSNRdetection);

    vector<double> combinerOutput;

    for(int j =0; j<MaxiSNRdetection.size(); j++){
        combinerOutput.insert(combinerOutput.end(), ReceiveMat[MaxiSNRdetection[j]][j]);

    }

    // PrintVectorDouble(combinerOutput);

    decodedBits = decisionBlock(combinerOutput);
    cout<<"decoded output"<<endl;
    PrintVectorDouble(decodedBits);
    cout<<endl;

    counterror = errorCalculation(sourceBits, decodedBits);
    P_error = counterror/one_million;

    cout<<"Error count: "<< counterror<<endl;
    cout<<"P_e:         "<< P_error<<endl;

    return 0;
}