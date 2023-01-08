/////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- Alamouti code (Space time codes)
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

#define one_million 1000000

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


// Function for multiplying fadding coff to particular antenna
// Inputs are the Transmitted energy and Rayleigh fadding coff 
// Output is the multiplication of element by element fadding and energy  
vector<double> RayleighOperation(vector<double> TxEnergy,
                                 vector<double> RayleighFaddingCoff)
{
    vector<double> Resultant;

    for(int j =0; j<TxEnergy.size(); j++){
        Resultant.insert(Resultant.end(), TxEnergy[j]*RayleighFaddingCoff[j]);
    }

    return Resultant;
}


// Function for gaussian noise for each tx bit to particular antenna
// Inputs are the Transmitted energy and Gnoise
// Output is the addition of element by element Gnoise and Tx  
vector<double> GaussianNoiseAdd(vector<double> TxEnergy,
                                 vector<double> Gnoise)
{
    vector<double> Resultant;

    for(int j =0; j<TxEnergy.size(); j++){
        Resultant.insert(Resultant.end(), TxEnergy[j]+Gnoise[j]);
    }

    return Resultant;
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
void datafile(vector<double> xindB, vector<double> Prob_error, char strName[])
{
    ofstream outfile;

    string filename = strName;

    outfile.open(filename +"."+"dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<xindB.size(); i++){
        outfile<< xindB[i] << " "<<"\t"<<" "<< Prob_error[i]<< endl;
    }

    outfile.close();
}


vector<double> CalculatedError(vector<double> SNR_dB, double L)
{
    vector<double> ProbError;

    double po, normalValue, inter;
    for (int k =0; k<SNR_dB.size(); k++){
        normalValue = pow(10, (SNR_dB[k]/10));
        po = 0.8/pow(normalValue,2);
        ProbError.insert(ProbError.end(), po);
    }

    return ProbError;
}


int main()
{
    // source defination
    vector<double> sourceBits;

    // Mapping of bits to symbols;
    vector<double> transmittedSymbol;

    // Noise definition
    vector<double> gnoise;

    //Rayleigh noise (Real guass and Img Guass)
    vector<double> realGaussian1;
    vector<double> imziGaussian1;
    vector<double> RayleighCoff1;
    vector<double> realGaussian2;
    vector<double> imziGaussian2;
    vector<double> RayleighCoff2;

    vector<double> Tx1, Tx2;
    vector<double> Rx1, Rx2;

    //Alamouti decoder
    vector<double> Dec1, Dec2;

    vector<double> decodedBits1;
    vector<double> decodedBits2;

    vector<double> mergeDecoded;

    double N_o =4;
    double p, stdnoise;
    stdnoise = sqrt(N_o);
    double counterror, P_error;

    double sigmaSquare = 0.5;
    double stddevRayleigh = sqrt(sigmaSquare);

    vector<double> SNR_dB;
    for(float i =0; i<=25; i=i+0.5)
    {
        SNR_dB.insert(SNR_dB.end(), i);
    }

    vector<double> energyOfSymbol;
    vector<double> Prob_error;
    double normalValue;

    for(int i =0; i<SNR_dB.size(); i++){

        normalValue = pow(10, (SNR_dB[i]/10));
        energyOfSymbol.insert(energyOfSymbol.end(), N_o*normalValue);
    }

    for(int step =0; step <energyOfSymbol.size(); step++){
        sourceBits = sourceVector();    //Source bit streame

        transmittedSymbol = bit_maps_to_symbol_of_energy_E(sourceBits, 
                                                            energyOfSymbol[step],
                                                            one_million);
        

        realGaussian1 = GnoiseVector(0.0, stddevRayleigh, one_million/2);
        imziGaussian1 = GnoiseVector(0.0, stddevRayleigh, one_million/2);

        RayleighCoff1 = RayleighFaddingCoff(realGaussian1, imziGaussian1);

        realGaussian2 = GnoiseVector(0.0, stddevRayleigh, one_million/2);
        imziGaussian2 = GnoiseVector(0.0, stddevRayleigh, one_million/2);

        RayleighCoff2 = RayleighFaddingCoff(realGaussian2, imziGaussian2);


        /*
        Tx1 = h1*s1+h2*s2
        Tx2 = h2*s1-h1*s2
        */
        for(int i =0; i<transmittedSymbol.size(); i=i+2 ){
            Tx1.insert(Tx1.end(), RayleighCoff1[i/2]*transmittedSymbol[i] + RayleighCoff2[i/2]*transmittedSymbol[i+1]);
            Tx2.insert(Tx2.end(),RayleighCoff2[i/2]*transmittedSymbol[i] - RayleighCoff1[i/2]*transmittedSymbol[i+1]);
        }

        gnoise = GnoiseVector(0.0, stdnoise, one_million);

        /*
        Rx1 = Tx1+gnoise;
        Rx2 = Tx2+gnoise;
        */
        for(int j =0; j<gnoise.size(); j =j+2){
            Rx1.insert(Rx1.end(), Tx1[j/2]+gnoise[j]);
            Rx2.insert(Rx2.end(), Tx2[j/2]+gnoise[j+1]);
        }

        for(int k =0; k<Tx1.size(); k++){
            Dec1.insert(Dec1.end(), RayleighCoff1[k]*Rx1[k]+RayleighCoff2[k]*Rx2[k]);
            Dec2.insert(Dec2.end(), RayleighCoff2[k]*Rx1[k]-RayleighCoff1[k]*Rx2[k]);
        }


        decodedBits1= decisionBlock(Dec1);
        decodedBits2 = decisionBlock(Dec2);

        for(int j =0; j<decodedBits1.size(); j++){
            mergeDecoded.insert(mergeDecoded.end(), decodedBits1[j]);
            mergeDecoded.insert(mergeDecoded.end(), decodedBits2[j]);
        }

        counterror = errorCalculation(sourceBits, mergeDecoded);
        P_error = counterror/one_million;
        Prob_error.insert(Prob_error.end(), P_error);

        cout<<endl;
        cout<<"Error count          :"<<counterror<<endl;
        cout<<"Probability of error :"<<P_error<<endl;
        cout<<endl;

        Tx1.clear();
        Tx2.clear();
        Rx1.clear();
        Rx2.clear();
        Dec1.clear();
        Dec2.clear();
        mergeDecoded.clear();

    }

    char NameofFile1[30] = "Alc1";
    char NameofFile2[30] = "Alcerr1";
    int L =2;

    datafile(SNR_dB, Prob_error, NameofFile1);

    vector<double> Error = CalculatedError(SNR_dB, 2);

    datafile(SNR_dB, Error, NameofFile2);

    return 0;
}