/////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- Performance of BPSK in Rayleigh fadding channel
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>
#include <fstream>
#include "BPSK.h"

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


// function for modelling wireless channel
// Inputs are the transmitted signal (X_k), rayleigh cofficient (H_k) and AWGN (N_k) component 
// Output is the H_k*X_k+N_k  
vector<double> ChannelOperation(vector<double> TransSignal, 
                                vector<double> RayleighCoff, 
                                vector<double> gnoise)
{
    vector<double> channelResult;
    double temp;
    for(int i =0; i<TransSignal.size(); i++){
        temp = (RayleighCoff[i]*TransSignal[i]) + gnoise[i];
        channelResult.insert(channelResult.end(), temp);
    } 

    return channelResult;
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

    outfile.open("fadeqchan2.dat");
    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<xindB.size(); i++){
        outfile<< xindB[i] << " "<<"\t"<<" "<< Prob_error[i]<< endl;
    }

    outfile.close();
}



// vector<double> Errorfunction(vector <double> SNR_dB)
// {
//     vector <double> Qvalue;
    
//     double po, normalValue, inter;
//     for (int k =0; k<SNR_dB.size(); k++){
//         normalValue = pow(10, (SNR_dB[k]/10));
//         // inter = 1+(2/normalValue);
//         // inter = sqrt(inter);
//         // po = 0.5*(1-(1/inter));
//         po = 1/(2*(normalValue));
//         Qvalue.insert(Qvalue.end(), po);
//     }

//     return Qvalue;
// }

// Function to store the data in the file (.dat)
// Input is the SNR per bit in dB and calculated Qfunction values
// Output is the nothing but in processing it is creating a file and writing data into it.
// void ErrorValueInFile(vector <double> SNR, vector <double> Qvalue)
// {
//     ofstream outfile;

//     outfile.open("EqualGain_Qvalue1.dat");

//     if(!outfile.is_open()){
//         cout<<"File opening error !!!"<<endl;
//         return;
//     }

//     for(int i =0; i<SNR.size(); i++){
//         outfile<< SNR[i] << " "<<"\t"<<" "<< Qvalue[i]<< endl;
//     }

//     outfile.close();
// }


vector<double> EqualGainCombining(vector<double> rec1, vector<double> rec2)
{
    vector<double> recev;
    for(int num =0; num<rec1.size(); num++){
        recev.insert(recev.end(), (rec1[num]+rec2[num]));
    }
    return recev;
}

int main()
{
    // source defination
    vector<double> sourceBits;
    
    // Mapping of bits to symbols;
    vector<double> transmittedSymbol;

    // Noise definition
    vector<double> gnoise;

    vector<double> realGaussian;
    vector<double> imziGaussian;

    vector<double> RayleighNoise;

    vector<double> receiveSignal;
    vector<double> receiveSignal1;
    vector<double> receiveSignal2;
    vector<double> receiveSignal3;

    vector<double> decodedBits;

    double sigmaSquare = 0.5;
    double stddevRayleigh = sqrt(sigmaSquare);
    double N_o =4;
    double p, stdnoise;
    double counterror, P_error;
    stdnoise = sqrt(N_o);
    
    vector<double> SNR_dB;
    for(float i =0; i<=35; i=i+0.5)
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
       
        sourceBits = sourceVector();

        transmittedSymbol=bit_maps_to_symbol_of_energy_E(sourceBits, energyOfSymbol[step],
                                                         one_million);

        realGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);
        imziGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);

        RayleighNoise = RayleighFaddingCoff(realGaussian, imziGaussian);
        
        gnoise = GnoiseVector(0.0, stdnoise, one_million);

        receiveSignal1 = ChannelOperation(transmittedSymbol, RayleighNoise, gnoise);

        realGaussian.clear();
        imziGaussian.clear();
        RayleighNoise.clear();
        gnoise.clear();

        realGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);
        imziGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);

        RayleighNoise = RayleighFaddingCoff(realGaussian, imziGaussian);
        
        gnoise = GnoiseVector(0.0, stdnoise, one_million);

        receiveSignal2 = ChannelOperation(transmittedSymbol, RayleighNoise, gnoise);

        realGaussian.clear();
        imziGaussian.clear();
        RayleighNoise.clear();
        gnoise.clear();

        receiveSignal= EqualGainCombining(receiveSignal1, receiveSignal2);

        decodedBits = decisionBlock(receiveSignal);
        counterror = errorCalculation(sourceBits, decodedBits);

        P_error = counterror/one_million;

        Prob_error.insert(Prob_error.end(), P_error);

        cout<<endl;

        cout<<"Energy of symbol     : "<<energyOfSymbol[step]<<endl;
        cout<<"Count errors         : "<<counterror<<endl;
        cout<<"Probability of error : "<<P_error<<endl;

        cout<<endl;
    }
    
    datafile(SNR_dB, Prob_error);
    // vector<double> qvalue = Errorfunction(SNR_dB);
    // ErrorValueInFile(SNR_dB, qvalue);

    return 0;
}