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


// Function for finding index that contain maximum value in every column
// Input is a matrix 
// Output is the vector with each element as index 
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



// function for factorial of n
// Input is a integer number greater than 0
// Output is the factorial result 
double factorial(double Num)
{
    if(Num==1){
        return 1;
    }
    else if(Num ==0){
        return 1;
    }
    else if(Num<0){
        cout<< "Worng Output"<<endl;
        return 0;
    }
    else{
        return (Num*factorial(Num-1));
    }
}


// Function for the combination l_C_k
// Inputs are the two integer numbers l and k 
// Output is the answer of the l choose K.
double l_Choose_K(double l, double k)
{
    double ans;
    if(l<k){
        cout<<"L should not be less than k!!!"<<endl;
        return 0;
    }
    else{
        ans = factorial(l)/(factorial(k)*factorial(l-k));
        return ans;
    }
}

double errorCalculation(double l, double snr)
{
    vector<double> calerror;
    double fac, val, error;

    for(int i =0; i<=l; i++){
        val = 1/pow((1+(2*i/snr)),0.5);
        fac = l_Choose_K(l,i);
        val = pow(-1, i)*(val/2);
        calerror.insert(calerror.end(), fac*val);
    }
    error =0;
    for(int i =0; i<calerror.size(); i++){
        error = error + calerror[i];
    }

    calerror.clear();

    return error;
}



vector<double> calculatedError(vector <double> SNR_dB, double l)
{

    vector <double> Qvalue;
    
    double po, normalValue, inter;
    for (int k =0; k<SNR_dB.size(); k++){
        normalValue = pow(10, (SNR_dB[k]/10));
        po = errorCalculation(l, normalValue);
        Qvalue.insert(Qvalue.end(), po);
    }

    return Qvalue;

}

// Function to store the data in the file (.dat)
// Input is the SNR per bit in dB and calculated Qfunction values
// Output is the nothing but in processing it is creating a file and writing data into it.
void ErrorValueInFile(vector <double> SNR, vector <double> Qvalue)
{
    ofstream outfile;

    outfile.open("SGcalL3.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<SNR.size(); i++){
        outfile<< SNR[i] << " "<<"\t"<<" "<< Qvalue[i]<< endl;
    }

    outfile.close();
}


// Function to store the data in the file (.dat)
// Input is the SNR per bit in dB and calculated probability of error
// Output is the nothing but in processing it is creating a file and writing data into it.
void datafile(vector<double> xindB, vector<double> Prob_error)
{
    ofstream outfile;

    outfile.open("SelectiveGainL3.dat");

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

    vector<double> realGaussian;
    vector<double> imziGaussian;
    vector<double> RayleighNoise;

    vector<double> MaxiSNRdetection;

    vector<double> combinerOutput;
    vector<double> decodedBits;

    double sigmaSquare = 0.5;
    double stddevRayleigh = sqrt(sigmaSquare);
    double N_o =4;
    double p, stdnoise;
    stdnoise = sqrt(N_o);
    double counterror, P_error;

    double L =3;

    vector<vector<double>> MultipleSignals;
    vector<vector<double>> RaleighMat;
    vector<vector<double>> SNRMat;
    vector<vector<double>> GnoiseMat;
    vector<vector<double>> NewMat;
    vector<vector<double>> ReceiveMat;

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


    for(int step =0; step<energyOfSymbol.size(); step++){
        
        sourceBits = sourceVector();

        transmittedSymbol = bit_maps_to_symbol_of_energy_E(sourceBits, energyOfSymbol[step], one_million);

        for(int i =0; i<L; i++){
            MultipleSignals.insert(MultipleSignals.end(), transmittedSymbol);
        }

        for(int i =0; i<L; i++){
            realGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);
            imziGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);

            RayleighNoise = RayleighFaddingCoff(realGaussian, imziGaussian);

            RaleighMat.insert(RaleighMat.end(), RayleighNoise);

            // realGaussian.clear();
            // imziGaussian.clear();
            // RayleighNoise.clear();
        }

        for(int i=0; i<L; i++){
            gnoise = GnoiseVector(0.0, stdnoise, one_million);
            GnoiseMat.insert(GnoiseMat.end(), gnoise);
            gnoise.clear();
        }

        NewMat = ElementWiseMultiplication(MultipleSignals, RaleighMat);

        ReceiveMat = ElementWiseAddition(NewMat, GnoiseMat);

        SNRMat = ElementWiseMultiplication(RaleighMat, RaleighMat);

        MaxiSNRdetection = MaximumSNRIndex(SNRMat);

        for(int j =0; j<MaxiSNRdetection.size(); j++){
            combinerOutput.insert(combinerOutput.end(), ReceiveMat[MaxiSNRdetection[j]][j]);
        }

        decodedBits = decisionBlock(combinerOutput);

        counterror = errorCalculation(sourceBits, decodedBits);

        P_error = counterror/one_million;
        Prob_error.insert(Prob_error.end(), P_error);

        cout<<"Error count: "<< counterror<<endl;
        cout<<"P_e:         "<< P_error<<endl;
        cout<<endl;

        MultipleSignals.clear();
        NewMat.clear();
        RaleighMat.clear();
        GnoiseMat.clear();
        ReceiveMat.clear();
        SNRMat.clear();
        combinerOutput.clear();
    }

    datafile(SNR_dB, Prob_error);
    vector<double> qvalue = calculatedError(SNR_dB, L);
    ErrorValueInFile(SNR_dB, qvalue);

    return 0;
}