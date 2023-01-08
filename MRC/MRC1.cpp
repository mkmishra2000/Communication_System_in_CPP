/////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- Performance of maximum ratio combiner using BPSK in Rayleigh fadding channel
/////////////////////////////////////////////////////////////////////////////////////////////////

/*
Flow of the program
SourceBIts --> Modulation --> Multiply with Rayleigh (L different branch)
--> Add gaussian noise in each branch --> Multiply with Rayleigh in each branch
--> Combine all values -->Put into BPSK decoder --> Count errors.
*/

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


// Function for sum of all raws in a matrix 
// Input is a Matrix 
// Output is a vector result as sum of all raws 
vector<double> RawWiseSum(vector<vector<double>> ReceiveMat)
{
    int numberOfRaws = ReceiveMat.size();

    vector<double> SumofRaws= ReceiveMat[0];

    for(int t=0; t<numberOfRaws-1; t++){
        SumofRaws = VectorSum(SumofRaws, ReceiveMat[t+1]);
    }

    return SumofRaws;
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


double ProbabilityOferror(double SNR, double L)
{
    
    double inter = pow((SNR)+1, L);
    double pe = 1/inter;
    return pe;
}


vector<double> CalculatedError(vector<double> SNR_dB, double L)
{
    vector<double> ProbError;

    double po, normalValue, inter;
    for (int k =0; k<SNR_dB.size(); k++){
        normalValue = pow(10, (SNR_dB[k]/10));
        po = ProbabilityOferror(normalValue, L);
        ProbError.insert(ProbError.end(), po);
    }

    return ProbError;
}


int main(){

    // source defination
    vector<double> sourceBits;

    // Mapping of bits to symbols;
    vector<double> transmittedSymbol;

    // Noise definition
    vector<double> gnoise;

    //Rayleigh noise (Real guass and Img Guass)
    vector<double> realGaussian;
    vector<double> imziGaussian;
    vector<double> RayleighNoise;

    //Combiner output
    vector<double> AfterCombining;

    //Output of ML detector
    vector<double> decodedBits;

    double L =2;
    double sigmaSquare = 0.5;
    double stddevRayleigh = sqrt(sigmaSquare);
    double N_o =4;
    double p, stdnoise;
    stdnoise = sqrt(N_o);
    double counterror, P_error;


    vector<vector<double>> RaleighMat;
    vector<vector<double>> GnoiseMat;
    vector<vector<double>> faddingOperation;
    vector<vector<double>> ReceiveMat;
    vector<vector<double>> AfterCophasing;


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

    for(int step = 0; step<energyOfSymbol.size(); step++){

        sourceBits = sourceVector();    //Source bit streame

        transmittedSymbol = bit_maps_to_symbol_of_energy_E(sourceBits, 
                                                           energyOfSymbol[step], 
                                                           one_million);

        //Rayleigh distributed fadding coffcients
        for(int i =0; i<L; i++){
            realGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);
            imziGaussian = GnoiseVector(0.0, stddevRayleigh, one_million);

            RayleighNoise = RayleighFaddingCoff(realGaussian, imziGaussian);

            RaleighMat.insert(RaleighMat.end(), RayleighNoise);

            // realGaussian.clear();
            // imziGaussian.clear();
            // RayleighNoise.clear();
        }

        //Gaussian noise (white noise)
        for(int i=0; i<L; i++){
            gnoise = GnoiseVector(0.0, stdnoise, one_million);
            GnoiseMat.insert(GnoiseMat.end(), gnoise);
            gnoise.clear();
        }

        //Fadding operation
        for(int k =0; k<L; k++){
            faddingOperation.insert(faddingOperation.end(),
                                    RayleighOperation(transmittedSymbol, RaleighMat[k]));
        }
        
        //Additive gaussian noise
        for(int k =0; k<L; k++){
            ReceiveMat.insert(ReceiveMat.end(), 
                              GaussianNoiseAdd(faddingOperation[k], GnoiseMat[k]));
        }

        // After co-phaseing multiply with H_i 
        for(int k =0; k<L; k++){
            AfterCophasing.insert(AfterCophasing.end(),
                                    RayleighOperation(ReceiveMat[k], RaleighMat[k]));
        }

        //combining operation
        AfterCombining = RawWiseSum(AfterCophasing);

        //ML decoder
        decodedBits = decisionBlock(AfterCombining);

        //Count error
        counterror = errorCalculation(sourceBits, decodedBits);

        //Probability of error
        P_error = counterror/one_million;
        Prob_error.insert(Prob_error.end(), P_error);

        cout<<"Error value : "<<counterror<<endl;
        cout<<"Probability of error: "<<P_error<<endl;
        cout<<endl;

        RaleighMat.clear();
        GnoiseMat.clear();
        faddingOperation.clear();
        ReceiveMat.clear();
        AfterCophasing.clear();
    }

    char NameofFile1[30] = "MRCL2";         //Name of file

    char NameofFile2[30] = "MRCL2Error";    //Name of file


    datafile(SNR_dB, Prob_error, NameofFile1);

    vector<double> Error = CalculatedError(SNR_dB, L);

    datafile(SNR_dB, Error, NameofFile2);

    return 0;
}