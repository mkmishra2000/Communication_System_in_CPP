/////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- Performance of maximum ratio combiner using BPSK in Rayleigh fadding channel
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>
#include <fstream>
#include "BFSK.h"
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
    
    double inter = pow(SNR+1, L);
    double pe = (L-1)/inter;
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

    vector<double> Symbols;

    // Mapping of bits to symbols first vector;
    vector<double> transmittedSymbol1;

    // Mapping of bits to symbols first vector;
    vector<double> transmittedSymbol2;

    // Noise definition
    vector<double> gnoise;
    

    //Rayleigh noise (Real guass and Img Guass)
    vector<double> realGaussian;
    vector<double> imziGaussian;
    vector<double> RayleighNoise;

    //Combiner output
    vector<double> AfterCombining1;

    //Combiner output
    vector<double> AfterCombining2;

    //Output of ML detector
    vector<double> decodedBits;

    double L = 8;
    int Base =2;
    double sigmaSquare = 0.5;
    double stddevRayleigh = sqrt(sigmaSquare);
    double N_o =4;
    double p, stdnoise;
    stdnoise = sqrt(N_o);
    double mean =0;
    double counterror, P_error;


    vector<vector<double>> RaleighMat1;
    vector<vector<double>> GnoiseMat1;
    vector<vector<double>> faddingOperation1;
    vector<vector<double>> ReceiveMat1;
    vector<vector<double>> AfterCophasing1;

    vector<vector<double>> RaleighMat2;
    vector<vector<double>> GnoiseMat2;
    vector<vector<double>> faddingOperation2;
    vector<vector<double>> ReceiveMat2;
    vector<vector<double>> AfterCophasing2;



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

        // Symbols = binaryToDecimalConversion(sourceBits, Base);

        transmittedSymbol1 = TransmittedSymbol_columnVector1(sourceBits, energyOfSymbol[step], one_million);
        transmittedSymbol2 = TransmittedSymbol_columnVector2(sourceBits, energyOfSymbol[step], one_million);

        //Rayleigh distributed fadding coffcients
        for(int i =0; i<L; i++){
            realGaussian = GnoiseVector(mean, stddevRayleigh, one_million);
            imziGaussian = GnoiseVector(mean, stddevRayleigh, one_million);

            RayleighNoise = RayleighFaddingCoff(realGaussian, imziGaussian);

            RaleighMat1.insert(RaleighMat1.end(), RayleighNoise);

            // realGaussian.clear();
            // imziGaussian.clear();
            // RayleighNoise.clear();
        }

        for(int i =0; i<L; i++){
            realGaussian = GnoiseVector(mean, stddevRayleigh, one_million);
            imziGaussian = GnoiseVector(mean, stddevRayleigh, one_million);

            RayleighNoise = RayleighFaddingCoff(realGaussian, imziGaussian);

            RaleighMat2.insert(RaleighMat2.end(), RayleighNoise);

            // realGaussian.clear();
            // imziGaussian.clear();
            // RayleighNoise.clear();
        }


        //Gaussian noise (white noise)
        for(int i=0; i<L; i++){
            gnoise = GnoiseVector(mean, stdnoise, one_million);
            GnoiseMat1.insert(GnoiseMat1.end(), gnoise);
            gnoise.clear();
        }

        //Gaussian noise (white noise)
        for(int i=0; i<L; i++){
            gnoise = GnoiseVector(mean, stdnoise, one_million);
            GnoiseMat2.insert(GnoiseMat2.end(), gnoise);
            gnoise.clear();
        }

        //Fadding operation
        for(int k =0; k<L; k++){
            faddingOperation1.insert(faddingOperation1.end(),
                                    RayleighOperation(transmittedSymbol1, RaleighMat1[k]));
        }
        
        //Fadding operation
        for(int k =0; k<L; k++){
            faddingOperation2.insert(faddingOperation2.end(),
                                    RayleighOperation(transmittedSymbol2, RaleighMat2[k]));
        }
        
        //Additive gaussian noise
        for(int k =0; k<L; k++){
            ReceiveMat1.insert(ReceiveMat1.end(), 
                              GaussianNoiseAdd(faddingOperation1[k], GnoiseMat1[k]));
        }

        //Additive gaussian noise
        for(int k =0; k<L; k++){
            ReceiveMat2.insert(ReceiveMat2.end(), 
                              GaussianNoiseAdd(faddingOperation2[k], GnoiseMat2[k]));
        }


        //Fadding operation
        for(int k =0; k<L; k++){
            AfterCophasing1.insert(AfterCophasing1.end(),
                                    RayleighOperation(ReceiveMat1[k], RaleighMat1[k]));
        }
        
        //Fadding operation
        for(int k =0; k<L; k++){
            AfterCophasing2.insert(AfterCophasing2.end(),
                                    RayleighOperation(ReceiveMat2[k], RaleighMat2[k]));
        }

        //combining operation
        AfterCombining1 = RawWiseSum(AfterCophasing1);

        //combining operation
        AfterCombining2 = RawWiseSum(AfterCophasing2);

        // RaleighMat1SQr = ElementWiseMultiplication(RaleighMat1, RaleighMat1);
        // RaleighMat2SQr = ElementWiseMultiplication(RaleighMat2, RaleighMat2);

        // SumofHi1 = RawWiseSum(RaleighMat1SQr);
        // SumofHi2 = RawWiseSum(RaleighMat2SQr);

        //ML decoder
        decodedBits = DecisionBlock(AfterCombining1, AfterCombining2);

        //Count error
        counterror = errorCalculation(sourceBits, decodedBits);

        //Probability of error
        P_error = counterror/sourceBits.size();
        Prob_error.insert(Prob_error.end(), P_error);

        cout<<"Error value : "<<counterror<<endl;
        cout<<"Probability of error: "<<P_error<<endl;
        cout<<endl;

        RaleighMat1.clear();
        GnoiseMat1.clear();
        faddingOperation1.clear();
        ReceiveMat1.clear();
        AfterCophasing1.clear();

        RaleighMat2.clear();
        GnoiseMat2.clear();
        faddingOperation2.clear();
        ReceiveMat2.clear();
        AfterCophasing2.clear();
    }

    char NameofFile1[30] = "MRCBFSKL8";
    // char NameofFile2[30] = "EGCL7Error";


    datafile(SNR_dB, Prob_error, NameofFile1);

    return 0;
}