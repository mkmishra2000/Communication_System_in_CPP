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

#define one_million 1000000

using namespace std;

// This class contains all functions/methods for BFSK
class BFSK
{
    public:
        vector<double> Source();
        vector<double> TransmittedSymbol_columnVector1(vector<double> sourceBits, double& energyOfbit);
        vector<double> TransmittedSymbol_columnVector2(vector<double> sourceBits, double& energyOfbit);
        vector<double> GnoiseVector(double& mean, double& stddev);
        vector<double> ChannelModel(vector<double> sourceSymbols, vector<double> AWGnoise);
        vector <double> DecisionBlock(vector<double> receive1, vector<double> receive2);
        int ErrorCount(vector <double> sourceSymbols, vector<double> decisionSymbols);

};



// Function for generating binary bits at source side. Each bit is equiprobable.
// Input is nothing
// Output is a vector that contains the binary bits of length one_million*1.
vector<double>BFSK::Source()
{
    vector<double> sourceBits;

    // Use current time as seed for random generator
    srand(time(0));
 
    for(int i = 0; i<one_million; i++){
        sourceBits.insert(sourceBits.end(), rand()%2);
    }

    return sourceBits;
}



// Function of generate vector form of signal, this function generate the upper value of the vector.
// Input is the source bit stream and energy of bit.
// Output is the vector containing all upper values based on the source bits. 
vector<double>BFSK::TransmittedSymbol_columnVector1(vector<double> sourceBits, double& energyOfbit)
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
vector<double>BFSK::TransmittedSymbol_columnVector2(vector<double> sourceBits, double& energyOfbit)
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
vector<double>BFSK::GnoiseVector(double& mean, double& stddev)
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




// Function for model the channel
// Input sourcesymbols and AWGN
// Output is a vector that represent the addition of two vectors
vector<double>BFSK::ChannelModel(vector<double> sourceSymbols, vector<double> AWGnoise)
{
    vector<double> received;

    for(int i=0; i<sourceSymbols.size(); i++){
        received.insert(received.end(), sourceSymbols[i]+AWGnoise[i]);
    }

    return received;
}



// Function for take decision after receiving the vector
// Input is the 2D receiver vector
// Output is the decision {0, 1} 
vector <double>BFSK::DecisionBlock(vector<double> receive1, vector<double> receive2)
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



// Function for counting errors in the symbols
// Input are the source symbols and decided bits 
// Output is the error count  
int BFSK::ErrorCount(vector <double> sourceSymbols, vector<double> decisionSymbols)
{
    int count =0;
    
    for(int i=0; i<sourceSymbols.size(); i++){
        if(sourceSymbols[i] != decisionSymbols[i]){
            count++;
        }
    }

    return count;
}



// Function to store the data in the file (.dat)
// Input is the SNR per bit in dB and calculated probability of error
// Output is the nothing but in processing it is creating a file and writing data into it.
void datafile(vector<double> xindB, vector<double> Prob_error)
{
    ofstream outfile;

    outfile.open("BFSK1.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<xindB.size(); i++){
        outfile<< xindB[i] << " "<<"\t"<<" "<< Prob_error[i]<< endl;
    }

    outfile.close();
}


// Function for calculate the Q function values.
// Input is any positive real number.
// Output is the result of erfc function (equal form of Q function).
double Qfunc(double x)
{
    double Qvalue = erfc(x/sqrt(2))/2;
    return Qvalue;
}



vector<double> Qfunction(vector <double> SNR_dB)
{
    vector <double> Qvalue;
    double po, normalValue;
    for (int k =0; k<SNR_dB.size(); k++){
        normalValue = pow(10, (SNR_dB[k]/10));
        po = Qfunc(sqrt(1*normalValue));
        Qvalue.insert(Qvalue.end(), po);
    }

    return Qvalue;
}



// Function to store the data in the file (.dat)
// Input is the SNR per bit in dB and calculated Qfunction values
// Output is the nothing but in processing it is creating a file and writing data into it.
void qvalueInFile(vector <double> SNR, vector <double> Qvalue)
{
    ofstream outfile;

    outfile.open("BFSK_Qvalue1.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<SNR.size(); i++){
        outfile<< SNR[i] << " "<<"\t"<<" "<< Qvalue[i]<< endl;
    }

    outfile.close();
}




int main(){


    // Object of class 
    BFSK bfsk1;

    // source defination
    vector<double> sourceBits;

    // Transmiited vectors
    vector<double> trans1;
    vector<double> trans2;

    // Noise vector  
    vector<double> gnoise1;
    vector<double> gnoise2;

    // Receive bits 
    vector<double> receivedBits1;
    vector<double> receivedBits2;

    // Decision block
    vector<double> decodedBits;

    double N_o = 8;
    double stddev = sqrt(N_o/2);
    double mean = 0.0;
    double errors=0;
    double pe;


    // SNR in dB 
    vector<double> SNR_dB;
    for(float i =0; i<=14; i=i+0.125)
    {
        SNR_dB.insert(SNR_dB.end(), i);
    }

    vector<double> energyOfBits;
    vector<double> Prob_error;
    double normalValue;

    for(int i =0; i<SNR_dB.size(); i++){

        normalValue = pow(10, (SNR_dB[i]/10));
        energyOfBits.insert(energyOfBits.end(), N_o*normalValue);
    }

    for(int step =0; step<energyOfBits.size(); step++){
        sourceBits = bfsk1.Source();

        trans1 = bfsk1.TransmittedSymbol_columnVector1(sourceBits, energyOfBits[step]);
        trans2 = bfsk1.TransmittedSymbol_columnVector2(sourceBits, energyOfBits[step]);

        gnoise1 = bfsk1.GnoiseVector(mean, stddev);
        gnoise2 = bfsk1.GnoiseVector(mean, stddev);

        receivedBits1 = bfsk1.ChannelModel(trans1, gnoise1);
        receivedBits2 = bfsk1.ChannelModel(trans2, gnoise2);

        decodedBits = bfsk1.DecisionBlock(receivedBits1, receivedBits2);

        errors = bfsk1.ErrorCount(sourceBits, decodedBits);

        pe = errors/one_million;

        Prob_error.insert(Prob_error.end(), pe);

        cout<< " Error count : "<<errors<< endl;
        cout<< " Pe : "<<pe << endl;

        cout<<endl;
    }
    
    datafile(SNR_dB, Prob_error);
    vector<double> qvalue = Qfunction(SNR_dB);
    qvalueInFile(SNR_dB, qvalue);

    return 0;
}