//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- MPSK MODEL 
/////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <cmath>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>
#include <fstream>
#include <bits/stdc++.h>
#include <vector>

#define one_million 1000000
#define pi 3.14159

using namespace std;


// Function for generating binary bits at source side. Each bit is equiprobable.
// Input is nothing
// Output is a vector that contains the binary bits of length one_million*1.
vector<double> Source()
{
    vector<double> sourceBits;

    // Use current time as seed for random generator
    srand(time(0));
 
    for(int i = 0; i<one_million; i++){
        sourceBits.insert(sourceBits.end(), rand()%2);
    }

    return sourceBits;
}


// Function to generate symbols from the binary bit stream.
// Input is the binary bit vector and base for symbol conversion
// Output is the symbol vector containing values from (0, 2^(base-1))
vector <double> binaryToDecimalConversion(vector<double> sourceBits, const int& base)
{
    vector <double> convertedBits;
    int start =0;

    if((sourceBits.size()% base) == 0 ){

        int finalSize = sourceBits.size()/base;
        start = 0;
        int conversion;


        for(int i=0; i<finalSize; i++){
            conversion =0;
            for(int j =base-1; j>-1; j--){
                conversion = conversion + (sourceBits[start])*(pow(2, j));
                start++;
            }
            convertedBits.insert(convertedBits.end(), conversion);
        }

    }
    else{
        int addedBitsNO =  base - (sourceBits.size()%base);

        for(int q=0;q<addedBitsNO;q++){
            sourceBits.insert(sourceBits.end(), 0);
        }

        int finalSize = sourceBits.size()/base;
        start = 0;
        int conversion;


        for(int i=0; i<finalSize; i++){
            conversion =0;
            for(int j =0; j<base; j++){
                conversion = conversion + (sourceBits[start])*(pow(2, j));
                start++;
            }
            convertedBits.insert(convertedBits.end(), conversion);
        }

    }

    return convertedBits;
}




// Functions SignalVector1 and SignalVector2 for convert sourceSymbols to 2D-vector
// Input are the source Symbols and energy of symbols for mapping
// Output is a vector represent signal component
vector<double> SignalVectors1(vector<double> sourceSymbols, double eneryOfSymbols, const double& M)
{
    vector<double> y1;
    double theta;
    double angDiv = (2*pi)/M;

    for(int i=0; i<sourceSymbols.size(); i++){
        theta = sourceSymbols[i]*angDiv;

        y1.insert(y1.end(), sqrt(eneryOfSymbols)*cos(theta));
    }

    return y1;
}


vector<double> SignalVectors2(vector<double> sourceSymbols, double eneryOfSymbols, const double& M)
{
    vector<double> y2;
    double theta;
    double angDiv = (2*pi)/M;

    for(int i=0; i<sourceSymbols.size(); i++){
        theta = sourceSymbols[i]*angDiv;

        y2.insert(y2.end(), sqrt(eneryOfSymbols)*sin(theta));
    }


    return y2;
}



// Function for generating random noise based on gaussian distribution N(mean, variance).
// Input mean and standard deviation.
// Output is the vector that contain gaussian noise as an element.
vector<double> GnoiseVector(double mean, double stddev)
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
vector<double> channelModel(vector<double> sourceSymbols, vector<double> AWGnoise)
{
    vector<double> received;

    for(int i=0; i<sourceSymbols.size(); i++){
        received.insert(received.end(), sourceSymbols[i]+AWGnoise[i]);
    }

    return received;
}


double arcTan(double yvalue, double xvalue)
{
    double theta;

    if(yvalue==0 && xvalue ==0){
        theta =0;
    }else{
        theta = atan2(yvalue, xvalue)*(180/pi);
    }

    if(yvalue <0){
        theta = theta +360;
    }

    theta = theta*(pi/180);
    return theta;
}



// Function for take decision after receiving the vector
// Input is the 2D receiver vector
// Output is the decision {0, 1, 2, ..., M-1} 
vector <double> decisionBlock(vector<double> receive1, vector<double> receive2, double& M)
{
    vector<double> decision;

    double theta;
    double upperLimit;
    double lowerLimit;
    double increment;

    increment = 360/M;
    lowerLimit = 180/M;


    for( int index =0; index < receive1.size(); index++){
        theta = arcTan(receive2[index], receive1[index]);

        theta = theta*(180/pi);

        increment = 360/M;
        lowerLimit = 180/M;
        upperLimit = 180/M;


        double j = 1;
        while(j<=M){
            if(j==1 && (theta<=lowerLimit||theta-360>-lowerLimit)){
                decision.insert(decision.end(), j-1);
                j = j+1;
                upperLimit = lowerLimit + increment; 
                continue;
            }
            if(theta<= upperLimit && theta> lowerLimit){
                decision.insert(decision.end(), j-1);
            }
            j = j+1;
            lowerLimit = upperLimit;
            upperLimit = lowerLimit+increment;


        }

    }
    cout <<endl;

    return decision;
}



// Function for counting errors in the symbols
// Input are the source symbols and decided bits 
// Output is the error count  
int errorCount(vector <double> sourceSymbols, vector<double> decisionSymbols)
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
void datafile(vector<double> SNR, vector<double> Prob_error)
{
    ofstream outfile;

    outfile.open("MPSK(m=8)Data.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<SNR.size(); i++){
        outfile<< SNR[i] << " "<<"\t"<<" "<< Prob_error[i]<< endl;
    }

    outfile.close();
}




// Function to compute Q function value.
// Input is the x.
// Output is Q function output.
double Qfunc(double x)
{
    double Qvalue = erfc(x/sqrt(2))/2;
    return Qvalue;
}



vector<double> Qfunction(vector <double> SNR_dB, const double& M)
{
    vector <double> Qvalue;
    double po, normalValue, extra;
    for (int k =0; k<SNR_dB.size(); k++){

        normalValue = pow(10, (SNR_dB[k]/10));
        extra = sin(pi/M);
        po = (2*Qfunc(sqrt(2*normalValue)*extra)); // Defination of Pe from calculations.
        Qvalue.insert(Qvalue.end(), po);
    }

    return Qvalue;
}


void qvalueInFile(vector <double> SNR, vector <double> Qvalue)
{
    ofstream outfile;

    outfile.open("MPSK(m=8)_Qvalue.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<SNR.size(); i++){
        outfile<< SNR[i] << " "<<"\t"<<" "<< Qvalue[i]<< endl;
    }

    outfile.close();
}



int main()
{
    // source defination
    vector<double> sourceBits;

    // Symbols formation
    vector<double> Symbols;

    // Base for conversion 
    int base =3; //3, 4, 5, 6 ...
    double M = pow(2, base);

    // Transmiited vectors
    vector<double> trans1;
    vector<double> trans2;

    int SizeOfTransmission = trans1.size();

    // Noise vector  
    vector<double> gnoise1;
    vector<double> gnoise2;

    // Receive bits 
    vector<double> receivedBits1;
    vector<double> receivedBits2;


    // Decision block
    vector<double> decodedBits;

    vector<double> EnergyVector;

    vector<double> p_error;

    double errors=0;
    double pe;
    double N_o = 8;
    double stddev = sqrt(N_o/2);


   // SNR in dB 
    vector<double> SNR_dB;
    for(float i =0; i<=20; i=i+0.125)
    {
        SNR_dB.insert(SNR_dB.end(), i);
    }


    // copy(begin(SNR), end(SNR), std::ostream_iterator<double>(std::cout, "   "));

    double normalValue;

    for(int i =0; i<SNR_dB.size(); i++){
        normalValue = pow(10, (SNR_dB[i]/10));
        EnergyVector.insert(EnergyVector.end(), normalValue*N_o);
    }

    // copy(begin(EnergyVector), end(EnergyVector), std::ostream_iterator<double>(std::cout, "   "));


    for(int step =0; step <SNR_dB.size(); step++){

        sourceBits = Source();

        Symbols = binaryToDecimalConversion(sourceBits, base);
        // cout<< "Done symbols"<<endl;

        trans1 = SignalVectors1(Symbols, EnergyVector[step], M);
        trans2 = SignalVectors2(Symbols, EnergyVector[step], M);

        gnoise1 = GnoiseVector(0.0, stddev);
        gnoise2 = GnoiseVector(0.0, stddev);

        receivedBits1 = channelModel(trans1, gnoise1);
        receivedBits2 = channelModel(trans2, gnoise2);


        decodedBits = decisionBlock(receivedBits1, receivedBits2, M);


        errors = errorCount(Symbols, decodedBits);
        double sizeOfsymbols = Symbols.size();
        pe = errors/ sizeOfsymbols;

        p_error.insert(p_error.end(), pe);

        cout<< errors<< "\n";
        cout<< pe << "\n";

        cout<< "\n";

    }


    datafile(SNR_dB,p_error);

    vector<double> qvalue = Qfunction(SNR_dB, M);
    qvalueInFile(SNR_dB, qvalue);

    return 0;


}