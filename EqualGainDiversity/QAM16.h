/////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- QAM MODEL AND "PROBABILITY OF Bit ERROR" vs "RATIO OF SIGNAL TO NOISE ENERGY (dB)"
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>
#include <fstream>

// #define one_million 1000000
#define pi 3.141592

using namespace std;



// This class contains all functions/methods for QAM
class QAM
{
    public:
        vector<double> Source();
        vector<double> TransSymbol_colVect1(vector<double> symbols, vector <double> energy, vector <double> angle);
        vector<double> TransSymbol_colVect2(vector<double> symbols, vector <double> energy, vector <double> angle);
        vector<double> GnoiseVector(double& mean, double& stddev, const int one_million);
        vector<double> ChannelModel(vector<double> sourceSymbols, vector<double> AWGnoise);
        vector <double> DecisionBlock_M4(vector<double> receive1, vector<double> receive2);
        vector <double> DecisionBlock_M8(vector<double> receive1, vector<double> receive2, double energy);
        vector <double> DecisionBlock_M16(vector<double> receive1, vector<double> receive2, double energy,vector<double> H1, vector<double> H2);
        
        int ErrorCount(vector <double> sourceSymbols, vector<double> decisionSymbols);

};




vector <double> EnergyVector(const double& M, double energy)
{
    vector <double> EnergyComponent;
    double magnitude;
    if(M==4)
    {
        magnitude = sqrt(2)*energy;
        for(int i =0; i<M; i++){
            EnergyComponent.insert(EnergyComponent.end(), magnitude);
        }
        return EnergyComponent;
    }
    if(M==8){
        magnitude = sqrt(2)*energy;
        for(int i =0; i<4; i++){
            EnergyComponent.insert(EnergyComponent.end(), magnitude);
        }
        magnitude = sqrt(10)*energy;
        for(int i =4; i<M; i++){
            EnergyComponent.insert(EnergyComponent.end(), magnitude);
        }
        return EnergyComponent;
    }
    if(M==16){
        magnitude = sqrt(2)*energy;
        for(int i =0; i<4; i++){
            EnergyComponent.insert(EnergyComponent.end(), magnitude);
        }
        magnitude = sqrt(10)*energy;
        for(int i =4; i<12; i++){
            EnergyComponent.insert(EnergyComponent.end(), magnitude);
        }
        magnitude = 3*sqrt(2)*energy;
        for(int i =12; i<M; i++){
            EnergyComponent.insert(EnergyComponent.end(), magnitude);
        }
        return EnergyComponent;
    }
    else{
        EnergyComponent.insert(EnergyComponent.end(), 0);
        return EnergyComponent;
    }
}



vector<double> AngleVector(const double& M )
{
    vector<double> angles;
    double angle;
    if(M == 4){
        angle = atan2(1,1);
        for(int i =0; i<M; i++){
            angles.insert(angles.end(), angle*180/pi);
        }
        return angles;
    }
    if(M == 8){
        angle = atan2(1,1);
        for(int i =0; i<4; i++){
            angles.insert(angles.end(), angle*180/pi);
        }
        angle = atan2(1,3);
        for(int i =4; i<M; i++){
            angles.insert(angles.end(), angle*180/pi);
        }
        return angles;
    }

    if(M == 16){
        angle = atan2(1,1);
        for(int i =0; i<4; i++){
            angles.insert(angles.end(), angle*180/pi);
        }
        angle = atan2(1,3);
        for(int i =4; i<8; i++){
            angles.insert(angles.end(), angle*180/pi);
        }
        angle = atan2(3,1);
        for(int i =8; i<12; i++){
            angles.insert(angles.end(), angle*180/pi);
        }
        angle = atan2(1,1);
        for(int i= 12; i<M; i++){
            angles.insert(angles.end(), angle*180/pi);
        }

        return angles;
    }
    else{
        angles.insert(angles.end(), 0);
        return angles;
    }
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



double radianConv(double degree){
    double rad = degree*pi/180;

    return rad;
}



vector<double>QAM::TransSymbol_colVect1(vector<double> symbols, vector <double> energy, vector <double> angle)
{
    vector<double> y1;
    double egy, ang;
    int index;
    for(int step = 0; step<symbols.size(); step++){
        index = symbols[step];
        egy = energy[index];

        // if((index)%4 == 0){
        //     ang = angle[index];
        // }
        // else if((index)%4 == 1){
        //     ang = 180- angle[index] ;
        // }
        // else if((index)%4 == 2){
        //     ang = angle[index] + 180;
        // }
        // else{
        //     ang = -1*angle[index] + 360;
        // }

        // ang = radianConv(ang);

        if(symbols[step] == 0){
            ang = 45;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));

        }else if(symbols[step] == 1){

            ang = 135;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));

        }else if(symbols[step] == 2){

            ang = 225;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));

            
        }else if(symbols[step] == 3){

            ang = 315;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 4){

            ang = 18.43;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 5){

            ang = 161.57;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 6){

            ang = 198.43;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 7){

            ang = 341.57;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 8){

            ang = 71.56;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 9){

            ang = 180-71.56;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 10){
            ang = 180+71.56;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 11){

            ang = 360-71.56;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 12){

            ang = 45;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 13){

            ang = 135;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 14){

            ang = 225;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
            
        }else if(symbols[step] == 15){
            
            ang = 315;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*cos(ang));
        }

        // y1.insert(y1.end(), egy*cos(ang));

    }

    return y1;

}



vector<double>QAM::TransSymbol_colVect2(vector<double> symbols, vector <double> energy, vector <double> angle)
{
    vector<double> y1;
    double egy, ang;
    int index;
    for(int step = 0; step<symbols.size(); step++){
        index = symbols[step];
        egy = energy[index];

        // if((index)%4 == 0){
        //     ang = angle[index];
        // }
        // else if((index)%4 == 1){
        //     ang = 180-angle[index];
        // }
        // else if((index)%4 == 2){
        //     ang = angle[index] + 180;
        // }
        // else{
        //     ang = -1*angle[index] +360;
        // }

        // ang = radianConv(ang);

        // y2.insert(y2.end(), egy*sin(ang));
         if(symbols[step] == 0){
            ang = 45;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));

        }else if(symbols[step] == 1){

            ang = 135;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));

        }else if(symbols[step] == 2){

            ang = 225;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));

            
        }else if(symbols[step] == 3){

            ang = 315;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 4){

            ang = 18.43;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 5){

            ang = 161.57;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 6){

            ang = 198.43;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 7){

            ang = 341.57;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 8){

            ang = 71.56;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 9){

            ang = 180-71.56;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 10){
            ang = 180+71.56;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 11){

            ang = 360-71.56;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 12){

            ang = 45;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 13){

            ang = 135;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 14){

            ang = 225;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
            
        }else if(symbols[step] == 15){
            
            ang = 315;
            ang = radianConv(ang);
            y1.insert(y1.end(), egy*sin(ang));
        }



    }

    return y1;

}




// Function for generating random noise based on gaussian distribution N(mean, variance).
// Input mean and standard deviation.
// Output is the vector that contain gaussian noise as an element.
vector<double>QAM::GnoiseVector(double& mean, double& stddev, const int one_million)
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
vector<double>QAM::ChannelModel(vector<double> sourceSymbols, vector<double> AWGnoise)
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
vector <double>QAM::DecisionBlock_M4(vector<double> receive1, vector<double> receive2)
{
    vector<double> decision;

    for(int j =0; j<receive1.size(); j++){

        if(receive1[j]>=0 && receive2[j]>=0){
            decision.insert(decision.end(), 0);
        }else if(receive1[j]< 0 && receive2[j]>=0){
            decision.insert(decision.end(), 1);
        }else if(receive1[j]< 0 && receive2[j]<0){
            decision.insert(decision.end(), 2);
        }else{
            decision.insert(decision.end(), 3);
        }
    }

    return decision;
}



// Function for take decision after receiving the vector
// Input is the 2D receiver vector
// Output is the decision {0, 1} 
vector <double>QAM::DecisionBlock_M8(vector<double> receive1, vector<double> receive2, double energy)
{
    vector<double> decision;

    for(int j =0; j<receive1.size(); j++){

        if(receive1[j]>=0 && receive1[j]<2*energy && receive2[j]>=0){
            decision.insert(decision.end(), 0);
        }else if(receive1[j]<0 && receive1[j]>-2*energy && receive2[j]>=0){
            decision.insert(decision.end(), 1);
        }else if(receive1[j]<0 && receive1[j]>-2*energy && receive2[j]<0){
            decision.insert(decision.end(), 2);
        }else if(receive1[j]>=0 && receive1[j]<2*energy && receive2[j]<0){
            decision.insert(decision.end(), 3);
        }else if(receive1[j]>=2*energy && receive2[j]>=0){
            decision.insert(decision.end(), 4);
        }else if(receive1[j]<= -2*energy && receive2[j]>=0){
            decision.insert(decision.end(), 5);
        }else if(receive1[j]<= -2*energy && receive2[j]<0){
            decision.insert(decision.end(), 6);
        }else{
            decision.insert(decision.end(), 7);
        }
    }

    return decision;
}



// Function for take decision after receiving the vector
// Input is the 2D receiver vector
// Output is the decision {0, 1} 
vector <double>QAM::DecisionBlock_M16(vector<double> receive1, vector<double> receive2, double energy, vector<double> H1, vector<double> H2)

{
    vector<double> decision;


    for(int j =0; j<receive1.size(); j++){

        if(receive1[j]>=0 && receive1[j]<2*energy && receive2[j]>=0 && receive2[j]<2*energy){
            decision.insert(decision.end(), 0);
        }else if(receive1[j]<0 && receive1[j]>-2*energy && receive2[j]>=0 && receive2[j]<2*energy){
            decision.insert(decision.end(), 1);
        }else if(receive1[j]<0 && receive1[j]>=-2*energy && receive2[j]<0 && receive2[j]>=-2*energy){
            decision.insert(decision.end(), 2);
        }else if(receive1[j]>=0 && receive1[j]<2*energy && receive2[j]<0 && receive2[j]>=-2*energy){
            decision.insert(decision.end(), 3);
        }else if(receive1[j]>=2*energy && receive2[j]>=0 && receive2[j]<2*energy){
            decision.insert(decision.end(), 4);
        }else if(receive1[j]>=0 &&receive1[j]<2*energy &&  receive2[j]>=2*energy){
            decision.insert(decision.end(), 8);
        }else if(receive1[j]<0 &&receive1[j]>=-2*energy &&  receive2[j]>=2*energy){
            decision.insert(decision.end(), 9);
        }else if(receive1[j]<-2*energy &&  receive2[j]>=0 &&receive2[j]<2*energy){
            decision.insert(decision.end(), 5);
        }else if(receive1[j]<-2*energy &&  receive2[j]<0 &&receive2[j]>=-2*energy){
            decision.insert(decision.end(), 6);
        }else if(receive1[j]>=-2*energy &&  receive1[j]<0 &&receive2[j]<-2*energy){
            decision.insert(decision.end(), 10);
        }else if(receive1[j]<2*energy &&  receive1[j]>=0 &&receive2[j]<-2*energy){
            decision.insert(decision.end(), 11);
        }else if(receive1[j]>=2*energy &&  receive2[j]<0 &&receive2[j]>=-2*energy){
            decision.insert(decision.end(), 7);
        }else if(receive1[j]>=2*energy &&  receive2[j]>=2*energy){
            decision.insert(decision.end(), 12);
        }else if(receive1[j]<-2*energy &&  receive2[j]>=2*energy){
            decision.insert(decision.end(), 13);
        }else if(receive1[j]<-2*energy &&  receive2[j]<-2*energy){
            decision.insert(decision.end(), 14);
        }else if(receive1[j]>=2*energy &&  receive2[j]<-2*energy){
            decision.insert(decision.end(), 15);
        }
    }

    return decision;
}




// Function for counting errors in the symbols
// Input are the source symbols and decided bits 
// Output is the error count  
int QAM::ErrorCount(vector <double> sourceSymbols, vector<double> decisionSymbols)
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

    outfile.open("QAM1_M8.dat");

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



// vector<double> Qfunction_M16(vector <double> SNR_dB)
// {
//     vector <double> Qvalue;
//     double po, normalValue;
//     for (int k =0; k<SNR_dB.size(); k++){
//         normalValue = pow(10, (SNR_dB[k]/10));
//         po = 3*Qfunc(sqrt(normalValue/5))-(9*pow(Qfunc(sqrt(normalValue)/5), 2)/4);
//         Qvalue.insert(Qvalue.end(), po);
//     }

//     return Qvalue;
// }




// Function for calculate the Q function values.
// Input is any positive real number.
// Output is the result of erfc function (equal form of Q function).
// double Qfunc(double x)
// {
//     double Qvalue = erfc(x/sqrt(2))/2;
//     return Qvalue;
// }



vector<double> Qfunction_M4(vector <double> SNR_dB)
{
    vector <double> Qvalue;
    double po, normalValue;
    for (int k =0; k<SNR_dB.size(); k++){
        normalValue = pow(10, (SNR_dB[k]/10));
        po = 2*Qfunc(sqrt(normalValue))-pow(Qfunc(sqrt(1*normalValue)),2);
        Qvalue.insert(Qvalue.end(), po);
    }

    return Qvalue;
}


vector<double> Qfunction_M8(vector <double> SNR_dB)
{
    vector <double> Qvalue;
    double po, normalValue;
    for (int k =0; k<SNR_dB.size(); k++){
        normalValue = pow(10, (SNR_dB[k]/10));
        po = ((5*(Qfunc(sqrt(normalValue/3))))-(3*pow(Qfunc(sqrt(normalValue/3)), 2)))/2;
        Qvalue.insert(Qvalue.end(), po);
    }

    return Qvalue;
}


vector<double> Qfunction_M16(vector <double> SNR_dB)
{
    vector <double> Qvalue;
    double po, normalValue, value;
    for (int k =0; k<SNR_dB.size(); k++){
        normalValue = pow(10, (SNR_dB[k]/10));
        value = Qfunc(sqrt(normalValue/5));
        po = ((12*value)-(9*pow(value, 2)))/4;
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

    outfile.open("QAM_Qvalue_M16.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<SNR.size(); i++){
        outfile<< SNR[i] << " "<<"\t"<<" "<< Qvalue[i]<< endl;
    }

    outfile.close();
}


