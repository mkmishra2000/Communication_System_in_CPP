#include <iostream>
#include <cmath>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>
#include <fstream>
#include <bits/stdc++.h>
#include <vector>

// #define one_million 100000
#define pi 3.14159

using namespace std;


// Function to generate symbols from the binary bit stream.
// Input is the binary bit vector and base for symbol conversion
// Output is the symbol vector containing values from (0, 2^(base-1))
vector <double> binaryToDecimalConversion(vector<double> sourceBits, int base)
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
vector<double> SignalVectors1(vector<double> sourceSymbols, double eneryOfSymbols, int M)
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


vector<double> SignalVectors2(vector<double> sourceSymbols, double eneryOfSymbols ,int M)
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
vector<double> GnoiseVector(double mean, double stddev, int sizeOfTransmission, const int one_million)
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
vector <double> decisionBlock(vector<double> receive1, vector<double> receive2, double Energy, vector<double> H1, vector<double> H2, int M)
{
    vector<double> decision;

    double theta;
    double Phi;

    // double values;
    // vector <double> angles;

    // for(int k =0; k<M; k++)
    // {
    //     values = (2*k+1)*(180/M);
    //     angles.insert(angles.end(), values);
    // }

    for( int index =0; index < receive1.size(); index++){
        theta = arcTan(receive2[index], receive1[index]);

        theta = theta*(180/pi);
        
        Phi = arcTan(H2[index]*sqrt(1*Energy), H1[index]*sqrt(1*Energy));

        Phi = Phi*(180/pi);
        Phi = Phi/2;

        // cout << theta << " degree  ";

        // for(int j =1; j<M; j++ ){
        //     if(j==1 && (theta<=(2*j-1)*180/M||theta-360>-(2*j-1)*180/M)){
        //         decision.insert(decision.end(), j-1);
        //     }else if(theta<=(2*j+1)*180/M && theta>(2*j-1)*180/M){
        //         decision.insert(decision.end(), j);
        //     }
        //     // else if(theta<=(2*j-1)*180/M && theta>(2*j-3)*180/M){
        //     //     decision.insert(decision.end(), j-1);
        //     // }
        //     // else if(theta<=(2*j-1)*180/M && theta>(2*j-3)*180/M){
        //     //      decision.insert(decision.end(), j-1);
        //     // }

        // }

        if(theta<=Phi|| theta-360>-Phi){
            decision.insert(decision.end(), 0);

        }else if(theta<=45+Phi && theta>Phi){
            decision.insert(decision.end(), 1);

        }else if(theta<=(90+Phi)&& theta>(45+Phi)){
            decision.insert(decision.end(), 2);

        }else if(theta<=(135+Phi) && theta>(90+Phi)){
            decision.insert(decision.end(), 3);

        }else if(theta<=(180+Phi)&& theta>(135+Phi)){
            decision.insert(decision.end(), 4);

        }else if(theta<=(225+Phi)&& theta>(180+Phi)){
            decision.insert(decision.end(), 5);

        }else if(theta<=(270+Phi)&&theta>(225+Phi)){
            decision.insert(decision.end(), 6);

        }else if(theta<=(315+Phi)&& theta>(270+Phi)){
            decision.insert(decision.end(), 7);

        }

    }
    // cout<<endl;

    return decision;
}
