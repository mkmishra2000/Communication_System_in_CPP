//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- MPSK MODEL AND "PROBABILITY OF SYMBOL ERROR" vs "RATIO OF SIGNAL TO NOISE ENERGY (dB)"
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


vector<double> SignalVectors2(vector<double> sourceSymbols, double eneryOfSymbols, int M)
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
vector<double> GnoiseVector(double mean, double stddev, int sizeOfTransmission)
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


vector <double> DistanceVector(double point1, double point2, double energy, vector<double> angle)
{
    vector<double> distance;
    double dis;
    double r;
    double theta = arcTan(point2, point1);

    // vector<double> constellationpoint1;
    // vector <double> constellationpoint2;

    // for(int i =0; i<angle.size(); i++){
    //     constellationpoint1.insert(constellationpoint1.end(), sqrt(energy)*cos(angle[i]));
    //     constellationpoint2.insert(constellationpoint2.end(), sqrt(energy)*sin(angle[i]));
    // }

    for(int j =0; j<angle.size(); j++){

        r = sqrt(pow(point1, 2)+pow(point2, 2));
        dis =energy + pow(r, 2) - sqrt(energy)*r*cos(theta-angle[j]);

        distance.insert(distance.end(), dis);
    }

    return distance;
}

vector<double> angleDifference(double x, double y, vector<double> angle)
{
    double theta = arcTan(y, x);
    double check;

    vector <double> angleDifference;

    for(int j =0; j<angle.size(); j++){
        check = angle[j] - theta;
        if(check<0){
            angleDifference.insert(angleDifference.end(), 4*pi);
        }else{
            angleDifference.insert(angleDifference.end(), check);
        }
    }


    return angleDifference;

}

int minElement_indexReturn (vector <double> distance)
{ 
    int minElementIndex = min_element(distance.begin(),distance.end()) - distance.begin();

    return minElementIndex;
}

// Function for take decision after receiving the vector
// Input is the 2D receiver vector
// Output is the decision {0, 1, 2, ..., M-1} 
vector <double> decisionBlock(vector<double> receive1, vector<double> receive2, double Energy, int M)
{
    vector<double> decision;

    double theta;

    double values;
    vector <double> angles;
    // int min_index_distance;
    // int min_index_angle;
    for(int k =0; k<M; k++)
    {
        values = (2*k+1)*(180/M);
        angles.insert(angles.end(), values);
    }

    // vector <double>distance;
    // vector<double> angledifference;
    

    // for(int i =0; i<receive1.size(); i++){
    //     distance = DistanceVector(receive1[i], receive2[i], Energy, angles);
    //     angledifference = angleDifference(receive1[i], receive2[i], angles);

    //     cout<< "distance"<< endl;
    //     copy(begin(distance), end(distance), std::ostream_iterator<double>(std::cout, "   "));
    //     cout<<endl;

    //     cout<< "Angle"<< endl;
    //     copy(begin(angledifference), end(angledifference), std::ostream_iterator<double>(std::cout, "   "));
    //     cout<<endl;


    //     min_index_distance =minElement_indexReturn(distance);
    //     min_index_angle = minElement_indexReturn(angledifference);

    //     if(min_index_angle == min_index_distance){
    //         decision.insert(decision.end(), min_index_distance);
    //     }else{

    //         decision.insert(decision.end(), min_index_distance);
    //     }

    // }
    // cout<< endl;

    for( int index =0; index < receive1.size(); index++){
        theta = arcTan(receive2[index], receive1[index]);
        // cout<< theta*(180/pi) << "degree  ";

        theta = theta*(180/pi);

        for(int j =1; j<=M; j++ ){
            if(j==1 && (theta<=(2*j-1)*180/M||theta-360>-(2*j-1)*180/M)){
                decision.insert(decision.end(), j-1);
            }else if(theta<=(2*j-1)*180/M && theta>(2*j-3)*180/M){
                decision.insert(decision.end(), j-1);
            }
            else if(theta<=(2*j-1)*180/M && theta>(2*j-3)*180/M){
                decision.insert(decision.end(), j-1);
            }
            else if(theta<=(2*j-1)*180/M && theta>(2*j-3)*180/M){
                 decision.insert(decision.end(), j-1);
            }

        }

        // if(theta<=(2*1-1)*180/M||theta-360>-(2*1-1)*180/M){
        //     decision.insert(decision.end(), 0);

        // }
        // else if(theta<=(2*2-1)*180/M && theta>(2*1-1)*180/M){
        //     decision.insert(decision.end(), 1);
        // }
        // else if(theta<=(2*3-1)*180/M && theta>(2*2-1)*180/M){
        //     decision.insert(decision.end(), 2);
        // }
        //  else if(theta<=(2*4-1)*180/M && theta>(2*3-1)*180/M){
        //      decision.insert(decision.end(), 3);
        // }

            // for(int j =0; j<angles.size()-1; j++){

            //     if(j ==0 ){
            //         if(theta<angles[0]||theta-360 > - angles[0]){
            //             decision.insert(decision.end(), 0);
            //         }
            //     }else{
            //         if(theta<=angles[j+1] && theta > angles[j]){
            //         decision.insert(decision.end(), j);
            //         }
            //     }

            // }
                



        // for(int j =0; j<M-1; j++){
        //     if(angles[j]<theta && angles[j+1]>= theta){
        //         decision.insert(decision.end(), j+1);
        //     }
        // }

        // if(theta <= angles[0] || theta>angles[M-1]){
        //     decision.insert(decision.end(), 0);
        // } 
    }
    cout<<endl;

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

    outfile.open("Mtest.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<SNR.size(); i++){
        outfile<< SNR[i] << " "<<"\t"<<" "<< Prob_error[i]<< endl;
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
    int base =2;
    int M = pow(2, base);

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

    vector <double> SNR;

    vector<double> EnergyVector;

    vector<double> p_error;

    double errors=0;
    double pe;
    double N_o = 2;
    double stddev = sqrt(N_o/2);


    for (double i =0.5; i<12; i=i+0.07){
        SNR.insert(SNR.end(), i);
       
    }
    // copy(begin(SNR), end(SNR), std::ostream_iterator<double>(std::cout, "   "));

    for(int i =0; i<SNR.size(); i++){
        EnergyVector.insert(EnergyVector.end(), SNR[i]*N_o);
    }

    // copy(begin(EnergyVector), end(EnergyVector), std::ostream_iterator<double>(std::cout, "   "));


    for(int step =0; step <SNR.size(); step++){

        sourceBits = Source();

        Symbols = binaryToDecimalConversion(sourceBits, base);
        // cout<< "Done symbols"<<endl;

        trans1 = SignalVectors1(Symbols, EnergyVector[step], M);
        trans2 = SignalVectors2(Symbols, EnergyVector[step], M);

        gnoise1 = GnoiseVector(0.0, stddev, SizeOfTransmission);
        gnoise2 = GnoiseVector(0.0, stddev, SizeOfTransmission);

        receivedBits1 = channelModel(trans1, gnoise1);
        receivedBits2 = channelModel(trans2, gnoise2);

        // cout<< "Done received bits"<<endl;

        // cout<< "Transmitted symbol 1"<<endl;
        // copy(begin(trans1), end(trans1), std::ostream_iterator<double>(std::cout, "   "));


        // cout<<endl;

        // cout<< "Transmitted symbol 2"<<endl;
        // copy(begin(trans2), end(trans2), std::ostream_iterator<double>(std::cout, "   "));
        // cout<<endl;


        // cout<< "Gaussian noise 1"<<endl;
        // copy(begin(gnoise1), end(gnoise1), std::ostream_iterator<double>(std::cout, "   "));
        // cout<<endl;


        // cout<< "Gaussian noise 2"<<endl;
        // copy(begin(gnoise2), end(gnoise2), std::ostream_iterator<double>(std::cout, "   "));
        // cout<<endl;


        // cout<< "Received symbol1"<<endl;
        // copy(begin(receivedBits1), end(receivedBits1), std::ostream_iterator<double>(std::cout, "   "));

        // cout<<endl;

        // cout<< "Received symbol 2"<<endl;
        // copy(begin(receivedBits2), end(receivedBits2), std::ostream_iterator<double>(std::cout, "   "));
        // cout<<endl;
        decodedBits = decisionBlock(receivedBits1, receivedBits2, EnergyVector[step], M);

        errors = errorCount(Symbols, decodedBits);
        double sizeOfsymbols = Symbols.size();
        pe = errors/ sizeOfsymbols;

        p_error.insert(p_error.end(), pe);

        cout<< errors<< "\n";
        cout<< pe << "\n";
        // cout<< EnergyVector[step]<<endl;

        // cout<<"Symbols\n";
        // copy(begin(Symbols), end(Symbols), std::ostream_iterator<double>(std::cout, "   "));
        // std::cout << "\n";

        // cout<<"Decoded \n";

        // copy(begin(decodedBits), end(decodedBits), std::ostream_iterator<double>(std::cout, "   "));
        // std::cout << "\n";


        cout<< "\n";

    }

    

    //     sourceBits = Source();

    //     Symbols = binaryToDecimalConversion(sourceBits, base);

    //     trans1 = SignalVectors1(Symbols, energyOfSymbols, M);
    //     trans2 = SignalVectors2(Symbols, energyOfSymbols, M);

    //     gnoise1 = GnoiseVector(0.0, stddev, SizeOfTransmission);
    //     gnoise2 = GnoiseVector(0.0, stddev, SizeOfTransmission);

    //     receivedBits1 = channelModel(trans1, gnoise1);
    //     receivedBits2 = channelModel(trans2, gnoise2);

    //     decodedBits = decisionBlock(receivedBits1, receivedBits2, M);

    //     errors = errorCount(Symbols, decodedBits);
    //     pe = errors/ Symbols.size();

    //     cout<<"Errors\n";

    //     cout<< errors<< "\n";
    //     cout<< pe << "\n";

    //     cout<< "\n";




    // cout<<"Symbols\n";
    // copy(begin(Symbols), end(Symbols), std::ostream_iterator<double>(std::cout, "   "));
    // std::cout << "\n";

    // cout<<"Decoded \n";

    // copy(begin(decodedBits), end(decodedBits), std::ostream_iterator<double>(std::cout, "   "));
    // std::cout << "\n";

    datafile(SNR,p_error);


    return 0;


}