// QPSK (Quadrature phase shift keying)


#include <iostream>
#include <math.h>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>

using namespace std;

#define oneMillion 1000000

// Function to generate message symbols (0, 1, 2, 3)
// Input is nothing  
// Output is the vector containing symbols 
vector<double> source_bits()
{
    vector<double> messageInSymbols;

    srand(time(0));

    for(int i=0; i<oneMillion; i++){
        messageInSymbols.insert(messageInSymbols.end(), rand()%2);
    }

    return messageInSymbols;

}

vector <double> binaryToDecimalConversion(vector<double> sourceBits, int base)
{
    vector <double> convertedBits;

    if(sourceBits.size()% base == 0 ){

        int finalSize = sourceBits.size()/base;
        int start = 0;
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
    else{
        int addedBitsNO =  base - (sourceBits.size()%base);

        for(int q=0;q<addedBitsNO;q++){
            sourceBits.insert(sourceBits.end(), 0);
        }

        int finalSize = sourceBits.size()/base;
        int start = 0;
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
vector<double> SignalVectors1(vector<double> sourceSymbols, double eneryOfSymbols)
{
    vector<double> y1;

    for(int i=0; i<sourceSymbols.size(); i++){
        if(sourceSymbols[i]==0)
        {
            y1.insert(y1.end(), sqrt(eneryOfSymbols));
        }
        else if (sourceSymbols[i]==1)
        {
            y1.insert(y1.end(), sqrt(eneryOfSymbols));
        }
        else if (sourceSymbols[i]==2)
        {
            y1.insert(y1.end(), -sqrt(eneryOfSymbols));
        }
        else
        {
            y1.insert(y1.end(), -sqrt(eneryOfSymbols));
        }
        
    }

    return y1;
}
vector<double> SignalVectors2(vector<double> sourceSymbols, double eneryOfSymbols)
{
    vector<double> y2;

    for(int i=0; i<sourceSymbols.size(); i++){
        if(sourceSymbols[i]==0)
        {
            y2.insert(y2.end(), sqrt(eneryOfSymbols));
        }
        else if (sourceSymbols[i]==1)
        {
            y2.insert(y2.end(), -sqrt(eneryOfSymbols));
        }
        else if (sourceSymbols[i]==2)
        {
            y2.insert(y2.end(), -sqrt(eneryOfSymbols));
        }
        else
        {
            y2.insert(y2.end(), sqrt(eneryOfSymbols));
        }
        
    }

    return y2;
}


// Function for generating random noise based on gaussian distribution N(mean, variance).
// Input mean and standard deviation.
// Output is the vector that contain gaussian noise as an element.
vector<double> GnoiseVector(double mean, double stddev)
{
    std::vector<double> noise;
    
    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    std::normal_distribution<double> dist(mean, stddev);

    // Add Gaussian noise
    for (int i =0; i<oneMillion; i++) {
        noise.insert(noise.end(),dist(generator));
    }

    return noise;
}


// Function for model the received symbols
// Input sourcesymbols and AWGN
// Output is a vector that represent the addition of two vectors
vector<double> receiveBits(vector<double> sourceSymbols, vector<double> AWGnoise)
{
    vector<double> received;

    for(int i=0; i<sourceSymbols.size(); i++){
        received.insert(received.end(), sourceSymbols[i]+AWGnoise[i]);
    }

    return received;
}


// Function for take decision after receiving the vector
// Input is the 2D receiver vector
// Output is the decision {0, 1, 2, 3} 
vector <double> decisionBlock(vector<double> receive1, vector<double> receive2)
{
    vector<double> decision;

    for(int i=0; i<receive1.size(); i++){
        if(receive1[i]>= 0 && receive2[i]>= 0){
            decision.insert(decision.end(), 0);

        }
        else if (receive1[i]>= 0 && receive2[i]<0)
        {
            decision.insert(decision.end(), 1);
        }
        else if (receive1[i] <0 && receive2[i]<0)
        {
            decision.insert(decision.end(), 2);
        }
        else{
            decision.insert(decision.end(), 3);
        }
        
    }

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


int main()
{

    vector<double> SourceBits;
    SourceBits = source_bits();

    vector<double> Symbols;
    int Base =2;
    Symbols = binaryToDecimalConversion(SourceBits, Base); 

    double Energy = 152;
    


    vector<double> Y1=SignalVectors1(Symbols, Energy);
    vector<double> Y2=SignalVectors2(Symbols, Energy) ;


    // Noise definition
    vector<double> gnoise1 = GnoiseVector(0.0, 2);
    vector<double> gnoise2 = GnoiseVector(0.0, 2);

    vector<double> receivedBits1 = receiveBits(Y1, gnoise1);
    vector<double> receivedBits2 = receiveBits(Y2, gnoise2);

    vector<double> decision = decisionBlock(receivedBits1, receivedBits2);

    int errors = errorCount(Symbols, decision); 

    // copy(begin(decision), end(decision), std::ostream_iterator<double>(std::cout, "   "));
    std::cout << "\n";

    cout<< errors;
    cout<<" \n";

    return 0;
}
