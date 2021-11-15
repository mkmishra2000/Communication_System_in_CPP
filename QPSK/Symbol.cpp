#include<iostream>
#include <math.h>
#include <iterator>
#include <random>
#include <chrono>
#include <time.h>

using namespace std;

#define oneMillion 10

vector<double> sourceSymbols()
{
    vector<double> messageInSymbols;

    srand(time(0));

    for(int i=0; i<oneMillion; i++){
        messageInSymbols.insert(messageInSymbols.end(), rand()%2);
    }

    return messageInSymbols;

}

vector <double> binaryToDecimalCOnversion(vector<double> sourceBits, int base)
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


int main()
{

    vector<double> sourceBits = sourceSymbols();
    copy(begin(sourceBits), end(sourceBits), std::ostream_iterator<double>(std::cout, "   "));

    cout<<endl;

    vector<double> convertedBits = binaryToDecimalCOnversion(sourceBits, 2);

    copy(begin(convertedBits), end(convertedBits), std::ostream_iterator<double>(std::cout, "   "));
    cout<<endl;

    cout << convertedBits.size();


    return 0;
}