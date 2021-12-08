//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organization:- IIITDM KANCHEEPURAM
// Topic:- Hamming distance based decoding.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <time.h>
#include <fstream>
#include "BPSK_BSC.h"

#define one_million 1000000

using namespace std;

// Class that contains all main components related to hamming distance
// based decoder 
class hamming
{
    public:

     // All possible messages vector.
    vector<vector<int>> Messages{{0,0,0,0},
                             {0,0,0,1},
                             {0,0,1,0},
                             {0,0,1,1},
                             {0,1,0,0},
                             {0,1,0,1},
                             {0,1,1,0},
                             {0,1,1,1},
                             {1,0,0,0},
                             {1,0,0,1},
                             {1,0,1,0},
                             {1,0,1,1},
                             {1,1,0,0},
                             {1,1,0,1},
                             {1,1,1,0},
                             {1,1,1,1}};
    // Generator Matrix. 
    vector<vector<int>> GenMAT{{1,0,0,0,1,1,0},
                            {0,1,0,0,0,1,1}, 
                            {0,0,1,0,1,1,1},
                            {0,0,0,1,1,0,1}};
    
    vector<int>SourceVector();
    
    vector<vector<int> >Encoder(vector<vector<int>> & MAT1, 
                                vector<vector<int>>&MAT2);
    
    vector<int> EncodedBits(vector<int> SourceBits, 
                            vector<vector<int> >GenMAT, 
                            const int K);
    
    vector<vector <int> > AllCodewordsList(vector<vector<int> > &MAT1,
                                        vector<vector<int> > &MAT2);
    
    vector<int> HammingDecoder(vector<double> decodedBits,
                             vector<vector<int> > MSGList,
                             vector<vector<int> > CodeList, 
                             const int n);
    
    double errorCalculation (vector<int> sourceBits, vector<int> decodedBits);
};


// Function for generating binary bits at source side. Each bit is equiprobable.
// Input is nothing
// Output is a vector that contains the binary bits of length one_million*1.
// one_million is defined by programmer of this code.
vector<int>hamming::SourceVector()
{
    vector<int> sourceBits;

    // Use current time as seed for random generator
    srand(time(0));
 
    for(int i = 0; i<one_million; i++){
        sourceBits.insert(sourceBits.end(), rand()%2);
    }

    return sourceBits;
}


// Function for sum of two binary number in galois filed.
// Input binary numbers either 0 or 1.
// Output is the sum of two binary numbers according to galois 2.
// 0 0 | 0
// 0 1 | 1
// 1 0 | 1
// 1 1 | 0 
int Galois2Sum(int a, int b){

    return (~a & b)| (a & ~b);
}


// Function for mul of two binary number in galois filed 2.
// Input binary numbers either 0 or 1.
// Output is the mul of two binary numbers according to galois 2.
// 0 0 | 0
// 0 1 | 0
// 1 0 | 0
// 1 1 | 1 
int Galois2Mul(int a, int b){

    return (a & b);

}


// Function of error message in matrices multiplications.
void errorMSG(){
    cout<<endl;
    cout<<"! Matrices size are not proper for multiplication !!! :-("<<endl;
    cout<<endl;
}


// Function for multiply two matrices by galois filed 2.
// Input two matries.
// Output is the resultant Matrix if input vectors are proper otherwise zero matrix return.
vector<vector<int> >MatrixMultiGF2(vector< vector<int> > &MAT1, vector< vector<int> > &MAT2)
{
    int Raws = MAT1.size(); //raws of the final vector.
    int Column = MAT2[0].size(); //Columns of the final vector.

    vector<int> finalRaw(Column, 0); //All zero entry in rows. 1*columns size.
    vector<vector<int>> resltantMAT(Raws, finalRaw); // raws*columns size.

    vector<int> inter1;
    vector<int> inter2;

    for(int i =0; i<resltantMAT.size(); i++){
        for(int j =0; j<resltantMAT[0].size(); j++){

            for(int col =0; col<MAT1[0].size(); col++){
                inter1.insert(inter1.end(), MAT1[i][col]);
            }
            for(int row = 0; row<MAT2.size(); row++){
                inter2.insert(inter2.end(), MAT2[row][j]);
            }

            // Columns of first matrixs should be same as raws of second matrix.
            if(inter1.size() == inter2.size()){
                for(int len =0; len<inter2.size(); len++){
                    resltantMAT[i][j] =Galois2Sum(resltantMAT[i][j], Galois2Mul(inter1[len], inter2[len]));
                }
            }else{
                errorMSG();
            }

            inter1.clear(); //Clear out intermediate vector
            inter2.clear(); //Clear out intermediate vector
            
        }
    }

    return resltantMAT;
}


// Function for encoding the messages
// Input is the code and generator matrix.
// Output is the encoded codeward.
vector<vector<int> >hamming::Encoder(vector<vector<int>> & MAT1, vector<vector<int>>&MAT2)
{
    if(MAT1[0].size()==MAT2.size()){
        return MatrixMultiGF2(MAT1, MAT2);
    }else{
        errorMSG();
        return {};
    }
    
}


//Function for encoding the source stream.
//Inputs are the Ssource bits, generator matrix, and K (Block size).
//Output is the encoded bits.
//It is taking k size block from source bit stream and doing encoding on it.
//encoding implies u*G=v. 
vector<int> hamming::EncodedBits(vector<int> SourceBits, 
                                vector<vector<int> > GenMAT, 
                                const int k)
{
    vector<int>inter; //for taking blocks of code from sourcebits.
    vector<int> encoded;//For final encoded bits.
    vector<vector<int> > MAT1; //For message bits.
    vector<vector<int> > MAT2; //For output of encoder.

    // If length of source bits is divisible by k (block size)
    if(SourceBits.size()%k == 0){
        //Encoding start (blocks of message).
        for(int i =0; i<SourceBits.size(); i = i+k){
            //generate the blocks of message.
            for(int j =i; j<k+i; j++){
                inter.insert(inter.end(), SourceBits[j]);
            }

            MAT1.push_back(inter);
            MAT2 = Encoder(MAT1, GenMAT); //Encoding operation.

            // Combining the all encoded bits into single vector 
            for(int i =0; i<MAT2[0].size(); i++){
                encoded.insert(encoded.end(), MAT2[0][i]);
            }

            //clear all extra vectors and matrices.
            MAT1.clear();
            MAT2.clear();
            inter.clear();
        }

        return encoded;
    }
    else{ //if source bits size is not divisible by k (block size).
        // For this case, do padding of extra zeros in the bit stream 
        int padding = SourceBits.size()%k;

        for(int i =0; i<k-padding; i++){
            SourceBits.insert(SourceBits.end(), 0);
        }

       //Encoding start (blocks of message).
        for(int i =0; i<SourceBits.size(); i = i+k){
            
            //generate the blocks of message.
            for(int j =i; j<k+i; j++){
                inter.insert(inter.end(), SourceBits[j]);
            }

            MAT1.push_back(inter);
            MAT2 = Encoder(MAT1, GenMAT); //Encoding operation.

            // Combining the all encoded bits into single vector 
            for(int i =0; i<MAT2[0].size(); i++){
                encoded.insert(encoded.end(), MAT2[0][i]);
            }

            //clear all extra vectors and matrices.
            MAT1.clear();
            MAT2.clear();
            inter.clear();
        }

        return encoded;
    }
}


// Function for generating all possible codewords for estimating the hamming distance
// Inputs are all possible message blocks and generator matrix.
// Output is the matrix of size (2^K * n) having all possible codewords as a row.. 
// Where k is block size of message and n is codeword length.
vector<vector <int> >hamming::AllCodewordsList(vector<vector<int> > &MAT1,
                                        vector<vector<int> > &MAT2)
{
    vector<vector<int> > CodewordList;
    vector<vector<int> > Inter;

    for(int i =0; i<MAT1.size(); i++){
        // Take raws of message matrix one by one. 
        Inter.push_back(MAT1[i]);

        // Calculate the codeword for the message
        Inter = MatrixMultiGF2(Inter, MAT2);

        // Put everycode into the codeword list
        CodewordList.push_back(Inter[0]);

        Inter.clear();
    }

    return CodewordList;
}


// Function of error, if vector length in hamming distance is not same. 
void hamerrorMSG()
{
    cout<<"Error vector length in not same... NO hamming distance"<<endl;
    cout<<endl;
}


// Function for printing the vector on the console output of int type
void PrintVectorInt(vector<int> vectr)
{
    std::copy(begin(vectr), end(vectr),  std::ostream_iterator<int>(std::cout, " "));
    cout<<endl;
}


// Function for hamming distance between two codewords.
// Inputs are two vectors V1 abnd V2.
// Output is the hamming distance as integer. 
// If codeword bits are same count zero, otherwise one.
int hamDistance(vector<int> v1, vector<int> v2)
{
    int hamDis=0;

    if(v1.size()==v2.size()){
        for(int i =0; i<v1.size(); i++){
            hamDis += (~v1[i] & v2[i])|(v1[i] & ~v2[i]); //XOR operation((~a & b)|(a & ~b))
        }

        return hamDis;
    }else{
        hamerrorMSG();
        return -1;
    }   
}


// Function for find minimum entry index in any array.
// Input is the array.
// Output is the index that contain minimum element.  
int minEntryIndex(vector<int> Array)
{
    int minIndex =0;

    for(int i =0; i<Array.size(); i++){
        if(Array[i]< Array[minIndex]){
            minIndex = i;
        }
    }

    return minIndex;
}


// Function for calculating the hamming distance between received codeword and set of possible codewords.
// After that calculating the minimum hamming distance amoung all hamming distances.
// Inputs are the received codeword and set of all codewords.
// Output is the index that has minimum hamming distance. That
// index also indicate the close messages index. 
int minHammingIndex(vector<int> vec, vector<vector<int> > CodeList)
{
    vector<int> hammingDistList;
    vector<int> interblock;

    for(int i =0; i<CodeList.size(); i++){
        for(int j =0 ; j<CodeList[i].size(); j++){
            interblock.insert(interblock.end(), CodeList[i][j]);
        }
        
        hammingDistList.insert(hammingDistList.end(), hamDistance(vec, interblock));
        interblock.clear();
    }

    int minIndex = minEntryIndex(hammingDistList);

    return minIndex;
}


// Function for decoding the received vector from BSC. It uses hamming distance methods.
// Inputs are decodedBits from BSC, all possible messages and codewords in same order,
// and size of codewords.
// Outputs are the decided message bits. 
vector<int>hamming::HammingDecoder(vector<double> decodedBits,
                             vector<vector<int> > MSGList,
                             vector<vector<int> > CodeList, 
                             const int n)
{
    vector<int> finalCode;
    vector<int> interblock;
    int index;

    for(int i=0; i<decodedBits.size(); i=i+n){
        for(int j =i; j<i+n; j++){
            interblock.insert(interblock.end(), decodedBits[j]);
        }
        index = minHammingIndex(interblock, CodeList);
        interblock.clear();
        interblock = MSGList[index];
        for(int len =0; len< interblock.size(); len++){
            finalCode.insert(finalCode.end(), interblock[len]);
        }
        interblock.clear();
    }

    return finalCode;

}


// Function for remove the padding of zero 
// Input are array and actual size without padding.
// Output is the array with removed zero.
vector<int> finalCut(vector<int> Array, int sizeofArray)
{
    int cut = Array.size()-sizeofArray;

    for(int i =0; i<cut; i++){
        Array.pop_back();
    }

    return Array;
}


// Function for printing the vector on the console output.
void PrintVectorDouble(vector<double> vectr)
{
    std::copy(begin(vectr), end(vectr),  std::ostream_iterator<double>(std::cout, " "));
    cout<<endl;
}


// Function to count number of errors in the received bits.
// Inputs are the sourcebits and decodedbits
// OUtput is the number of error in received bits.
// error: if sourcebit  != receivebit
double hamming::errorCalculation (vector<int> sourceBits, vector<int> decodedBits)
{
    double countError =0;
    for(int i =0; i<sourceBits.size();i++){
        if(sourceBits[i]!= decodedBits[i]){
            countError++;
        }
    }

    return countError;
}


// Function to count number of errors in the received bits without hamming distance decoder.
// Inputs are the sourcebits and decodedbits
// OUtput is the number of error in received bits.
// error: if sourcebit  != receivebit
double errorCalculationcheck(vector<int> sourceBits, vector<double> decodedBits)
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
void datafile(vector<double> xindB, vector<double> Prob_error)
{
    ofstream outfile;

    outfile.open("Ham1.dat");

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

    // Object of class 
    hamming ham;

    // source definition
    vector<int> sourceBits;

    // Encoded bits
    vector<int> encodedBits;

    // ALL possible codewords (Codeword table)
    vector <vector<int> > CodewordTable;

    // Transmitted bits
    vector<double> Trans;

    // Noise definition
    vector<double> gnoise;

    // Receive bits
    vector<double> recevBits;

    // Decision block
    vector<double> decodedBits;

    //Hamming distance decoder output
    vector<int> finalBits;

    // Energy values 
    vector<double> energyOfBits;

    // Probability of error vector 
    vector<double> Prob_error;
    
    double N_o =4;
    double pe, stdnoise, errorCount, normalValue;
    int k =4;
    int n =7;

    stdnoise = sqrt(N_o/2); // std of noise.

    // SNR in dB 
    vector<double> SNR_dB;

    for(float i =0; i<=10; i=i+0.125)
    {
        SNR_dB.insert(SNR_dB.end(), i);
    }

    for(int i =0; i<SNR_dB.size(); i++){

        normalValue = pow(10, (SNR_dB[i]/10));
        energyOfBits.insert(energyOfBits.end(), N_o*normalValue);
    }

    CodewordTable =ham.AllCodewordsList(ham.Messages, ham.GenMAT); //All possible codewords.

    for(int step=0; step<energyOfBits.size(); step++){
        
        sourceBits = ham.SourceVector();

        encodedBits = ham.EncodedBits(sourceBits, ham.GenMAT, k);

        Trans = bit_maps_to_symbol_of_energy_E(encodedBits, energyOfBits[step]);

        gnoise = GnoiseVector(0.0, stdnoise, encodedBits.size());

        recevBits = receiveBits(Trans, gnoise);

        decodedBits = decisionBlock(recevBits);

        finalBits = ham.HammingDecoder(decodedBits, ham.Messages, CodewordTable, n);

        if(finalBits.size()!= sourceBits.size()){
            finalBits= finalCut(finalBits, sourceBits.size()); 
        }

        errorCount = errorCalculationcheck(encodedBits, decodedBits);
        cout<<"Error without hamming distance "<<errorCount<<endl;

        errorCount = ham.errorCalculation(sourceBits, finalBits);
        pe = errorCount/(one_million);

        cout<<"Errors : "<<errorCount<<endl;
        cout<<"Pe     : "<<pe<<endl;
        cout<<endl;

        Prob_error.insert(Prob_error.end(), pe);
        cout<<endl;

    }
    

    datafile(SNR_dB, Prob_error);
    cout<<endl;
    return 0;
}