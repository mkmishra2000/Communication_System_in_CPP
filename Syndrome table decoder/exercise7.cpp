//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organization:- IIITDM KANCHEEPURAM
// Topic:- Syndrome table decoding.
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

class SyndromeDecoding
{
    public:

    // Generator Matrix. 
    vector<vector<int>> GenMAT{{1,0,0,0,1,1,0},
                               {0,1,0,0,0,1,1}, 
                               {0,0,1,0,1,1,1},
                               {0,0,0,1,1,0,1}};

    // parity check matrix 
    vector<vector<int> > H_mat{{1,0,1,1,1,0,0},
                               {1,1,1,0,0,1,0},
                               {0,1,1,1,0,0,1}};

    // syndrome 
    vector<vector<int> > Syndromes{{0,0,0},
                                   {0,0,1},
                                   {0,1,0},
                                   {1,0,0},
                                   {0,1,1},
                                   {1,1,0},
                                   {1,1,1},
                                   {1,0,1}};
    // Coset ladder 
    vector<vector <int> > ErrorTable{{0,0,0,0,0,0,0},
                                     {0,0,0,0,0,0,1},
                                     {0,0,0,0,0,1,0},
                                     {0,0,0,0,1,0,0},
                                     {0,0,0,1,0,0,0},
                                     {0,0,1,0,0,0,0},
                                     {0,1,0,0,0,0,0},
                                     {1,0,0,0,0,0,0}};
    
    vector<int>SourceVector();

    vector<vector<int> >Encoder(vector<vector<int>> & MAT1, 
                                vector<vector<int>>&MAT2);
    
    vector<int> EncodedBits(vector<int> SourceBits, 
                            vector<vector<int> >GenMAT, 
                            const int K);
    
    vector<int> SyndromDecoder(vector<double> demodBits,
                                              vector<vector<int> > H_mat,
                                              vector<vector<int> > SyndromeTable,
                                              vector<vector<int> > errorPattern,
                                              const int &n);

    double errorCalculation (vector<int> sourceBits, vector<int> decodedBits);

};



// Function for generating binary bits at source side. Each bit is equiprobable.
// Input is nothing
// Output is a vector that contains the binary bits of length one_million*1.
// one_million is defined by programmer of this code.
vector<int>SyndromeDecoding::SourceVector()
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


// Function for taking transpose of the given matrix. 
// Input is the matrix for some m*n dimension.
// Output is the matrix with n*m dimension having transpose of actual matrix 
vector<vector <int> > Transpose_MAT(vector<vector <int> > MAT)
{
    vector<vector<int> > TransMAT;
    vector <int> interMAT;

    for(int i =0; i<MAT[0].size(); i++){
        for(int j =0; j<MAT.size(); j++){
            interMAT.insert(interMAT.end(), MAT[j][i]);
        }
        TransMAT.push_back(interMAT);
        interMAT.clear();
    }

    return TransMAT;
}


// Function for encoding the messages
// Input is the code and generator matrix.
// Output is the encoded codeward.
vector<vector<int> >SyndromeDecoding::Encoder(vector<vector<int>> & MAT1, vector<vector<int>>&MAT2)
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
vector<int> SyndromeDecoding::EncodedBits(vector<int> SourceBits, 
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


// Function for add two codewords according to galois field 2.
// Inputs are two vectors in GF2. 
// Output is sum of two vector (Output type <int>). 
vector<int> AddVectorsinGF2(vector<int> A,
                       vector<int> B)
{

    if(A.size()== B.size()){
        vector<int> Add;
        for(int mem =0; mem<A.size(); mem++){
            Add.insert(Add.end(), Galois2Sum(A[mem], B[mem]));
        }

        return Add;
    }else{
        cout<<"Error vector length is not same !!!"<<endl;
        return {};
    }

}


// Function for printing the matrix on console.
void PrintMat(vector<vector<int> > & MAT)
{
    for(int j =0; j<MAT.size();j++){
        for(int k =0; k<MAT[j].size(); k++){
            cout<<MAT[j][k]<< "   ";
        }
        cout<<endl;
    }
}


// Function for decoding the received bits from BSC channel.
// Inputs are following 1.Bits from BSC 2. Parity check matrix, 3. Symdrome list,
// 4. Error pattern list (coset ladder), 5. BLock size.
// Output is the decoded bits   
vector<int> SyndromeDecoding::SyndromDecoder(vector<double> demodBits,
                                              vector<vector<int> > H_mat,
                                              vector<vector<int> > SyndromeTable,
                                              vector<vector<int> > errorPattern, 
                                              const int &n)
{
    vector<int> decodedVector;
    vector<vector <int> > interMat;
    vector<int> block;
    vector<int> error;
    vector<vector<int> > Trans=Transpose_MAT(H_mat);

    for(int i =0; i<demodBits.size(); i=i+n){

        // Extract block of the codewords.
        for(int j =i; j<n+i; j++){
            block.insert(block.end(), demodBits[j]);
            
        }

        interMat.push_back(block);

        // Calculate the syndrome for that block.
        interMat = MatrixMultiGF2(interMat, Trans);

        // Error patterns
        for(int i =0; i<SyndromeTable.size();i++){

            if(interMat[0] == SyndromeTable[i]){
                
                error = errorPattern[i];
            }
        }

        // Corrected codewords. 
        block = AddVectorsinGF2(block, error);

        for(int mem =0; mem<block.size(); mem++){
            decodedVector.insert(decodedVector.end(), block[mem]);
        }
        
        block.clear();
        interMat.clear();
        error.clear();

    }

    return decodedVector;
}


// Function to count number of errors in the received bits.
// Inputs are the sourcebits and decodedbits
// OUtput is the number of error in received bits.
// error: if sourcebit  != receivebit
double SyndromeDecoding::errorCalculation (vector<int> sourceBits, vector<int> decodedBits)
{
    double countError =0;
    for(int i =0; i<sourceBits.size();i++){
        if(sourceBits[i]!= decodedBits[i]){
            countError++;
        }
    }

    return countError;
}


// Function for printing the vector on the console output.
void PrintVectorDouble(vector<double> vectr)
{
    std::copy(begin(vectr), end(vectr),  std::ostream_iterator<double>(std::cout, " "));
    cout<<endl;
}


// Function for printing the vector on the console output of int type
void PrintVectorInt(vector<int> vectr)
{
    std::copy(begin(vectr), end(vectr),  std::ostream_iterator<int>(std::cout, " "));
    cout<<endl;
}


// Function to store the data in the file (.dat)
// Input is the SNR per bit in dB and calculated probability of error
// Output is the nothing but in processing it is creating a file and writing data into it.
void datafile(vector<double> xindB, vector<double> Prob_error)
{
    ofstream outfile;

    outfile.open("Synd1.dat");

    if(!outfile.is_open()){
        cout<<"File opening error !!!"<<endl;
        return;
    }

    for(int i =0; i<xindB.size(); i++){
        outfile<< xindB[i] << " "<<"\t"<<" "<< Prob_error[i]<< endl;
    }

    outfile.close();
}



int main()
{
    // Object of the class 
    SyndromeDecoding syn1;

    // source definition
    vector<int> sourceBits;

    // Encoded bits
    vector<int> encodedBits;

    // Transmitted bits
    vector<double> Transmitted;

    // Noise definition
    vector<double> gnoise;

    // Receive bits
    vector<double> recevBits;

    // Decision block
    vector<double> decodedBits;

    // Syndrome decoder. 
    vector<int> finalout;

    // Energy values 
    vector<double> energyOfBits;

    // Probability of error vector 
    vector<double> Prob_error;
    

    int k =4, n=7;
    double energy = 11;
    double N_o =4;
    double pe, stdnoise, errorCount, normalValue;

    stdnoise = sqrt(N_o/2); // std of noise.


    // SNR in dB 
    vector<double> SNR_dB;

    for(float i =0; i<=12; i=i+0.125)
    {
        SNR_dB.insert(SNR_dB.end(), i);
    }

    for(int i =0; i<SNR_dB.size(); i++){

        normalValue = pow(10, (SNR_dB[i]/10));
        energyOfBits.insert(energyOfBits.end(), N_o*normalValue);
    }

    for(int step =0; step<energyOfBits.size(); step++){

        sourceBits = syn1.SourceVector();

        encodedBits = syn1.EncodedBits(sourceBits, syn1.GenMAT, k);

        ///////////   BSC    /////////

        Transmitted = bit_maps_to_symbol_of_energy_E(encodedBits, energyOfBits[step]);


        gnoise = GnoiseVector(0.0,stdnoise, encodedBits.size());


        recevBits = receiveBits(Transmitted, gnoise);


        decodedBits = decisionBlock(recevBits);

        ///////////   BSC    /////////

        finalout = syn1.SyndromDecoder(decodedBits, syn1.H_mat, syn1.Syndromes, syn1.ErrorTable, n);


        errorCount = syn1.errorCalculation(encodedBits, finalout);

        pe = errorCount / (encodedBits.size());

        cout<<" Error : "<< errorCount<<endl;
        cout<<" Pe    : "<< pe<< endl;

        Prob_error.insert(Prob_error.end(), pe);
        cout<<endl;
    }
    // sourceBits = syn1.SourceVector();

    // cout<<"Source bits  "<<endl;
    // PrintVectorInt(sourceBits);

    // encodedBits = syn1.EncodedBits(sourceBits, syn1.GenMAT, k);

    // cout<<"Encoded bits "<<endl;
    // PrintVectorInt(encodedBits);

    // ///////////   BSC    /////////

    // Transmitted = bit_maps_to_symbol_of_energy_E(encodedBits, energy);
    // cout<<"Transmitted symbols "<<endl;
    // PrintVectorDouble(Transmitted);

    // gnoise = GnoiseVector(0.0,stdnoise, encodedBits.size());
    // cout<<"Gnoise "<<endl;
    // PrintVectorDouble(gnoise);

    // recevBits = receiveBits(Transmitted, gnoise);
    // cout<<"Received symbols "<<endl;
    // PrintVectorDouble(recevBits);

    // decodedBits = decisionBlock(recevBits);
    // cout<<"Decoded symbols"<<endl;
    // PrintVectorDouble(decodedBits);

    // ///////////   BSC    /////////

    // finalout = syn1.SyndromDecoder(decodedBits, syn1.H_mat, syn1.Syndromes, syn1.ErrorTable, n);
    // cout<<"Final symbols "<<endl;
    // PrintVectorInt(finalout);

    // errorCount = syn1.errorCalculation(encodedBits, finalout);

    // pe = errorCount / (encodedBits.size());

    // cout<<" Error : "<< errorCount<<endl;
    // cout<<" Pe    : "<< pe<< endl;

    datafile(SNR_dB, Prob_error);
    cout<<endl;


    return 0;
}
