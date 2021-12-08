//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organization:- IIITDM KANCHEEPURAM
// Topic:- Encoder design for generating the codewords.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

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
vector<vector<int> > Encoder(vector<vector<int>> & MAT1, vector<vector<int>>&MAT2)
{
    if(MAT1[0].size()==MAT2.size()){
        return MatrixMultiGF2(MAT1, MAT2);
    }else{
        errorMSG();
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

vector<vector <int> > Tranpose_MAT(vector<vector <int> > MAT)
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



int main(){

    // messages vector.
    vector<vector<int>> MAT1{{1,0,0,1}};

    // Generator Matrix. 
    vector<vector<int>> GenMAT{{1,0,0,0,1,1,0},
                            {0,1,0,0,0,1,1}, 
                            {0,0,1,0,1,1,1},
                            {0,0,0,1,1,0,1}};

    // parity check matrix 
    vector<vector<int> > H_mat{{1,0,1,1,1,0,0},
                              {1,1,1,0,0,1,0},
                              {0,1,1,1,0,0,1}};
    
    cout<< "Message vector"<<endl;
    PrintMat(MAT1);

    cout<< endl;

    cout<< "Generator Matrix "<<endl;
    PrintMat(GenMAT);

    // Encoder Output
    // vector<vector<int>> ReMAT =Encoder(MAT1, MAT2);

    cout<< endl;

    cout<< "Final Output (codeward)"<<endl;
    vector<vector<int> > codewords;
    codewords = MatrixMultiGF2(MAT1, GenMAT);
    PrintMat(codewords);

    // codewords[0][1]+=1;

    vector<vector<int> > H_trans = Tranpose_MAT(H_mat);
    cout<<"Transpose of Parity check matrix "<<endl;
    PrintMat(H_trans);

    vector<vector<int> >Syndrome = MatrixMultiGF2(codewords, H_trans);

    cout<<"Syndrome values "<<endl;
    PrintMat(Syndrome);

    // PrintMat(ReMAT);

    cout<<endl;
    return 0;
}