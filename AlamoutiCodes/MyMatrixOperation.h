#include <iostream>
#include <vector>
#include <iterator>

using namespace std;


// Function for sum of two double numbers.
// Input is double numbers.
// Output is the sum of two double numbers.
double SumValue(double a, double b)
{
    return (a+b);
}


// Function for mul of two binary number in galois filed 2.
// Input binary numbers either 0 or 1.
// Output is the mul of two binary numbers according to galois 2.
double MulValue(double a, double b){

    return (a*b);

}


// Function of error message in matrices multiplications.
void errorMSG(){
    cout<<endl;
    cout<<"! Matrices size are not proper for multiplication !!! :-("<<endl;
    cout<<endl;
}

// Function for printing the matrix on console.
void PrintMat(vector<vector<double> > & MAT)
{
    for(int j =0; j<MAT.size();j++){
        for(int k =0; k<MAT[j].size(); k++){
            cout<<MAT[j][k]<< "   ";
        }
        cout<<endl;
    }
}


// Function for taking transpose of the given matrix. 
// Input is the matrix for some m*n dimension.
// Output is the matrix with n*m dimension having transpose of actual matrix 
vector<vector <double> > Transpose_MAT(vector<vector <double> > MAT)
{
    vector<vector<double> > TransMAT;
    vector <double> interMAT;

    for(int i =0; i<MAT[0].size(); i++){
        for(int j =0; j<MAT.size(); j++){
            interMAT.insert(interMAT.end(), MAT[j][i]);
        }
        TransMAT.push_back(interMAT);
        interMAT.clear();
    }

    return TransMAT;
}



// Function for multiply two matrices by galois filed 2.
// Input two matries.
// Output is the resultant Matrix if input vectors are proper otherwise zero matrix return.
vector<vector<double> >MatrixMulti(vector< vector<double> > &MAT1, vector< vector<double> > &MAT2)
{
    int Raws = MAT1.size(); //raws of the final vector.
    int Column = MAT2[0].size(); //Columns of the final vector.

    vector<double> finalRaw(Column, 0); //All zero entry in rows. 1*columns size.
    vector<vector<double>> resltantMAT(Raws, finalRaw); // raws*columns size.

    vector<double> inter1;
    vector<double> inter2;

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
                    resltantMAT[i][j] =SumValue(resltantMAT[i][j], MulValue(inter1[len], inter2[len]));
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


// Function to fix the issue of vector and matrixes
// Input is the vector signal 
// Output is the Matrix that contain vector as it's first raw
vector<vector<double>> convertVectorToMatrix(vector<double> Vec)
{
    vector<vector<double>> Mat;
    Mat.insert(Mat.end(), Vec);
    return Mat;
}


// Function for Multiplying two matrixs element by elements 
// Input are two matrix of same dimension 
// Output is another matrix that contain each element
//  as multiplication of each element.
vector<vector<double>> ElementWiseMultiplication(vector<vector<double>> Mat1, vector<vector<double>> Mat2)
{
    vector<vector<double>> MatResult;

    vector<double> multi;

    for(int i =0; i<Mat1.size();i++){
        for(int j=0; j<Mat1[0].size(); j++){
            multi.insert(multi.end(), Mat1[i][j]*Mat2[i][j]);
        }

        MatResult.insert(MatResult.end(), multi);
        multi.clear();
    }

    return MatResult;
}


// Function for Adding two matrixs element by elements 
// Input are two matrix of same dimension 
// Output is another matrix that contain each element
//  as Addition of both element.
vector<vector<double>> ElementWiseAddition(vector<vector<double>> Mat1, vector<vector<double>> Mat2)
{

    vector<vector<double>> MatResult;

    vector<double> Add;

    for(int i =0; i<Mat1.size();i++){
        for(int j=0; j<Mat1[0].size(); j++){
            Add.insert(Add.end(), Mat1[i][j]+Mat2[i][j]);
        }

        MatResult.insert(MatResult.end(), Add);
        Add.clear();
    }

    return MatResult;

}



// function for sum of two vectors 
// Inputs are two vectors 
// Output is the sum of two vectors.
vector<double> VectorSum(vector<double> Vec1, vector<double> Vec2)
{
    vector<double> SumResult;

    if(Vec1.size()==Vec2.size()){
        for(int k=0; k<Vec1.size(); k++){
            SumResult.insert(SumResult.end(), Vec1[k]+Vec2[k]);
        }

        return SumResult;
    }else{
        cout<<"Error !!!! Vector size should be same!!!"<<endl;
        return SumResult;
    }
}