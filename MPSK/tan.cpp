#include<iostream>
#include<cmath>
#include <random>
#include <chrono>
#include <time.h>

#define pi 2*acos(0.0)
#define one_million 10

using namespace std;



// Function for generating binary bits at source side. Each bit is equiprobable.
// Input is nothing
// Output is a vector that contains the binary bits of length one_million*1.
vector<double> Source()
{
    vector<double> sourceBits;
    double value;

    // Use current time as seed for random generator
    srand(time(0));
 
    for(int i = 0; i<one_million; i++){
        value = (rand()%4) -2;
        sourceBits.insert(sourceBits.end(), value);
    }

    return sourceBits;
}

vector<double> Source2()
{
    vector<double> sourceBits;
    double value;

    // Use current time as seed for random generator
    srand(time(0)+2);
 
    for(int i = 0; i<one_million; i++){
        value = (rand()%4)-1;
        sourceBits.insert(sourceBits.end(), value);
    }

    return sourceBits;
}


int main(){

    vector <double> x;
    vector <double> y;

    x = Source2();
    y = Source();

    vector <double> result;

    for(int i =0; i<x.size(); i++ ){
        cout << x[i] << " "<< y[i] << endl;


        if(y[i]==0 && x[i]==0){
             result.insert(result.end(), 0); 
        }
        else{
             result.insert(result.end(), atan2(y[i], x[i])); 
             
        }

        result[i] = result[i]*(180/3.14159);

        if(y[i]<0){
            result[i] = result[i] + 360;
        }

        cout<< result[i]<< "degree" << endl;
    }


    vector<double> angles;
    int M =4;
    double values;
    for(int k =0; k<M; k++)
    {
        values = (2*k+1)*(180/M);
        angles.insert(angles.end(), values);

        cout<<angles[k]<<" ";
    }

    cout<< endl;

    for( int index =0; index < result.size(); index++){
        for(int j =0; j<M-1; j++){
            if(angles[j]<result[index] && angles[j+1]>= result[index]){
                cout << j+1 <<" ";
            }
        }
        if(result[index] <= angles[0] || result[index]>angles[M-1]){
            cout<< 0 <<" ";
        } 
    }
    return 0;

}
