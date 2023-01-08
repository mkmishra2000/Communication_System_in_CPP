#include<iostream>

using namespace std;


// function for factorial of n
// Input is a integer number greater than 0
// Output is the factorial result 
double factorial(int Num)
{
    if(Num==1){
        return 1;
    }
    else if(Num ==0){
        return 1;
    }
    else if(Num<0){
        cout<< "Worng Output"<<endl;
        return 0;
    }
    else{
        return (Num*factorial(Num-1));
    }
}


// Function for the combination l_C_k
// Inputs are the two integer numbers l and k 
// Output is the answer of the l choose K.
double l_Choose_K(int l, int k)
{
    double ans;
    if(l<k){
        cout<<"L should not be less than k"<<endl;
        return 0;
    }
    else{
        ans = factorial(l)/(factorial(k)*factorial(l-k));
        return ans;
    }
}


int main(){

    int l=12;
    int k =2;

    // cout<<factorial(l)<<endl;
    cout<<l_Choose_K(l, k)<<endl;
    return 0;
}