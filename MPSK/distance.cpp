#include<iostream>
#include<cmath>
#include <random>
#include <chrono>
#include <time.h>
#include <bits/stdc++.h>

#define one_million 20

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
        value = (rand()%4);
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
        value = (rand()%4);
        sourceBits.insert(sourceBits.end(), value);
    }

    return sourceBits;
}

int minElement_indexReturn (vector <double> distance)
{
    // double minValue = *min_element(distance.begin(), distance.end()); 
    int minElementIndex = min_element(distance.begin(),distance.end()) - distance.begin();

    return minElementIndex;
}

int main()
{
    vector <double> x;
    vector <double> y;

    x = Source2();
    y = Source();

    vector <double> result;

    copy(begin(x), end(x), std::ostream_iterator<double>(std::cout, "   "));

    cout<<endl;

    cout<<minElement_indexReturn(x)<<endl;
    

    return 0;
}