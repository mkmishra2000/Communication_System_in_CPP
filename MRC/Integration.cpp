#include <iostream>
#include <vector>
#include <cmath>
#define pi 3.14179

using namespace std;

double DefiniteIntegration(vector<double> y_axis, 
                           double upperLimit, 
                           double lowerLimit, 
                           int numOfDiv)
{

    double width = (upperLimit-lowerLimit)/numOfDiv;

    double area, integration;
    integration=0;

    for(int k =0; k<numOfDiv; k++){
        area = y_axis[k]*width;

        integration = integration +area;
    }
    return integration;

}

int main()
{
    double a =0;
    double b =pi/2;
    double num =200;
    double width = (b-a)/num;
    double temp = a;

    vector<double> y_axisValues;
    vector<double> x_axisValues;

    for(int j=0; j<num; j++){
        x_axisValues.insert(x_axisValues.end(), temp);
        temp = temp + width;
    }

    for(int i =0; i<num; i++){
        y_axisValues.insert(y_axisValues.end(), sin(x_axisValues[i]));
    }

    double area, integration;
    integration=0;

    for(int k =0; k< num; k++){
        area = y_axisValues[k]*width;

        integration = integration +area;
    }

    cout<<integration<<endl;
    cout<<"My function"<<DefiniteIntegration(y_axisValues, b, a, num);
    return 0;
}