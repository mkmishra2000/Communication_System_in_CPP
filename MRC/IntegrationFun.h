/////////////////////////////////////////////////////////////////////////////////////////////////
// Author:- MANAS KUMAR MISHRA
// Organisation:- IIITDM KANCHEEPURAM
// Topic:- Definite integration [a, b]
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <cmath>

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
