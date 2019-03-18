#include <iostream>
#include <cmath>
#include <vector>

const double EPSILON = 10e-20; 
const double STEP = 10e-2;

double f(double x){ 
    return ((std::pow(x, 2) - 1) * std::pow(sinh(x), 3)); 
}

double falsi_f(double x1, double x2){ 
    return (f(x1) * x2 - f(x2) * x1) / (f(x1) - f(x2)); 
}

void falsi(double start_point, double end_point){ 
    double x1 = start_point; 
    double x2 = x1 + STEP; 
    std::vector<double> results; 

    while (x1 < end_point){ 
        if (f(x1) * f(x2) <= 0){ 
            double x3 = falsi_f(x1, x2);
            while (abs(f(x3)) > EPSILON){

                if(f(x1) * f(x3) <= 0) { 
                    x3 = falsi_f(x1, x3); 
                }

                if (f(x3) * f(x2) <= 0){ 
                    x3 = falsi_f(x3, x2); 
                }
            }   
            results.push_back(x3); 
        }
        x1 += STEP; 
        x2 += STEP; 
    }

    if(results.size() != 0){
        std::cout << "Roots of function f(x) = ((x^2 - 1) * (sinh x)^3): " << std::endl;
        for(int i = 0; i < results.size(); i++){ 
            std::cout << "x" << i << " = " << results[i] << " " << std::endl; 
        }
    } 
    else { 
        std::cout << "This non-liner equation has no roots." <<std::endl; 
    }

}


int main (int argc, char* argv[]){
    falsi(-100.0, 100.0); 
}
