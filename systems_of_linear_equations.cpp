#include <iostream>
#include <string>
#include <cmath>
#include <set>
using std::string;
const double EPSILON = 10e-8;

class Equation { 
    public: 
        Equation(); 
        void results(double start, double end, double x_start, double y_start);
    private:
        double step;
        typedef std::pair<double, double> pairs;
        std::set<pairs> roots;    
        double sigma(double x, double y, string coor);
        void Newton_method(double x, double y);  
};

Equation::Equation(){
    this->step = 0.01; 
}

double Equation::sigma(double x, double y, string coor){
    if(coor == "x") return (-2 + pow(x, 2) + 8 * y - 8 * x * y + 2 * pow(x, 2) *  y - 2 * pow(y, 2))/(2 * (x - 4 * y + 2 * x * y));
    if (coor == "y") return ((2 + pow(x, 2) - 2 * pow(y, 2) + x * (-3 + y + pow(y, 2)))/(x - 4 * y + 2 * x * y));
    return -1; 
}

void Equation::Newton_method(double x, double y){ 
    double x_start =  - (std::numeric_limits<double>::infinity());
    double y_start =  - (std::numeric_limits<double>::infinity());

    while (abs(x - x_start) > EPSILON || abs(y - y_start) > EPSILON) {
        x_start = x;
        y_start = y;
        x -= sigma(x_start, y_start, "x");
        y -= sigma(x_start, y_start, "y");
    }
    roots.insert(std::make_pair(((int)(x / EPSILON) * EPSILON), ((int)(y / EPSILON) * EPSILON)));
}

void Equation::results(double from_x, double to_x, double from_y, double to_y) {
    for (double x = from_x; x < to_x; x += step) {
        for (double y = from_y; y < to_y; y += step) {
            Newton_method(x, y);
        }
    }
    
    for (auto itr = roots.begin(); itr != roots.end(); itr++) {
        std::cout << "(" <<itr->first << ", " << itr->second << ")"<< std::endl; 
    }
}

int main(int argc, char* argv[]){ 
    Equation e; 
    e.results(-2, 2, -2, 2); 
}

