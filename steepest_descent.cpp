#include <iostream>
#include <cmath>
#include <random>

class Extreme {

    public:
        Extreme();
        int size; 
        double* result;
        double* steepest_descent(double x_from, double y_from); 
        double generate_function(double x, double y);  
    private:
        double* data;  
        void generate_data(); 
        
};

Extreme::Extreme(): size(1000) {  
    data = new double[size];
    result = new double[2];  
    generate_data(); 
}

void Extreme::generate_data(){
    int i = 0; 
    for (int i = 0; i < 250; i++) {
        data[i] = ((double) 1 / (double) (i + 1));
        data[250 + i] = (double) i;
    }
    for (int i = 0; i < 500; i++){
        data[500 + i] = -data[i];
    } 
}

double Extreme::generate_function(double x, double y){ 
    return pow((2-x), 2) + 46 * pow(y - pow(x, 2), 2); 
}

double* Extreme::steepest_descent(double x_from, double y_from){
    double curr_x = x_from;
    double curr_y = y_from;
    while (true) {
        double current_func_value = generate_function(curr_x, curr_y);
        double tempX = curr_x, tempY = curr_y;
        double func_value = current_func_value;
        double tX, tY; 
        for (int i = 0; i < size; i++) {
            tX = data[i]; 
            for (int j = 0; j < size; j++) {
                tY = data[j]; 
                double nX = curr_x + tX;
                double nY = curr_y + tY;
                if (generate_function(nX, nY) < func_value){
                    func_value = generate_function(nX, nY);
                    tempX = nX;
                    tempY = nY;
                }
            }
        }
        if (current_func_value == func_value) {
            result[0] = curr_x; 
            result[1] = curr_y;
            return result;
        }
        curr_x = tempX;
        curr_y = tempY;
    }
}

int main(int argc, char*argv[]){ 
    srand(time(NULL)); 
    int i = 0; 
    std::random_device random;
    std::mt19937 generator(random());
    std::uniform_real_distribution<double> dist(0, 1);

    Extreme e; 
    double* zeroP = e.steepest_descent(0, 0);
    double minX = zeroP[0];
    double minY = zeroP[1];
    double minVal = e.generate_function(minX, minY);

    while (i < 10) {
        double randX = pow(dist(generator) * 10, dist(generator) * 5);
        double randY = pow(dist(generator) * 10, dist(generator) * 5);
        double* res = e.steepest_descent(randX, randY);
        if (e.generate_function(res[0], res[1]) < minVal) {
            minX = res[0];
            minY = res[1];
            minVal = e.generate_function(res[0], res[1]);
        }
        i++; 
    }

    printf("x_min = %f\n",  minX);
    printf("y_min = %f\n", minY);
    printf("F(x_min,y_min) = %f\n", minVal);
}

