#include <iostream> 
#include <cmath>
const int SIZE = 6; 
const double epsilon = 1e-15;
const double X = (double)1/(double)12; 
double data[][6] = {{19, 13, 10, 10, 13, -17}, 
                    {13, 13, 10, 10, -11, 13}, 
                    {10, 10, 10, -2, 10, 10}, 
                    {10, 10, -2, 10, 10, 10}, 
                    {13, -11, 10, 10, 13, 13}, 
                    {-17, 13, 10, 10, 13, 19}};


double** load_matrix(){ 
    double **arr = new double*[SIZE];    
    for(int i = 0; i < SIZE; i ++) { 
        arr[i] = new double[SIZE]; 
    }
    for (int i = 0; i < SIZE; i ++){ 
        for (int j = 0; j < SIZE; j ++) {
            arr[i][j] = data[i][j]; 
        }
    }
    return arr;  
}

double* load_vector() { 
    double* vector = new double[SIZE]; 
    for (int i = 0, x = 200; i < SIZE; i++) { 
        if (i < 3) { 
            vector[i] = x;
            x /= 2; 
        }
        if (i == 3) x = 200; 
        if ( i >= 3) { 
            vector[i] = x; 
            x /= 2; 
        }
    }
    return vector; 
}

double scalar(double*a, double*b) { 

    double scalar = 0; 
    for( int i = 0; i < SIZE;  i++) { 
        scalar += a[i] * b[i]; 
    }

    return scalar; 
}

void multiply_arr(double** arr) { 
    for (int i = 0; i < SIZE; i++){ 
        for (int j = 0; j < SIZE; j++) {
            arr[i][j] *= X; 
        }
    }
}
void multiply_v(double* v, double x) { 
    for (int i = 0; i < SIZE; i++){ 
        v[i] *= x; 
    }
}

void print_matrix(double** arr) { 
    for(int i = 0; i <  SIZE; i++) {
        for(int j = 0; j <  SIZE; j++) {
            printf("\x1B[3m%.5f  \x1B[0m", arr[i][j]);
        }
        printf("\n");
    }
}

void print_vector(double* v) { 
    std::cout << "[" <<std::flush;
    for(int i = 0; i <  SIZE; i++) {
        printf("\x1B[3m%.17f\x1B[0m", v[i]);
        if (i < SIZE - 1) { 
            std::cout << ", " << std::flush; 
        }
    }
    std::cout << "]" << std::endl;
}

double norm(double* vector){ 
    double norm = 0; 
    for (int i = 0; i < SIZE; i++) { 
        norm += pow(vector[i], 2);
    }
    return sqrt(norm); 
}

void normalize(double* vector){ 
    double norm_ = norm(vector); 
    for (int i = 0; i < SIZE; i++){
        vector[i] /= norm_; 
    }
}

double* multiply_by_v(double** arr, double* v) { 
    double* arr_by_v = new double[SIZE]; 
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++){
            arr_by_v[i] += arr[i][j] * v[j]; 
        }
    }
    return arr_by_v; 
}

double* sub_v(double* v1, double* v2){ 
    double *sub = new double[SIZE]; 
    for (int i = 0; i < SIZE; i++){ 
        sub[i] = v1[i] - v2[i]; 
    }
    return sub; 
}

double* copying(double* v){ 
    double* v_copy = new double[SIZE];
    for(int i = 0; i < SIZE; i++) {    
        v_copy[i] = v[i];        
    }
    return v_copy; 
}

double* perpendicular(double* v) { 
    double* p = new double[SIZE]; 
    double s = 0; 
    for (int i = 0; i < SIZE - 1; i++) { 
        p[i] = floor(rand() * 1000);
        s += p[i] * v[i]; 
    }
    p[SIZE - 1] = -s * ((double)1 / v[SIZE - 1]);
    return p;  
}


int main() {
    double** arr = load_matrix();
    double* vector = load_vector();  
    double* vector2 = new double[SIZE]; 
    double* vector2_norm = new double[SIZE]; 
    double* eig_v1 = new double[SIZE]; 
    double* eig1_c = new double[SIZE]; 
    double* eig_v2 = new double[SIZE]; 
    double* sub = new double[SIZE];
    double* sub_temp = new double[SIZE];
    double l1, l2; 
    
    multiply_arr(arr); 
    normalize(vector); 
    while (1) { 
        vector2_norm = multiply_by_v(arr, vector); 
        vector2 = copying(vector2_norm); 
        normalize(vector2_norm);     
        sub = sub_v(vector, vector2_norm);  
        if(epsilon > norm(sub)) { 
            eig_v1 = copying(vector2_norm); 
            l1 = norm(vector2); 
            break; 
        }
        vector = copying(vector2_norm); 
    }
    std::cout << "Dominant eigenvalue: " << l1 << std::endl; 
    std::cout << "Dominant eigenvector: " << std::flush; 
    print_vector(eig_v1); 
    vector = perpendicular(eig_v1); 

    while(1) { 
        eig1_c = copying(eig_v1); 
        vector2_norm = multiply_by_v(arr, vector);  
        multiply_v(eig1_c, scalar(eig_v1, vector2_norm)); 
        sub_temp = sub_v(vector2_norm, eig1_c); 
        vector2_norm = sub_temp; 
        vector2 = copying(vector2_norm); 
        normalize(vector2_norm); 
        sub = sub_v(vector, vector2_norm); 
        if(epsilon > norm(sub)) { 
            eig_v2 = copying(vector2_norm); 
            l2 = norm(vector2); 
            break; 
        }
        vector = copying(vector2_norm); 
    }

    std::cout << "Second dominant eigenvalue: " << l2 << std::endl; 
    std::cout << "Second dominant eigenvector: " << std::flush; 
    print_vector(eig_v2); 

}