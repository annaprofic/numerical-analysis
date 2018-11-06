#include <iostream> 
#include <cmath>

void print_matrix(int size, double** arr) { 
    
    for(int i = 0; i < size; i++) {

        for(int j = 0; j < size; j++) {
            printf("\x1B[3m%.5f  \x1B[0m", arr[i][j]);
        }
        printf("\n");
    }
}

double** load_matrix(int size){ 
    double **arr = new double*[size]; 
    
    for(int i = 0; i < size; i ++) { 
        arr[i] = new double[size]; 
    }

    for (int i = 0; i < size; i ++){ 
        for (int j = 0; j < size; j ++) {
            arr[i][j] = 0; 
        }
    }
    std::cout << "\nPlease, enter the " << size << "x" << size << " matrix: " << std::endl;
    for (int i = 0; i < size; i ++ ){ 
        for( int j = 0; j < size; j ++) { 
            std::cin >> arr[i][j];
        }
        
    }

    return arr; 
}

double** copying(int size, double** arr){ 
    double** arr_copy = new double*[size];

    for(int i = 0; i < size; i ++) { 
        arr_copy[i] = new double[size]; 
    }

    for(int i = 0; i < size; i++) { 
        for (int j = 0; j < size; j++) { 
            arr_copy[i][j] = arr[i][j]; 
        }
    }
    return arr_copy; 
}

double **createL(int size, double** arr) { 
    double **L = new double*[size]; 
    
    for(int i = 0; i < size; i ++) { 
        L[i] = new double[size]; 
    }
    int i, j; 
    for (i = 0; i < size; i ++){ 
        for ( j = 0; j < size; j ++) {

            L[i][j] = 0; 
            L[i][i] = 1; 
            if (j < i) { 
                L[i][j] = arr[i][j]; 
            }
           
        }
    }

    return L; 
}

double **createU(int size, double** arr) { 
    double **U = new double*[size]; 
    
    for(int i = 0; i < size; i ++) { 
        U[i] = new double[size]; 
    }

    for (int i = 0; i < size; i ++){ 
        for (int j = 0; j < size; j ++) {
            U[i][j] = 0; 
            if (i <= j) { 
                U[i][j] = arr[i][j]; 
            }
        }
    }

    return U; 
}


void permute(int size, double** arr, int a, int b) { 
    double* temp = new double[size]; 
    for (int i = 0; i < size; i ++) { 
        temp[i] = arr[a][i]; 
        arr[a][i] = arr[b][i]; 
        arr[b][i] = temp[i]; 
    }   
}

double** permutation_matrix(int size, double *p) { 
    double** p_ = new double *[size]; 

      for(int i = 0; i < size; i ++) { 
        p_[i] = new double[size]; 
    }

    for (int i = 0; i < size; i ++){ 
        for (int j = 0; j < size; j ++) {
            p_[i][j] = 0; 
            p_[i][(int)p[i]] = 1; 
        }
    }
    return p_; 
}


void gauss(int size,  double** arr) { 

    int i = 0, j = 0, max = 0; 

    double *p = new double[size]; 
    for ( int i = 0; i < size; i ++) { 
        p[i] = i; 
    }

    for(i = 0; i < size; i++) { 
        max = i; 
        for (j = i + 1; j < size; j ++) { 
            if (abs(arr[max][i]) < abs(arr[j][i])) { 
                max = j; 
            }
        }
        permute(size, arr, i, max); 
        std::swap(p[i], p[max]); 

        if(i+1 < size) {
            for(int it = i + 1; it < size; it ++) {
                arr[it][i] /= arr[i][i]; 
                double div = arr[it][i]; 
                for(int j = i + 1; j < size; j++) {
                    arr[it][j] -= arr[i][j]*div;
                }
            } 
        }
    }
    double ** p_ = permutation_matrix(size,p); 
    std::cout << "\n\x1B[1mP\x1B[0m matrix: " <<std::endl;
    print_matrix(size, p_); 

}

void multiply(int size, double **L, double **U ){ 
    double **mult = new double*[size]; 
    
    for(int i = 0; i < size; i ++) { 
        mult[i] = new double[size]; 
    }

    int i, j, k; 
    for (i = 0; i < size; i++) { 
        for (j = 0; j < size; j++) { 
            mult[i][j] = 0; 
            for (k = 0; k < size; k++) 
                mult[i][j] += L[i][k] * U[k][j]; 
        } 
    }
    std::cout << "\x1B[1mL * U\x1B[0m matrix: " <<std::endl;
    print_matrix(size, mult); 
}



int main (int args, char *argv[]) { 

    //ustawienie rozmiaru macierzy 
    int size; 
    std::cout << "\nPlease, enter the \x1B[1msize\x1B[0m of your matrix : " << std::flush;
    std::cin >> size; 

    //laduje macierz uzytownika
    double** arr = load_matrix(size); 
    double** arr_copy = copying(size, arr);

    std::cout << "\x1B[3m\nFor any n × n nonsingular matrix A, there exists a permutation matrix P such that" 
              <<  "\nPA has an LU factorization.\nThat is, such that \n\x1B[1mPA = LU.\x1B[0m\n" << std::endl;

    //eliminacja Gaussa 
    gauss(size, arr); 
    //z macierzy wyjsciowej po eliminacji Gaussa tworzymy macierze L i U 
    double ** L = createL(size, arr); 
    double ** U = createU(size, arr); 

    std::cout << "\n\x1B[1mA\x1B[0m matrix: " << std::endl;
    print_matrix(size, arr_copy);
    std::cout << "\n\x1B[1mL\x1B[0m matrix: " << std::endl;
    print_matrix(size, L);
    std::cout << "\n\x1B[1mU\x1B[0m matrix: " << std::endl;
    print_matrix(size, U);

    std::cout << "\x1B[3m\nNow we can multiply \x1B[1mL\x1B[0m by \x1B[1mU\x1B[0m.\x1B[0m\n" << std::endl;
    multiply(size, L, U); 
    std::cout << "\x1B[3m\n\x1B[1mL * U\x1B[0m\x1B[3m matrix should be equal to \x1B[1mP * A.\x1B[0m\n" << std::endl;

} 
