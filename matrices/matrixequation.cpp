#include <iostream> 

double* gen_res(int size) {
     double *res = new double[size];

     for(int i = 0; i < size; i++) {
          res[i] = i+1;
     }
     return res;
}

double** gen_matrix(int size, int c) { 
    double** tab = new double*[size]; 

    for(int i = 0; i < size; i ++) { 
        tab[i] = new double[size]; 
    }

    for (int i = 0; i < size; i ++){ 
        for (int j = 0; j < size; j ++) {
            tab[i][j] = 0; 
        }
    }
    
    tab[0][0] = 6; 
    tab[1][0] = 1; 
    tab[size-1][size-1] = 6; 
    tab[size-2][size-1] = 1; 
    tab[0][size-1] = c; 
    tab[size-1][0] = c; 

    for (int i = 1; i < size - 1; i ++) { 
        tab[i][i] = 6; 
        tab[i+1][i] = 1; 
        tab[i-1][i] = 1;  
    }


    return tab; 
} 

void gauss(int size,  double** tab, double* res) { 

    for(int i = 0; i < size; i++) {
        double div = tab[i][i];
        for(int j = i; j < size; j++) {
            tab[i][j] /= div;
        }
        res[i] /= div;
        
        if(i+1 < size) {
            
            for(int it = i + 1; it < size; it ++) {
                double div2 = tab[i][i]*tab[it][i]; 
                for(int j = i; j < size; j++) {
                    tab[it][j] -= tab[i][j] * div2;
                }
                res[it] -= res[i] * div2; 
            } 
        }
    }

}

void print_matrix(int size, double** tab, double*res) { 

  for(int i = 0; i < size; i++) {
          for(int j = 0; j < size; j++) {
               printf("%.5f ", tab[i][j]);
          }
          printf("|| %.5f", res[i]);
          printf("\n");
     }
}

void back_substitution(int size,  double** tab, double*res) { 

    int i; 
    for (int j = size - 1; j > 0; j --){ 
        for (i = j -1 ; i >= 0; i --) { 
            res[i] -= tab[i][j] * res[j];
        }
    }
     
} 

void print_results(int size, double*res) { 
    int i = 0; 
    int j = 0; 
    std::cout<<std::endl;
    for (i = 0; i < size; i ++) {
        printf("x%d = %.8f  ", i+1, res[i]); 
    }
}

int main(int args, char *argv[]) { 

    for(int c = 0; c < 2; c++) {

        for (int size = 4; size <= 10; size++) {

            std::cout<<std::endl;
            double **matrix = gen_matrix(size, c);
            double *res = gen_res(size); 

            printf("\nInput matrix, size %dx%d for c = %d\n", size, size, c);
            print_matrix(size, matrix, res);  

            std::cout<< "\nMatrix after Gaussian elimination: " <<std::endl;
            gauss(size, matrix, res); 
            print_matrix(size, matrix, res);  

            std::cout<< "\nEquation results: " <<std::flush;
            back_substitution(size, matrix, res);
            print_results(size, res); 
        }
    }
    return 0; 
} 