#include <iostream> 
#include <cmath>


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

double **createQ(int size) { 
    double **Q = new double*[size]; 
    
    for(int i = 0; i < size; i ++) { 
        Q[i] = new double[size]; 
    }

    for (int i = 0; i < size; i ++){ 
        for (int j = 0; j < size; j ++) {
            Q[i][j] = 0; 
        }
    }
    return Q; 
}

double **createR(int size) { 
    double **R = new double*[size]; 
    
    for(int i = 0; i < size; i ++) { 
        R[i] = new double[size]; 
    }

    for (int i = 0; i < size; i ++){ 
        for (int j = 0; j < size; j ++) {
            R[i][j] = 0; 
        }
    }
    return R; 
}

void print_matrix(int size, double** tab) { 

    printf("Matrix %dx%d\n", size, size);  
    for(int i = 0; i < size; i++) {
          for(int j = 0; j < size; j++) {
               printf("%.5f  ", tab[i][j]);
          }
          printf("\n");
    }
}

double** transposition(int size, double** tab) {

    double** trans = new double*[size]; 

    for(int i = 0; i < size; i ++) { 
        trans[i] = new double[size]; 
    }

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            trans[j][i] = tab[i][j];
        }
    }
    return trans; 
}

double scalar(int size, double*a, double*b) { 

    double scalar = 0; 
    for( int i = 0; i < size;  i++) { 
        scalar += a[i] * b[i]; 
    }

    return scalar; 
}

double norma(int size, double*a){ 
    double norma = 0; 
    for (int i = 0; i < size; i++) { 
        norma += pow(a[i], 2);
    }
    return sqrt(norma); 
}


double*sub (int size, double*a, double *b) { 

    for ( int i = 0; i < size; i ++) { 
        a[i] -= b[i]; 
    }
    return a; 
}

double *e (int size, double *u) { 

    double *e = new double[size]; 
    double n = norma(size, u);
    for (int i = 0; i < size; i ++){
        u[i] = u[i]/n; 
    }
    return u; 
}

double *mult(int size, double *a, double mul) { 
    double *m = new double[size]; 
    for (int i = 0; i < size; i ++){
        m[i] = a[i] * mul; 
    }
    return m; 
}

void multiply(int size, double** tab, double **R, double **Q ){ 


    int i, j, k; 
    for (i = 0; i < size; i++) { 
        for (j = 0; j < size; j++) { 
            tab[i][j] = 0; 
            for (k = 0; k < size; k++) 
                tab[i][j] += R[i][k] * Q[k][j]; 
        } 
    }

}


void found(double** tab, double** Q, double** R, int size){ 

    double *u = new double[size]; 
    double *ev = new double[size]; 



    for ( int i = 0; i < size; i++ ) { 

        for (int k = 0; k < size; k++) { 
            u[k] = tab[i][k]; 
        }
 
        for( int j = 0; j < i; j ++) { 
            u = sub(size, u, mult(size, Q[j], scalar(size, u , Q[j])));
        }

        ev = e(size, u); 
     
        for(int j = 0; j < size; j++) {
            Q[i][j] = ev[j]; 
        }

        std::cout <<std::endl;
        for(int j = i; j < size; j++) { 
            R[i][j] = scalar(size, tab[j], Q[i]);
        }
        
    }  

    Q = transposition(size, Q);   
    multiply(size, tab, R, Q); 

}


int main(int args, char *argv[]) { 

    for(int c = 1; c < 2; c++) {
        for (int size = 4; size <= 10; size++) {
            double **arr = gen_matrix(size, c); 
            double **Q = createQ(size);
            double **R = createR(size);  
            for(int i = 1; i < 70; i++){ 
                found(arr, Q, R,  size); 
                printf("\nc = %d\n", c); 
                print_matrix(size, arr); 
            }  
            std::cout<<std::endl;
        }
    }
    return 0; 
}