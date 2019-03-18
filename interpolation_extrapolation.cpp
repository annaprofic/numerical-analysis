#include <random>
#include <iostream>

const int size = 10; 
const int v = 4; 
const int h = 1; 

void gen_pairs(std::pair<int, double>* pairs){ 
    std::random_device random;
    std::mt19937 generator(random());
    std::uniform_real_distribution<double> dist(0, 1);
    for (int i = 0; i < size; i++){
        pairs[i] = std::pair<int, double>(i + 1, dist(generator));
    }
}

double* gen_u(std::pair<int, double>* pairs) {
    double *b = new double[size]; 
    double *u = new double[size];

    for (int i = 0; i < size; i++){ 
        b[i] = 1.f / h * (pairs[i + 1].second - pairs[i].second);
    }
    u[0] = 0; 
    for(int i = 1; i < size; i++) {
        u[i] = 6 * (b[i] - b[i - 1]);
    }
    return u;
}

double** gen_matrix() { 
    double** arr = new double*[size]; 
    for(int i = 0; i < size; i ++) { 
        arr[i] = new double[size]; 
    }
    for (int i = 0; i < size; i ++){ 
        for (int j = 0; j < size; j ++) {
            arr[i][j] = 0; 
        }
    }
    
    arr[0][0] = v; 
    arr[1][0] = h; 
    arr[size-1][size-1] = v; 
    arr[size-2][size-1] = h; 

    for (int i = 1; i < size - 1; i ++) { 
        arr[i][i] = v; 
        arr[i+1][i] = h; 
        arr[i-1][i] = h;  
    }
    return arr; 
} 

void gauss(double** arr, double* z) { 

    for(int i = 0; i < size - 2; i++) {
        double div = arr[i][i];
        for(int j = i; j < size; j++) {
            arr[i][j] /= div;
        }
        z[i + 1] /= div;
        
        if(i + 1 < size) {
            
            for(int it = i + 1; it < size - 2; it ++) {
                double div2 = arr[i][i] * arr[it][i]; 
                for(int j = i; j < size; j++) {
                    arr[it][j] -= arr[i][j] * div2;
                }
                
                z[it + 1] -= z[i + 1] * div2; 
            } 
        }
    }
}

void back_substitution(double** arr, double*z) { 
    int i; 
    for (int j = size - 3; j > 0; j --){ 
        for (i = j - 1 ; i >= 0; i --) { 
            z[i + 1] -= arr[i][j] * z[j + 1];
        }
    }    
} 

void splined(double* z, std::pair<int, double>* pairs){ 
    for (int i = 0; i < size - 1; i++){ 
        std::cout << "S" << i << "[x] = " << z[i+1]/(6*h) << "(x - " << pairs[i].first << ")^3 + " << z[i]/(6*h) << 
        "(" << pairs[i + 1].first << " - x)^3 + " << pairs[i + 1].second / h - z[i + 1] / 6 << "(x - " << pairs[i].first << 
        ") + " << pairs[i].second / h - z[i] / 6 << "(" << pairs[i + 1].first << " - x)" << std::endl; 
    }
}

void extrapolation(int left, int right, std::pair<int, double>* pairs, int func){ 
    
    std::cout << "L" << func << "[x] = " << std::flush;
    std::string add; 
    int n = 0; 
    for (int i = left; i < right; i++) {
        add.append(std::to_string(pairs[i].second) + "("); 
        for(int j = left; j < right; j++) {
            if(j != i) {
                n = pairs[i].first - pairs[j].first; 
                add.append("(x - " + std::to_string(pairs[j].first) + ") / " + std::to_string(n)); 
            } else continue; 
        }
        add.append(")"); 
        if ((right - 1) != i) { 
            add.append("+"); 
        }
    }
    std::cout << add << std::endl; 
}

int main(int argc, char* argv[]) { 
    std::pair<int, double>* pairs = new std::pair<int, double>[size]; 
    gen_pairs(pairs); 
    double** m = gen_matrix(); 
    double* z = gen_u(pairs); 
    gauss(m, z); 
    back_substitution(m, z); 
    splined(z, pairs); 

    std::cout << "points={" << std::endl;
    for (int i = 0; i < size; i++){
        std::cout << "{" << pairs[i].first << "," << pairs[i].second << "}";
        if (i != size - 1)
            std::cout << ",";
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;

    std::cout << "Show[";
    std::cout << "Plot[Piecewise[{";

    for (int i = 0; i < size - 1; i++) {
        std::cout << "{S" << i << "[x]," << pairs[i].first
                  << "<x<" << pairs[i + 1].first << "}";
        if (i != size - 2)
            std::cout << ",";
        std::cout << std::endl;
    }
    std::cout << "}], {x, 0, 12}]]" << std::endl;

    extrapolation(0, 3, pairs, 1); 
    extrapolation(6, 10, pairs, 2); 

}