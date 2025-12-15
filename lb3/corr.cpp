#include <iostream>
#include <vector>
#include <cmath>

double corr(std::vector<int>& a, std::vector<int>& b) {
    if (a.size() != b.size()) {
        printf("Векторы должны быть одинакового размера\n");
        return 0;
    }
    double result = 0;
    for (int i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    
    return result;
}

double norm_corr(std::vector<int>& a, std::vector<int>& b) {
    double unnormal_cor = corr(a, b);
    double norm_coeff_a = 0;
    double norm_coeff_b = 0;
    
    for (int i = 0; i < a.size(); ++i) {
        norm_coeff_a += std::pow(a[i], 2);
        norm_coeff_b += std::pow(b[i], 2);
    }
    
    return unnormal_cor / std::sqrt(norm_coeff_a * norm_coeff_b);
}

int main() {
    std::vector<int> a = {6, 2, 3, -2, -4, -4, 1, 1};
    std::vector<int> b = {3, 1, 5, 0, -3, -4, 2, 3};
    std::vector<int> c = {-1, -1, 3, -9, 2, -8, 4, -4};
    
    std::vector<std::vector<int>> all_vec = {a, b, c};
    std::vector<char> ls = {'a', 'b', 'c'};
    
    printf("   Ненормированная корреляция\n");
    printf("\t a \t b \t c\n");
    printf("-----------------------------\n");
    for (int i = 0; i < 3; ++i) {
        printf("%c\t", ls[i]);
        for (int j = 0; j < 3; ++j) {
            if (all_vec[i] != all_vec[j])
            {
            printf("%.0f\t", corr(all_vec[i], all_vec[j]));
            }
            else
            {
                printf("\t");
            }
        }
        printf("\n");
        printf("-----------------------------\n");
    }
    
    printf("\n   Нормированная корреляция\n");
    printf("\t a \t b \t c\n");
    printf("-----------------------------\n");
    for (int i = 0; i < 3; ++i) {
        printf("%c\t", ls[i]);
        for (int j = 0; j < 3; ++j) {
            printf("%.2f\t", norm_corr(all_vec[i], all_vec[j]));
        }
        printf("\n");
        printf("-----------------------------\n");
    }
    
    return 0;
}