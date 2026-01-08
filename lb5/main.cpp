#include <iostream>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <ctime>

template <typename T>
std::vector<T> rm_front(const std::vector<T>& vec){
    std::vector<T> result;
    for(int i = 1; i < vec.size(); ++i){
        result.push_back(vec[i]);
    }
    return result;
}

template <typename T>
std::vector<T> push_front(const std::vector<T>& vec, T value){
    std::vector<T> result(vec.size() + 1);
    for(int i = 1; i < vec.size()+1; ++i){
        result[i] = vec[i-1];
    }
    result[0] = value;
    return result;
}

template <typename T>
void print_vec(const std::vector<T>& vec){
    for(T el : vec){
        printf("%d", el);
    }
    printf("\n");
}

template <typename T>
std::vector<T> crc(std::vector<T> data, std::vector<T> poly) {
    for (int i = 0; i < poly.size() - 1; ++i)
        data.push_back(0);
    
    std::vector<T> result = data;
    
    while (result.size() >= poly.size()) {
        for (int i = 0; i < poly.size(); ++i)
            result[i] ^= poly[i];
        
        while (!result.empty() && result[0] == 0){
            result = rm_front(result);
        }
    }
    
    while (result.size() < poly.size() - 1){
        result = push_front(result, static_cast<T>(0));
    }
    
    return result;
}

template <typename T>
bool check_crc(const std::vector<T>& data_with_crc, const std::vector<T>& poly){
    std::vector<T> result = data_with_crc;
    
    while (result.size() >= poly.size()) {
        for (int i = 0; i < poly.size(); ++i)
            result[i] ^= poly[i];
        
        while (!result.empty() && result[0] == 0){
            result = rm_front(result);
        }
    }
    
    for(T bit : result){
        if(bit != 0){
            return false;
        }
    }
    
    return true;
}

template <typename T>
std::vector<T> generate_data(int N){
    std::vector<T> result(N);
    for(int i = 0; i < N; ++i){
        result[i] = rand() % 2;
    }
    return result;
}

template <typename T>
void crc_test(const std::vector<T>& data, const std::vector<T>& polynome, const char* test_name){
    printf("\n%s\n", test_name);
    for(int i = 0; i < strlen(test_name); ++i) printf("-");
    printf("\n");
    
    printf("Данные (%d бит): ", (int)data.size());
    print_vec(data);
    
    printf("Полином (%d бит):   ", (int)polynome.size());
    print_vec(polynome);
    
    std::vector<T> tx_crc = crc(data, polynome);
    printf("CRC на передатчике: ");
    print_vec(tx_crc);
    
    std::vector<T> data_with_crc;
    data_with_crc.reserve(data.size() + tx_crc.size());
    data_with_crc.insert(data_with_crc.begin(), data.begin(), data.end());
    data_with_crc.insert(data_with_crc.end(), tx_crc.begin(), tx_crc.end());    
    
    printf("Данные + CRC:       ");
    print_vec(data_with_crc);
    
    bool check_result = check_crc(data_with_crc, polynome);
    printf("Проверка на приемнике: %s\n", check_result ? "Ошибок нет" : "Обнаружена ошибка");
    
    if(data_with_crc.size() > 0){
        std::vector<T> corrupted = data_with_crc;
        corrupted[0] = corrupted[0] ^ 1;
        bool error_check = check_crc(corrupted, polynome);
        printf("Тест с ошибкой:       %s\n", error_check ? "Ошибок нет" : "Обнаружена ошибка");
    }
}

template <typename T>
void test_single_errors(const std::vector<T>& data, const std::vector<T>& polynome){
    
    std::vector<T> tx_crc = crc(data, polynome);
    
    std::vector<T> full_data = data;
    for(T bit : tx_crc){
        full_data.push_back(bit);
    }
    
    int total_bits = full_data.size();
    int errors_detected = 0;
    int errors_undetected = 0;
    
    printf("Всего бит для тестирования: %d\n", total_bits);
    
    for(int i = 0; i < total_bits; ++i){
        std::vector<T> corrupted = full_data;
        corrupted[i] = corrupted[i] ^ 1;
        bool error_detected = !check_crc(corrupted, polynome);
        
        if(error_detected){
            errors_detected++;
        } else {
            errors_undetected++;
        }
    }
    
    printf("Всего протестировано бит: %d\n", total_bits);
    printf("Ошибок обнаружено:        %d\n", errors_detected);
    printf("Ошибок НЕ обнаружено:     %d\n", errors_undetected);
    printf("Эффективность обнаружения: %.2f%%\n", 
           (float)errors_detected / total_bits * 100);
    
    if(errors_undetected == 0){
        printf("\nВЫВОД: Все одиночные ошибки успешно обнаружены!\n");
    } else {
        printf("\nВНИМАНИЕ: Найдены необнаруженные ошибки!\n");
    }
}

int main(){
    srand(time(0));
    
    int num_in_journal = 10;
    
    std::vector<int8_t> polynome_v10 = {1, 1, 0, 1, 1, 1, 1, 0};
    
    printf("G = x^7 + x^6 + x^4 + x^3 + x^2 + x\n");
    printf("Бинарное представление: ");
    print_vec(polynome_v10);
    printf("Длина CRC: %d бит\n", (int)polynome_v10.size() - 1);

    
    int N1 = 20 + num_in_journal;
    std::vector<int8_t> data1 = generate_data<int8_t>(N1);
    crc_test(data1, polynome_v10, "ТЕСТ 1: N = 30 бит (20 + номер в журнале)");
    
    int N2 = 250;
    std::vector<int8_t> data2 = generate_data<int8_t>(N2);
    crc_test(data2, polynome_v10, "ТЕСТ 2: N = 250 бит");
    
    test_single_errors(data2, polynome_v10);

    return 0;
}