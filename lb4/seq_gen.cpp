#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>

const int INIT_NUM = 10;      // Порядковый номер в журнале
const int BIT_WIDTH = 5;     // Разрядность чисел
const int SEQ_LENGTH = 31;   // 2^5 - 1

// Вспомогательная функция для отображения последовательности
void display_sequence(const std::vector<int>& seq) {
    for (int val : seq) {
        std::cout << val << " ";
    }
}

// Циклический сдвиг последовательности вправо
std::vector<int> rotate_sequence(const std::vector<int>& seq, int offset) {
    std::vector<int> result = seq;
    offset %= seq.size();
    
    if (offset > 0) {
        std::rotate(result.rbegin(), result.rbegin() + offset, result.rend());
    }
    
    return result;
}

// Генератор М-последовательности на основе полинома
std::vector<int> generate_m_sequence(int seed, const std::vector<int>& feedback_taps) {
    std::vector<int> reg(BIT_WIDTH);
    std::vector<int> output(SEQ_LENGTH);
    
    // Инициализация регистра
    for (int i = 0; i < BIT_WIDTH; ++i) {
        reg[i] = (seed >> i) & 1;
    }
    
    for (int step = 0; step < SEQ_LENGTH; ++step) {
        output[step] = reg[0];
        
        int feedback_bit = 0;
        for (int tap : feedback_taps) {
            feedback_bit ^= reg[tap];
        }
        
        // Сдвиг регистра
        for (int i = 0; i < BIT_WIDTH - 1; ++i) {
            reg[i] = reg[i + 1];
        }
        reg[BIT_WIDTH - 1] = feedback_bit;
    }
    
    return output;
}

// Подсчет совпадающих битов
int count_matches(const std::vector<int>& seq1, const std::vector<int>& seq2) {
    int match_count = 0;
    for (size_t i = 0; i < seq1.size(); ++i) {
        if (seq1[i] == seq2[i]) match_count++;
    }
    return match_count;
}

// Генерация последовательности Голда
std::vector<int> generate_gold_sequence(const std::vector<int>& m_seq1, 
                                       const std::vector<int>& m_seq2) {
    std::vector<int> gold_seq(SEQ_LENGTH);
    
    for (int i = 0; i < SEQ_LENGTH; ++i) {
        gold_seq[i] = m_seq1[i] ^ m_seq2[i];
    }
    
    return gold_seq;
}

// Вычисление корреляции
double compute_correlation(const std::vector<int>& seq1, const std::vector<int>& seq2) {
    int matches = count_matches(seq1, seq2);
    return (2.0 * matches - SEQ_LENGTH) / SEQ_LENGTH;
}

int main() {
    // Исходные данные
    int base_seq1 = INIT_NUM;            // x = 6 (00110)
    int base_seq2 = INIT_NUM + 7;        // y = 13 (01101)
    int base_seq3 = INIT_NUM + 1;        // x' = 7 (00111)
    int base_seq4 = INIT_NUM + 2;        // y' = 8 (01000)
    
    // Полиномы обратной связи (индексы от 0)
    std::vector<int> poly1 = {0, 2};      // x^5 + x^3 + 1
    std::vector<int> poly2 = {0, 1, 3, 4}; // x^5 + x^4 + x^2 + x + 1
    
    // Генерация М-последовательностей
    std::vector<int> m_seq_a = generate_m_sequence(base_seq1, poly1);
    std::vector<int> m_seq_b = generate_m_sequence(base_seq2, poly2);
    std::vector<int> m_seq_c = generate_m_sequence(base_seq3, poly1);
    std::vector<int> m_seq_d = generate_m_sequence(base_seq4, poly2);
    
    // Поиск оптимального сдвига
    int best_shift = 0;
    double best_corr = -1.0;
    
    for (int shift = 0; shift < SEQ_LENGTH; ++shift) {
        std::vector<int> shifted_seq = rotate_sequence(m_seq_b, shift);
        std::vector<int> temp_gold = generate_gold_sequence(m_seq_a, shifted_seq);
        
        int ones_count = 0;
        for (int bit : temp_gold) ones_count += bit;
        
        // Ищем последовательность, близкую к сбалансированной
        if (abs(ones_count - SEQ_LENGTH/2) <= 1) {
            best_shift = shift;
            break;
        }
    }
    
    // Генерация последовательностей Голда
    std::vector<int> shifted_m_seq_b = rotate_sequence(m_seq_b, best_shift);
    std::vector<int> gold_seq_1 = generate_gold_sequence(m_seq_a, shifted_m_seq_b);
    
    std::vector<int> shifted_m_seq_d = rotate_sequence(m_seq_d, best_shift);
    std::vector<int> gold_seq_2 = generate_gold_sequence(m_seq_c, shifted_m_seq_d);
    
    // Вывод результатов
    std::cout << "=======================================\n";
    std::cout << "Лабораторная работа №4: Генерация PN-последовательностей\n";
    std::cout << "=======================================\n\n";
    
    std::cout << "Исходные параметры:\n";
    std::cout << "Номер в журнале: " << INIT_NUM << "\n";
    std::cout << "Полином 1: x^5 + x^3 + 1\n";
    std::cout << "Полином 2: x^5 + x^4 + x^2 + x + 1\n\n";
    
    std::cout << "Последовательность Голда 1:\n";
    display_sequence(gold_seq_1);
    
    // Проверка баланса
    int ones_count = 0;
    for (int bit : gold_seq_1) ones_count += bit;
    std::cout << "\nБаланс: " << ones_count << " единиц, " 
              << (SEQ_LENGTH - ones_count) << " нулей\n";
    
    // Автокорреляционный анализ
    std::cout << "\nАнализ автокорреляции (первые 5 сдвигов):\n";
    for (int lag = 0; lag <= 5; ++lag) {
        std::vector<int> shifted = rotate_sequence(gold_seq_1, lag);
        double corr = compute_correlation(gold_seq_1, shifted);
        std::cout << "Сдвиг " << lag << ": " << corr << "\n";
    }
    
    std::cout << "\nПоследовательность Голда 2:\n";
    display_sequence(gold_seq_2);
    
    ones_count = 0;
    for (int bit : gold_seq_2) ones_count += bit;
    std::cout << "\nБаланс: " << ones_count << " единиц, " 
              << (SEQ_LENGTH - ones_count) << " нулей\n";
    
    // Взаимная корреляция
    double cross_corr = compute_correlation(gold_seq_1, gold_seq_2);
    std::cout << "\nВзаимная корреляция между последовательностями: " 
              << cross_corr << "\n";
    
    // Детальный анализ автокорреляции
    std::cout << "\n=======================================\n";
    std::cout << "Детальный анализ автокорреляции Seq1:\n";
    std::cout << "Сдвиг | Корреляция | Совпадения\n";
    std::cout << "---------------------------------------\n";
    
    for (int lag = 0; lag < 10; ++lag) {
        std::vector<int> shifted = rotate_sequence(gold_seq_1, lag);
        int matches = count_matches(gold_seq_1, shifted);
        double corr = compute_correlation(gold_seq_1, shifted);
        std::cout << lag << "     | " << corr << "     | " 
                  << matches << "/" << SEQ_LENGTH << "\n";
    }
    
    return 0;
}