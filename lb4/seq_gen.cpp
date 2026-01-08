#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>

const int j_num = 10;      
const int bit_count = 5;     
const int seq_length = 31;   

void display_sequence(const std::vector<int>& seq) {
    for (int val : seq) {
        std::cout << val << " ";
    }
}

std::vector<int> rotate_sequence(const std::vector<int>& seq, int offset) {
    std::vector<int> result = seq;
    offset %= seq.size();
    
    if (offset > 0) {
        std::rotate(result.rbegin(), result.rbegin() + offset, result.rend());
    }
    
    return result;
}

std::vector<int> generate_m_sequence(int seed, const std::vector<int>& feedback_taps) {
    std::vector<int> reg(bit_count);
    std::vector<int> output(seq_length);
    
    for (int i = 0; i < bit_count; ++i) {
        reg[i] = (seed >> i) & 1;
    }
    
    for (int step = 0; step < seq_length; ++step) {
        output[step] = reg[0];
        
        int feedback_bit = 0;
        for (int tap : feedback_taps) {
            feedback_bit ^= reg[tap];
        }
        
        for (int i = 0; i < bit_count - 1; ++i) {
            reg[i] = reg[i + 1];
        }
        reg[bit_count - 1] = feedback_bit;
    }
    
    return output;
}

int count_matches(const std::vector<int>& seq1, const std::vector<int>& seq2) {
    int match_count = 0;
    for (size_t i = 0; i < seq1.size(); ++i) {
        if (seq1[i] == seq2[i]) match_count++;
    }
    return match_count;
}

std::vector<int> generate_gold_sequence(const std::vector<int>& m_seq1, 
                                       const std::vector<int>& m_seq2) {
    std::vector<int> gold_seq(seq_length);
    
    for (int i = 0; i < seq_length; ++i) {
        gold_seq[i] = m_seq1[i] ^ m_seq2[i];
    }
    
    return gold_seq;
}

double compute_correlation(const std::vector<int>& seq1, const std::vector<int>& seq2) {
    int matches = count_matches(seq1, seq2);
    return (2.0 * matches - seq_length) / seq_length;
}

int main() {
    int base_seq1 = j_num;            
    int base_seq2 = j_num + 7;        
    int base_seq3 = j_num + 1;        
    int base_seq4 = j_num + 2;       
    
    std::vector<int> poly1 = {0, 2};     
    std::vector<int> poly2 = {0, 1, 3, 4};
    
    std::vector<int> m_seq_a = generate_m_sequence(base_seq1, poly1);
    std::vector<int> m_seq_b = generate_m_sequence(base_seq2, poly2);
    std::vector<int> m_seq_c = generate_m_sequence(base_seq3, poly1);
    std::vector<int> m_seq_d = generate_m_sequence(base_seq4, poly2);
    
    int best_shift = 0;
    double best_corr = -1.0;
    
    for (int shift = 0; shift < seq_length; ++shift) {
        std::vector<int> shifted_seq = rotate_sequence(m_seq_b, shift);
        std::vector<int> temp_gold = generate_gold_sequence(m_seq_a, shifted_seq);
        
        int ones_count = 0;
        for (int bit : temp_gold) ones_count += bit;
        
        if (abs(ones_count - seq_length/2) <= 1) {
            best_shift = shift;
            break;
        }
    }
    
    std::vector<int> shifted_m_seq_b = rotate_sequence(m_seq_b, best_shift);
    std::vector<int> gold_seq_1 = generate_gold_sequence(m_seq_a, shifted_m_seq_b);
    
    std::vector<int> shifted_m_seq_d = rotate_sequence(m_seq_d, best_shift);
    std::vector<int> gold_seq_2 = generate_gold_sequence(m_seq_c, shifted_m_seq_d);
    
    std::cout << "Исходные параметры:\n";
    std::cout << "Номер в журнале: " << j_num << "\n";
    std::cout << "Полином 1: x^5 + x^3 + 1\n";
    std::cout << "Полином 2: x^5 + x^4 + x^2 + x + 1\n\n";
    
    std::cout << "Последовательность Голда 1:\n";
    display_sequence(gold_seq_1);
    
    int ones_count = 0;
    for (int bit : gold_seq_1) ones_count += bit;
    std::cout << "\nБаланс: " << ones_count << " единиц, " 
              << (seq_length - ones_count) << " нулей\n";
    
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
              << (seq_length - ones_count) << " нулей\n";
    
    double cross_corr = compute_correlation(gold_seq_1, gold_seq_2);
    std::cout << "\nВзаимная корреляция между последовательностями: " 
              << cross_corr << "\n";
    
    std::cout << "Детальный анализ автокорреляции Seq1:\n";
    std::cout << "Сдвиг | Корреляция | Совпадения\n";
    
    for (int lag = 0; lag < 10; ++lag) {
        std::vector<int> shifted = rotate_sequence(gold_seq_1, lag);
        int matches = count_matches(gold_seq_1, shifted);
        double corr = compute_correlation(gold_seq_1, shifted);
        std::cout << lag << "     | " << corr << "     | " 
                  << matches << "/" << seq_length << "\n";
    }
    
    return 0;
}