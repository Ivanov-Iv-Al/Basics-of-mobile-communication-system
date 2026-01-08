j_num = 10;         
bit_count = 5;          
seq_len = 2^bit_count - 1;  


base_seq_x = de2bi(j_num, bit_count, 'left-msb');
base_seq_y = de2bi(j_num + 7, bit_count, 'left-msb');
base_seq_x2 = de2bi(j_num + 1, bit_count, 'left-msb');
base_seq_y2 = de2bi(j_num + 2, bit_count, 'left-msb');

fprintf('Параметры генерации\n');
fprintf('Номер в журнале: %d\n', j_num);
fprintf('Длина последовательности: %d\n', seq_len);
fprintf('\nБазовые последовательности (5 бит):\n');
fprintf('x:  %s\n', num2str(base_seq_x));
fprintf('y:  %s\n', num2str(base_seq_y));
fprintf('x'': %s\n', num2str(base_seq_x2));
fprintf('y'': %s\n', num2str(base_seq_y2));

poly_set1 = [5, 3];        
poly_set2 = [5, 4, 2, 1];  

m_seq_x = create_m_sequence(base_seq_x, poly_set1);
m_seq_y = create_m_sequence(base_seq_y, poly_set2);
m_seq_x2 = create_m_sequence(base_seq_x2, poly_set1);
m_seq_y2 = create_m_sequence(base_seq_y2, poly_set2);

time_axis = 0:seq_len-1;

figure('Position', [100, 100, 800, 600]);
subplot(2,2,1);
stem(time_axis, m_seq_x, 'filled');
title('М-последовательность X');
xlabel('Отсчеты'); ylabel('Амплитуда');
grid on;

subplot(2,2,2);
stem(time_axis, m_seq_y, 'filled');
title('М-последовательность Y');
xlabel('Отсчеты'); ylabel('Амплитуда');
grid on;

subplot(2,2,3);
stem(time_axis, m_seq_x2, 'filled');
title('М-последовательность X''');
xlabel('Отсчеты'); ylabel('Амплитуда');
grid on;

subplot(2,2,4);
stem(time_axis, m_seq_y2, 'filled');
title('М-последовательность Y''');
xlabel('Отсчеты'); ylabel('Амплитуда');
grid on;

shift_val = find_best_shift(m_seq_x, m_seq_y);
gold_seq1 = xor(m_seq_x, circshift(m_seq_y, shift_val));
gold_seq2 = xor(m_seq_x2, circshift(m_seq_y2, shift_val));

fprintf('\nПоследовательности Голда\n');
fprintf('Оптимальный сдвиг: %d\n', shift_val);
fprintf('\nGold Seq1 (первые 20 бит):\n');
fprintf('%d ', gold_seq1(1:20)');
fprintf('\n\nGold Seq2 (первые 20 бит):\n');
fprintf('%d ', gold_seq2(1:20)');

analyze_balance(gold_seq1, 'Последовательность Голда 1');
analyze_balance(gold_seq2, 'Последовательность Голда 2');

fprintf('\nПроверка PN-свойств:\n');
fprintf('\nПоследовательность Голда 1:\n');
check_pn_properties(gold_seq1);
fprintf('\nПоследовательность Голда 2:\n');
check_pn_properties(gold_seq2);
fprintf('\nСлучайная последовательность:\n');
random_seq = randi([0 1], 1, seq_len);
check_pn_properties(random_seq);

figure('Position', [100, 100, 800, 400]);

subplot(1,2,1);
stem(time_axis, gold_seq1, 'filled');
title('Последовательность Голда 1');
xlabel('Отсчеты'); ylabel('Амплитуда');
xlim([0, 50]); grid on;

subplot(1,2,2);
stem(time_axis, gold_seq2, 'filled');
title('Последовательность Голда 2');
xlabel('Отсчеты'); ylabel('Амплитуда');
xlim([0, 50]); grid on;

auto_corr_gold1 = zeros(seq_len, 1);
auto_corr_gold2 = zeros(seq_len, 1);

for shift = 0:seq_len-1
    shifted1 = circshift(gold_seq1, shift);
    shifted2 = circshift(gold_seq2, shift);
    
    auto_corr_gold1(shift+1) = calc_normalized_corr(gold_seq1, shifted1);
    auto_corr_gold2(shift+1) = calc_normalized_corr(gold_seq2, shifted2);
end

cross_corr = zeros(seq_len, 1);
for shift = 0:seq_len-1
    shifted2 = circshift(gold_seq2, shift);
    cross_corr(shift+1) = calc_normalized_corr(gold_seq1, shifted2);
end

lag_range = -(seq_len-1):(seq_len-1);

figure('Position', [100, 100, 1000, 800]);

subplot(3,2,1);
stem(0:seq_len-1, auto_corr_gold1, 'filled');
title('Автокорреляция Gold Seq1');
xlabel('Сдвиг'); ylabel('Корреляция');
grid on; ylim([-1, 1]);

subplot(3,2,2);
stem(0:seq_len-1, auto_corr_gold2, 'filled');
title('Автокорреляция Gold Seq2');
xlabel('Сдвиг'); ylabel('Корреляция');
grid on; ylim([-1, 1]);

subplot(3,2,[3,4]);
stem(-(seq_len-1):(seq_len-1), [flipud(cross_corr(2:end)); cross_corr], 'filled');
title('Взаимная корреляция Gold Seq1 и Gold Seq2');
xlabel('Сдвиг'); ylabel('Корреляция');
grid on; ylim([-1, 1]);

subplot(3,2,5);
plot(0:seq_len-1, auto_corr_gold1, 'b-', 'LineWidth', 1.5);
hold on;
plot(0:seq_len-1, auto_corr_gold2, 'r--', 'LineWidth', 1.5);
title('Сравнение автокорреляционных функций');
xlabel('Сдвиг'); ylabel('Корреляция');
legend('Gold Seq1', 'Gold Seq2');
grid on; ylim([-1, 1]);

subplot(3,2,6);
histogram(auto_corr_gold1, 20, 'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(auto_corr_gold2, 20, 'FaceColor', 'r', 'FaceAlpha', 0.5);
title('Распределение значений корреляции');
xlabel('Корреляция'); ylabel('Частота');
legend('Gold Seq1', 'Gold Seq2');
grid on;

fprintf('\nАвтокорреляция Gold Seq1 (сдвиг 0): %.4f\n', auto_corr_gold1(1));
fprintf('Автокорреляция Gold Seq2 (сдвиг 0): %.4f\n', auto_corr_gold2(1));
fprintf('Взаимная корреляция (сдвиг 0): %.4f\n', cross_corr(1));
min_auto1 = min(auto_corr_gold1(2:end));
min_auto2 = min(auto_corr_gold2(2:end));
max_cross = max(abs(cross_corr(2:end)));

fprintf('\nМинимальная автокорреляция (ненулевые сдвиги):\n');
fprintf('  Gold Seq1: %.4f\n', min_auto1);
fprintf('  Gold Seq2: %.4f\n', min_auto2);
fprintf('Максимальная взаимная корреляция: %.4f\n', max_cross);

analyze_cycles(gold_seq1, 'Gold Seq1');
analyze_cycles(gold_seq2, 'Gold Seq2');

function m_seq = create_m_sequence(init_state, poly_taps)
    reg_len = length(init_state);
    m_seq = zeros(2^reg_len - 1, 1);
    reg = init_state;
    
    for idx = 1:length(m_seq)
        m_seq(idx) = reg(end);
        
        fb_bit = 0;
        for tap = poly_taps
            fb_bit = xor(fb_bit, reg(tap));
        end
        
        reg = circshift(reg, 1);
        reg(1) = fb_bit;
    end
end

function optimal_shift = find_best_shift(seq1, seq2)
    sequence_length = length(seq1);  
    optimal_shift = 0;
    best_balance = sequence_length;  
    
    for s = 0:sequence_length-1
        temp_gold = xor(seq1, circshift(seq2, s));
        ones_count = sum(temp_gold);
        
        balance_diff = abs(2*ones_count - sequence_length);
        if balance_diff < best_balance
            best_balance = balance_diff;
            optimal_shift = s;
        end
    end
end

function analyze_balance(sequence, name)
    ones_cnt = sum(sequence);
    zeros_cnt = length(sequence) - ones_cnt;
    
    fprintf('\n%s:\n', name);
    fprintf('  Единиц: %d (%.1f%%)\n', ones_cnt, 100*ones_cnt/length(sequence));
    fprintf('  Нулей:  %d (%.1f%%)\n', zeros_cnt, 100*zeros_cnt/length(sequence));
    fprintf('  Баланс: %+d\n', ones_cnt - zeros_cnt);
end

function corr_result = calc_normalized_corr(sig1, sig2)
    if length(sig1) ~= length(sig2)
        error('Сигналы должны иметь одинаковую длину');
    end
    
    norm_factor = sqrt(sum(sig1.^2) * sum(sig2.^2));
    if norm_factor == 0
        corr_result = 0;
    else
        corr_result = sum(sig1 .* sig2) / norm_factor;
    end
end

function analyze_cycles(sequence, name)
    fprintf('\n%s - анализ циклов:\n', name);
    
    current_bit = sequence(1);
    cycle_len = 1;
    cycle_counts = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    
    for i = 2:length(sequence)
        if sequence(i) == current_bit
            cycle_len = cycle_len + 1;
        else
            if ~isKey(cycle_counts, cycle_len)
                cycle_counts(cycle_len) = 0;
            end
            cycle_counts(cycle_len) = cycle_counts(cycle_len) + 1;
            
            current_bit = sequence(i);
            cycle_len = 1;
        end
    end
    
    if ~isKey(cycle_counts, cycle_len)
        cycle_counts(cycle_len) = 0;
    end
    cycle_counts(cycle_len) = cycle_counts(cycle_len) + 1;
    
    keys_arr = cell2mat(keys(cycle_counts));
    [sorted_keys, idx] = sort(keys_arr);
    
    for k = 1:length(sorted_keys)
        key_val = sorted_keys(k);
        fprintf('  Длина %d: %d циклов\n', key_val, cycle_counts(key_val));
    end
end

function check_pn_properties(seq)
    N = length(seq);
    
    ones_count = sum(seq);
    zeros_count = N - ones_count;
    balance_diff = abs(ones_count - zeros_count);
    
    runs = [];
    current_run = 1;
    for i = 2:N
        if seq(i) == seq(i-1)
            current_run = current_run + 1;
        else
            runs = [runs, current_run];
            current_run = 1;
        end
    end
    runs = [runs, current_run];
    
    max_run = max(runs);
    
    seq_bipolar = 2*seq - 1;
    autocorr_sidelobes = zeros(1, N-1);
    for tau = 1:N-1
        shifted_seq = circshift(seq_bipolar, tau);
        autocorr_sidelobes(tau) = abs(sum(seq_bipolar .* shifted_seq) / N);
    end
    max_sidelobe = max(autocorr_sidelobes);
    
    balance_ok = balance_diff <= 1;
    run_ok = max_run <= 5;
    autocorr_ok = max_sidelobe <= 0.3;
    
    fprintf('  Баланс: %d единиц, %d нулей (разница: %d)\n', ones_count, zeros_count, balance_diff);
    fprintf('  Макс.длина серии: %d\n', max_run);
    fprintf('  Макс.бок.лепесток: %.4f\n', max_sidelobe);
    fprintf('  Критерии: баланс(%d), серии(%d), автокорр(%d)\n', balance_ok, run_ok, autocorr_ok);
    
    is_pn = balance_ok && run_ok && autocorr_ok;
    fprintf('  PN-статус: %s\n', string(is_pn));
end
