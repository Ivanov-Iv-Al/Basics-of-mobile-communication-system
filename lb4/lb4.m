%% Лабораторная работа №4: Анализ PN-последовательностей
clear; clc; close all;

%% Параметры системы
journal_num = 6;         % Номер в журнале
bit_count = 5;          % Разрядность
seq_len = 2^bit_count - 1;  % Длина последовательности

%% Базовые последовательности
base_seq_x = de2bi(journal_num, bit_count, 'left-msb');
base_seq_y = de2bi(journal_num + 7, bit_count, 'left-msb');
base_seq_x2 = de2bi(journal_num + 1, bit_count, 'left-msb');
base_seq_y2 = de2bi(journal_num + 2, bit_count, 'left-msb');

fprintf('=== Параметры генерации ===\n');
fprintf('Номер в журнале: %d\n', journal_num);
fprintf('Длина последовательности: %d\n', seq_len);
fprintf('\nБазовые последовательности (5 бит):\n');
fprintf('x:  %s\n', num2str(base_seq_x));
fprintf('y:  %s\n', num2str(base_seq_y));
fprintf('x'': %s\n', num2str(base_seq_x2));
fprintf('y'': %s\n', num2str(base_seq_y2));

%% Полиномы обратной связи
poly_set1 = [5, 3];        % x^5 + x^3 + 1
poly_set2 = [5, 4, 2, 1];  % x^5 + x^4 + x^2 + x + 1

%% Генерация М-последовательностей
m_seq_x = create_m_sequence(base_seq_x, poly_set1);
m_seq_y = create_m_sequence(base_seq_y, poly_set2);
m_seq_x2 = create_m_sequence(base_seq_x2, poly_set1);
m_seq_y2 = create_m_sequence(base_seq_y2, poly_set2);

%% Визуализация М-последовательностей
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

%% Генерация последовательностей Голда
shift_val = find_best_shift(m_seq_x, m_seq_y);
gold_seq1 = xor(m_seq_x, circshift(m_seq_y, shift_val));
gold_seq2 = xor(m_seq_x2, circshift(m_seq_y2, shift_val));

fprintf('\n=== Последовательности Голда ===\n');
fprintf('Оптимальный сдвиг: %d\n', shift_val);
fprintf('\nGold Seq1 (первые 20 бит):\n');
fprintf('%d ', gold_seq1(1:20)');
fprintf('\n\nGold Seq2 (первые 20 бит):\n');
fprintf('%d ', gold_seq2(1:20)');

analyze_balance(gold_seq1, 'Последовательность Голда 1');
analyze_balance(gold_seq2, 'Последовательность Голда 2');

%% Визуализация последовательностей Голда
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

%% Автокорреляция
auto_corr_gold1 = zeros(seq_len, 1);
auto_corr_gold2 = zeros(seq_len, 1);

for shift = 0:seq_len-1
    shifted1 = circshift(gold_seq1, shift);
    shifted2 = circshift(gold_seq2, shift);
    
    auto_corr_gold1(shift+1) = calc_normalized_corr(gold_seq1, shifted1);
    auto_corr_gold2(shift+1) = calc_normalized_corr(gold_seq2, shifted2);
end

%% Взаимная корреляция
cross_corr = zeros(seq_len, 1);
for shift = 0:seq_len-1
    shifted2 = circshift(gold_seq2, shift);
    cross_corr(shift+1) = calc_normalized_corr(gold_seq1, shifted2);
end

%% Визуализация корреляционных функций
lag_range = -(seq_len-1):(seq_len-1);

figure('Position', [100, 100, 1000, 800]);

% Автокорреляция Seq1
subplot(3,2,1);
stem(0:seq_len-1, auto_corr_gold1, 'filled');
title('Автокорреляция Gold Seq1');
xlabel('Сдвиг'); ylabel('Корреляция');
grid on; ylim([-1, 1]);

% Автокорреляция Seq2
subplot(3,2,2);
stem(0:seq_len-1, auto_corr_gold2, 'filled');
title('Автокорреляция Gold Seq2');
xlabel('Сдвиг'); ylabel('Корреляция');
grid on; ylim([-1, 1]);

% Взаимная корреляция
subplot(3,2,[3,4]);
stem(-(seq_len-1):(seq_len-1), [flipud(cross_corr(2:end)); cross_corr], 'filled');
title('Взаимная корреляция Gold Seq1 и Gold Seq2');
xlabel('Сдвиг'); ylabel('Корреляция');
grid on; ylim([-1, 1]);

% Сравнение автокорреляций
subplot(3,2,5);
plot(0:seq_len-1, auto_corr_gold1, 'b-', 'LineWidth', 1.5);
hold on;
plot(0:seq_len-1, auto_corr_gold2, 'r--', 'LineWidth', 1.5);
title('Сравнение автокорреляционных функций');
xlabel('Сдвиг'); ylabel('Корреляция');
legend('Gold Seq1', 'Gold Seq2');
grid on; ylim([-1, 1]);

% Гистограмма корреляционных значений
subplot(3,2,6);
histogram(auto_corr_gold1, 20, 'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(auto_corr_gold2, 20, 'FaceColor', 'r', 'FaceAlpha', 0.5);
title('Распределение значений корреляции');
xlabel('Корреляция'); ylabel('Частота');
legend('Gold Seq1', 'Gold Seq2');
grid on;

%% Вывод статистики
fprintf('\n=== Корреляционный анализ ===\n');
fprintf('Автокорреляция Gold Seq1 (сдвиг 0): %.4f\n', auto_corr_gold1(1));
fprintf('Автокорреляция Gold Seq2 (сдвиг 0): %.4f\n', auto_corr_gold2(1));
fprintf('Взаимная корреляция (сдвиг 0): %.4f\n', cross_corr(1));

% Нахождение минимальной корреляции (кроме сдвига 0)
min_auto1 = min(auto_corr_gold1(2:end));
min_auto2 = min(auto_corr_gold2(2:end));
max_cross = max(abs(cross_corr(2:end)));

fprintf('\nМинимальная автокорреляция (ненулевые сдвиги):\n');
fprintf('  Gold Seq1: %.4f\n', min_auto1);
fprintf('  Gold Seq2: %.4f\n', min_auto2);
fprintf('Максимальная взаимная корреляция: %.4f\n', max_cross);

analyze_cycles(gold_seq1, 'Gold Seq1');
analyze_cycles(gold_seq2, 'Gold Seq2');

fprintf('\n=== Анализ завершен ===\n');

%% ========== ФУНКЦИИ (должны быть в конце файла) ==========

%% Функция генерации М-последовательности
function m_seq = create_m_sequence(init_state, poly_taps)
    reg_len = length(init_state);
    m_seq = zeros(2^reg_len - 1, 1);
    reg = init_state;
    
    for idx = 1:length(m_seq)
        % Выходной бит
        m_seq(idx) = reg(end);
        
        % Вычисление обратной связи
        fb_bit = 0;
        for tap = poly_taps
            fb_bit = xor(fb_bit, reg(tap));
        end
        
        % Сдвиг регистра
        reg = circshift(reg, 1);
        reg(1) = fb_bit;
    end
end

%% Функция поиска оптимального сдвига для последовательностей Голда
function optimal_shift = find_best_shift(seq1, seq2)
    sequence_length = length(seq1);  % Используем локальную переменную
    optimal_shift = 0;
    best_balance = sequence_length;  % Инициализируем максимально возможным значением
    
    for s = 0:sequence_length-1
        temp_gold = xor(seq1, circshift(seq2, s));
        ones_count = sum(temp_gold);
        
        % Ищем наиболее сбалансированную последовательность
        balance_diff = abs(2*ones_count - sequence_length);
        if balance_diff < best_balance
            best_balance = balance_diff;
            optimal_shift = s;
        end
    end
end

%% Функция анализа баланса
function analyze_balance(sequence, name)
    ones_cnt = sum(sequence);
    zeros_cnt = length(sequence) - ones_cnt;
    
    fprintf('\n%s:\n', name);
    fprintf('  Единиц: %d (%.1f%%)\n', ones_cnt, 100*ones_cnt/length(sequence));
    fprintf('  Нулей:  %d (%.1f%%)\n', zeros_cnt, 100*zeros_cnt/length(sequence));
    fprintf('  Баланс: %+d\n', ones_cnt - zeros_cnt);
end

%% Функция вычисления нормированной корреляции
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

%% Функция анализа циклических свойств
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
    
    % Последний цикл
    if ~isKey(cycle_counts, cycle_len)
        cycle_counts(cycle_len) = 0;
    end
    cycle_counts(cycle_len) = cycle_counts(cycle_len) + 1;
    
    % Вывод результатов
    keys_arr = cell2mat(keys(cycle_counts));
    [sorted_keys, idx] = sort(keys_arr);
    
    for k = 1:length(sorted_keys)
        key_val = sorted_keys(k);
        fprintf('  Длина %d: %d циклов\n', key_val, cycle_counts(key_val));
    end
end