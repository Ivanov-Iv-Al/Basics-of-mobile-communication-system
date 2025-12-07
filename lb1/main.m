
clear all; close all; clc;
disp('=' * 60);
disp('ПОЛНЫЙ АНАЛИЗ СИГНАЛА И ОБРАБОТКА АУДИО');
disp('=' * 60);


disp('1. ГЕНЕРАЦИЯ И ВИЗУАЛИЗАЦИЯ СИГНАЛА');
disp('-' * 40);

f = 5;  
t_continuous = linspace(0, 1, 10000);  

% Сигнал y(t) = cos(6πft) + sin(4πft) + sin(10t)
y_continuous = cos(6*pi*f*t_continuous) + sin(4*pi*f*t_continuous) + sin(10*t_continuous);

figure;
plot(t_continuous, y_continuous, 'b-', 'LineWidth', 1.5);
title('Непрерывный сигнал: y(t) = cos(6πft) + sin(4πft) + sin(10t)');
xlabel('Время (с)');
ylabel('Амплитуда');
grid on;
axis tight;


disp('2. ОПРЕДЕЛЕНИЕ МАКСИМАЛЬНОЙ ЧАСТОТЫ');
disp('-' * 40);



max_freq = 15;  % Гц
fprintf('Максимальная частота в спектре: %d Гц\n', max_freq);


disp('3. ЧАСТОТА ДИСКРЕТИЗАЦИИ ПО ТЕОРЕМЕ КОТЕЛЬНИКОВА');
disp('-' * 40);

min_fs = 2 * max_freq;
fprintf('Минимальная частота дискретизации: %d Гц\n', min_fs);


disp('4. ОЦИФРОВКА СИГНАЛА');
disp('-' * 40);

fs = 30;  % Гц
t_discrete = linspace(0, 1, fs + 1);
y_discrete = cos(6*pi*f*t_discrete) + sin(4*pi*f*t_discrete) + sin(10*t_discrete);

fprintf('Частота дискретизации: %d Гц\n', fs);
fprintf('Количество отсчетов: %d\n', length(y_discrete));
disp('Отсчеты:');
disp(y_discrete');


disp('5. ДИСКРЕТНОЕ ПРЕОБРАЗОВАНИЕ ФУРЬЕ');
disp('-' * 40);

N = length(y_discrete);
Y = fft(y_discrete);
frequencies = (0:N-1)*(fs/N);  % Частотная ось

amplitude_spectrum = (2/N) * abs(Y(1:N/2));

figure;
plot(frequencies(1:N/2), amplitude_spectrum, 'r-', 'LineWidth', 2);
title('Амплитудный спектр (Fs = 30 Гц)');
xlabel('Частота (Гц)');
ylabel('Амплитуда');
grid on;


threshold = 0.05 * max(amplitude_spectrum);
significant_indices = amplitude_spectrum > threshold;
spectrum_width = max(frequencies(significant_indices));
fprintf('Ширина спектра: %.2f Гц\n', spectrum_width);


memory_size_bytes = length(y_discrete) * 8;  % double
memory_size_kb = memory_size_bytes / 1024;
fprintf('Объем памяти: %d байт (%.2f КБ)\n', memory_size_bytes, memory_size_kb);


disp('6. ВОССТАНОВЛЕНИЕ СИГНАЛА');
disp('-' * 40);

figure;
plot(t_continuous, y_continuous, 'b-', 'LineWidth', 2, 'DisplayName', 'Оригинальный');
hold on;
plot(t_discrete, y_discrete, 'ro-', 'MarkerSize', 6, 'DisplayName', 'Восстановленный');
title('Сравнение оригинального и восстановленного сигналов');
xlabel('Время (с)');
ylabel('Амплитуда');
legend;
grid on;
hold off;


disp('7. УВЕЛИЧЕНИЕ ЧАСТОТЫ ДИСКРЕТИЗАЦИИ В 4 РАЗА');
disp('-' * 40);

fs_high = 4 * fs;  % 120 Гц
t_discrete_high = linspace(0, 1, fs_high + 1);
y_discrete_high = cos(6*pi*f*t_discrete_high) + sin(4*pi*f*t_discrete_high) + sin(10*t_discrete_high);

fprintf('Новая частота дискретизации: %d Гц\n', fs_high);
fprintf('Количество отсчетов: %d\n', length(y_discrete_high));


N_high = length(y_discrete_high);
Y_high = fft(y_discrete_high);
frequencies_high = (0:N_high-1)*(fs_high/N_high);
amplitude_spectrum_high = (2/N_high) * abs(Y_high(1:N_high/2));

figure;
plot(frequencies_high(1:N_high/2), amplitude_spectrum_high, 'g-', 'LineWidth', 2);
title('Амплитудный спектр (Fs = 120 Гц)');
xlabel('Частота (Гц)');
ylabel('Амплитуда');
grid on;


figure;
plot(t_continuous, y_continuous, 'b-', 'LineWidth', 1, 'DisplayName', 'Оригинальный');
hold on;
plot(t_discrete_high, y_discrete_high, 'g-', 'LineWidth', 1, 'DisplayName', 'Восстановленный (120 Гц)');
title('Сравнение с повышенной частотой дискретизации');
xlabel('Время (с)');
ylabel('Амплитуда');
legend;
grid on;
hold off;


disp('8-12. АНАЛИЗ АУДИОФАЙЛА');
disp('-' * 40);

analyze_audio_file();


disp('13. ВЛИЯНИЕ РАЗРЯДНОСТИ АЦП НА СПЕКТР');
disp('-' * 40);

analyze_adc_quantization(y_discrete, fs);

disp('=' * 60);
disp('АНАЛИЗ ЗАВЕРШЕН!');
disp('=' * 60);