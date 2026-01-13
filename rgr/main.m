jor_num = 10;
char_bits = 8;
samp_per_bit = 10;

user_name = input("Мне нужно имя, тока на английском: ", "s");

bin_seq = string_to_binary(user_name);

fprintf('Битовая последовательность имени:\n');
disp(bin_seq');

figure;
stairs(0:length(bin_seq)-1, bin_seq, 'LineWidth', 1.5);
xlabel('Время (биты)');
ylabel('Значение бита');
title('Битовое представление сообщения');
grid on;

crc_poly = [1,0,0,0,0,0,1,1,1];
check_bits = calculate_check_sum(bin_seq, crc_poly);

fprintf('Контрольная сумма CRC:\n');
disp(check_bits');

msg_with_check = [bin_seq; check_bits];

base_seq1 = dec2bin(jor_num, char_bits) - '0';
base_seq2 = dec2bin(jor_num + 7, char_bits) - '0';

feed_poly1 = [5, 3];
feed_poly2 = [5, 4, 2, 1];

m_seq_one = create_m_sequence(base_seq1, feed_poly1);
m_seq_two = create_m_sequence(base_seq2, feed_poly2);

sync_sequence = xor(m_seq_one, m_seq_two);
complete_frame = [sync_sequence; msg_with_check];

fprintf('Синхропоследовательность Голда (%d бит):\n', length(sync_sequence));
disp(sync_sequence');

figure;
stem(0:length(sync_sequence)-1, sync_sequence, 'filled', 'MarkerSize', 5);
xlabel('Позиция');
ylabel('Бит');
title('Последовательность Голда для синхронизации');
grid on;

sample_data = expand_bits(complete_frame, samp_per_bit);
impulse_shape = ones(1, samp_per_bit);
shaped_signal = filter_signal(sample_data, impulse_shape);

figure;
plot(0:length(shaped_signal)-1, shaped_signal);
xlabel('Отсчеты');
ylabel('Амплитуда');
title('Сформированный сигнал для передачи');
grid on;

empty_channel = zeros(2*length(shaped_signal), 1);
insert_point = str2double(input("Введите позицию для вставки пакета: ", "s"));

while true
    if insert_point + length(shaped_signal) > length(empty_channel)
        insert_point = str2double(input("Недопустимая позиция. Введите снова: ", "s"));
        continue;
    end
    break;
end

for k = 1:length(shaped_signal)
    empty_channel(insert_point + k) = shaped_signal(k);
end

noise_level = str2double(input("Введите дисперсию шума: ", "s"));
channel_noise = normrnd(0, noise_level, length(empty_channel), 1);
received_signal = empty_channel + channel_noise;

figure;
plot(0:length(channel_noise)-1, channel_noise);
xlabel('Отсчеты');
ylabel('Уровень шума');
title('Шум в канале связи');
grid on;

figure;
plot(0:length(received_signal)-1, received_signal);
xlabel('Отсчеты');
ylabel('Амплитуда');
title('Принятый сигнал с шумом');
grid on;

ref_sync_signal = filter_signal(expand_bits(sync_sequence, samp_per_bit), impulse_shape);
correlation_values = find_sync_position(received_signal, ref_sync_signal);

figure;
plot(0:length(correlation_values)-1, correlation_values);
xlabel('Отсчеты');
ylabel('Корреляция');
title('Функция корреляции');
grid on;

sync_start = find_max_correlation(correlation_values);
fprintf('Начало синхросигнала: %d\n', sync_start);

extracted_segment = received_signal(sync_start:sync_start+length(sample_data)-1);
detected_bits = convert_samples(extracted_segment, samp_per_bit);
detected_bits = detected_bits(length(sync_sequence)+1:end);

figure;
stairs(0:length(detected_bits)-1, detected_bits, 'LineWidth', 1.5);
xlabel('Биты');
ylabel('Значение');
title('Детектированные биты');
grid on;

rx_message_bits = detected_bits(1:end-length(check_bits));
rx_check_bits = detected_bits(end-length(check_bits)+1:end);
verify_check = calculate_check_sum(rx_message_bits, crc_poly);

if isequal(verify_check, rx_check_bits)
    fprintf('Передача успешна, ошибок нет\n');
    recovered_text = binary_to_string(rx_message_bits);
    fprintf('Восстановленный текст: %s\n', recovered_text);
else
    fprintf('Обнаружены ошибки передачи\n');
end

figure;
[spectrum1, freq1] = compute_signal_spectrum(shaped_signal, samp_per_bit);
plot(freq1, 10*log10(spectrum1), 'LineWidth', 1.5);
hold on;

long_symbols = filter_signal(expand_bits(complete_frame, samp_per_bit*2), ones(1, samp_per_bit*2));
[spectrum2, freq2] = compute_signal_spectrum(long_symbols, samp_per_bit*2);
plot(freq2, 10*log10(spectrum2), 'LineWidth', 1.5);

short_symbols = filter_signal(expand_bits(complete_frame, floor(samp_per_bit/2)), ones(1, floor(samp_per_bit/2)));
[spectrum3, freq3] = compute_signal_spectrum(short_symbols, floor(samp_per_bit/2));
plot(freq3, 10*log10(spectrum3), 'LineWidth', 1.5);

xlabel('Нормированная частота');
ylabel('Спектральная плотность (дБ)');
title('Спектры при разной длительности символов');
legend({'Базовая', 'Удлиненные', 'Короткие'}, 'Location', 'best');
grid on;
hold off;

sigmas = 0.05:0.05:0.5;
samp_factors = [5, 10, 20];
num_trials = 50;
max_corr_results = zeros(length(sigmas), length(samp_factors));

for sigma_idx = 1:length(sigmas)
    sigma = sigmas(sigma_idx);
    
    for samp_idx = 1:length(samp_factors)
        current_samp_factor = samp_factors(samp_idx);
        
        ref_sync_samples = expand_bits(sync_sequence, current_samp_factor);
        impulse_shape_current = ones(1, current_samp_factor);
        ref_sync_signal = filter_signal(ref_sync_samples, impulse_shape_current);
        
        max_corr_values = zeros(num_trials, 1);
        
        for trial = 1:num_trials
            test_signal_length = length(ref_sync_signal) * 2;
            test_signal = zeros(test_signal_length, 1);
            
            insert_pos = randi([1, length(test_signal) - length(ref_sync_signal)]);
            test_signal(insert_pos:insert_pos+length(ref_sync_signal)-1) = ref_sync_signal;
            
            noise = normrnd(0, sigma, size(test_signal));
            noisy_signal = test_signal + noise;
            
            correlation = find_sync_position(noisy_signal, ref_sync_signal);
            
            [max_val, ~] = max(correlation);
            max_corr_values(trial) = max_val;
        end
        
        max_corr_results(sigma_idx, samp_idx) = mean(max_corr_values);
    end
end

figure('Position', [100, 100, 900, 600]);
colors = lines(length(samp_factors));
hold on;

for samp_idx = 1:length(samp_factors)
    plot(sigmas, max_corr_results(:, samp_idx), ...
         'LineWidth', 2.5, ...
         'Marker', 'o', ...
         'MarkerSize', 8, ...
         'MarkerFaceColor', colors(samp_idx,:), ...
         'DisplayName', sprintf('%d отсчетов/бит', samp_factors(samp_idx)));
end

xlabel('Дисперсия шума (σ)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Максимум корреляции', 'FontSize', 12, 'FontWeight', 'bold');
title('Зависимость максимума корреляции от уровня шума для разной длительности символа', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on;
box on;

hold off;

function bin_data = string_to_binary(input_text)
    text_len = length(input_text);
    total_bits = text_len * 8;
    bin_data = zeros(total_bits, 1);

    for char_idx = 1:text_len
        ascii_val = double(input_text(char_idx));
        for bit_idx = 1:8
            position = (char_idx-1)*8 + bit_idx;
            bin_data(position) = bitand(bitshift(int8(ascii_val), -(8-bit_idx)), 1);
        end
    end
end

function check_result = calculate_check_sum(data_bits, poly_vector)
    data_vect = data_bits(:);
    check_len = length(poly_vector) - 1;
    ext_data = [data_vect; zeros(check_len, 1)];

    for idx = 1:length(data_vect)
        if ext_data(idx) == 1
            ext_data(idx:idx+check_len) = xor(ext_data(idx:idx+check_len), poly_vector(:));
        end
    end

    check_result = ext_data(end-check_len+1:end);
end

function m_sequence = create_m_sequence(init_reg, feedback_pattern)
    reg_size = length(init_reg);
    seq_size = 2^reg_size - 1;
    m_sequence = zeros(seq_size, 1);
    shift_reg = init_reg(:);

    for cycle = 1:seq_size
        fb_value = 0;
        for tap_pos = feedback_pattern
            fb_value = xor(fb_value, shift_reg(tap_pos));
        end

        m_sequence(cycle) = shift_reg(end);
        shift_reg = circshift(shift_reg, 1);
        shift_reg(1) = fb_value;
    end
end

function samples_out = expand_bits(bit_stream, factor)
    samples_out = zeros(length(bit_stream) * factor, 1);

    for n = 1:length(bit_stream)
        sample_pos = (n-1)*factor + 1;
        samples_out(sample_pos) = bit_stream(n);
    end
end

function filtered_out = filter_signal(input_signal, kernel)
    sig_len = length(input_signal);
    ker_len = length(kernel);
    filtered_out = zeros(sig_len, 1);

    for n = 1:sig_len
        accum = 0;
        for m = 1:ker_len
            if n - m > 0
                accum = accum + input_signal(n-m) * kernel(m);
            end
        end
        filtered_out(n) = accum;
    end
end

function correlation_map = find_sync_position(signal, pattern)
    sig_len = length(signal);
    pat_len = length(pattern);
    correlation_map = zeros(sig_len, 1);

    pattern_norm = sqrt(sum(pattern.^2));

    for pos = 1:sig_len
        window_end = min(pos + pat_len - 1, sig_len);
        window_data = signal(pos:window_end);
        
        if length(window_data) == pat_len
            signal_norm = sqrt(sum(window_data.^2));
            correlation_map(pos) = sum(window_data .* pattern) / (signal_norm * pattern_norm);
        end
    end
end

function max_pos = find_max_correlation(corr_vector)
    max_val = -100;
    max_pos = 1;
    
    for idx = 1:length(corr_vector)
        if corr_vector(idx) > max_val
            max_val = corr_vector(idx);
            max_pos = idx;
        end
    end
end

function bit_stream = convert_samples(samples, factor)
    bit_count = floor(length(samples) / factor);
    bit_stream = zeros(bit_count, 1);

    for k = 1:bit_count
        start_pos = (k-1)*factor + 1;
        end_pos = k*factor;
        segment = samples(start_pos:end_pos);
        bit_stream(k) = mean(segment) > 0.5;
    end
end

function text_out = binary_to_string(bit_array)
    char_count = length(bit_array) / 8;
    text_out = "";
    
    for ch = 1:char_count
        start_bit = (ch-1)*8 + 1;
        end_bit = ch*8;
        char_bits = bit_array(start_bit:end_bit);
        
        char_val = 0;
        for bit = 1:8
            char_val = char_val + char_bits(bit) * 2^(8-bit);
        end
        
        text_out = text_out + char(char_val);
    end
end


function [spectrum, freq] = compute_signal_spectrum(signal, samp_per_bit)
    N = length(signal);
    window = hann(N);
    signal = signal .* window;
    fft_result = fft(signal, N);
    spectrum = abs(fft_result(1:floor(N/2)+1)).^2 / N;
    Fs = 1000 * samp_per_bit;
    freq = (0:floor(N/2)) * Fs / N;
end

