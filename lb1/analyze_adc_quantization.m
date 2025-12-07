function analyze_adc_quantization(signal, fs)
    
    function quantized = quantize_signal(signal, bits)
        max_val = max(abs(signal));
        levels = 2^bits - 1;
        quantized = round((signal + max_val) / (2*max_val) * levels);
        quantized = (quantized / levels) * 2*max_val - max_val;
    end

    bits_range = [3, 4, 5, 6];
    errors = zeros(size(bits_range));
    
    figure;
    
    
    N = length(signal);
    Y_original = fft(signal);
    freq = (0:N-1)*(fs/N);
    amp_original = (2/N) * abs(Y_original(1:N/2));
    
    subplot(2,1,1);
    plot(freq(1:N/2), amp_original, 'k-', 'LineWidth', 2, 'DisplayName', 'Оригинал');
    hold on;
    
    colors = ['r', 'g', 'b', 'm'];
    for i = 1:length(bits_range)
        bits = bits_range(i);
        quantized = quantize_signal(signal, bits);
        
        
        errors(i) = mean((signal - quantized).^2);
        
        
        Y_quant = fft(quantized);
        amp_quant = (2/N) * abs(Y_quant(1:N/2));
        
        plot(freq(1:N/2), amp_quant, [colors(i) '--'], ...
            'DisplayName', sprintf('%d бит, ошибка: %.6f', bits, errors(i)));
    end
    
    title('Влияние разрядности АЦП на спектр');
    xlabel('Частота (Гц)'); ylabel('Амплитуда');
    legend; grid on; hold off;
    
    
    subplot(2,1,2);
    semilogy(bits_range, errors, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Зависимость ошибки квантования от разрядности АЦП');
    xlabel('Разрядность АЦП (бит)'); ylabel('Среднеквадратическая ошибка');
    grid on;
    
    
    disp('ОШИБКИ КВАНТОВАНИЯ:');
    disp('-------------------');
    for i = 1:length(bits_range)
        fprintf('%d бит: %.8f\n', bits_range(i), errors(i));
    end
end