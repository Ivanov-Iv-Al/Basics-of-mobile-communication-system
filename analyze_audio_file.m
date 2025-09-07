function analyze_audio_file()
    %Чтение аудиофайла
    filename = 'voice.wav';
    
    if exist(filename, 'file') == 2
        %Чтение файла
        [y, Fs] = audioread(filename);
        fprintf('Файл %s успешно прочитан\n', filename);
        
        
        if size(y, 2) > 1
            y = mean(y, 2);
            disp('Преобразовано stereo → mono');
        end
        
        %Расчет частоты дискретизации
        duration = length(y) / Fs;
        calculated_Fs = length(y) / duration;
        
        fprintf('Заявленная Fs: %d Гц\n', Fs);
        fprintf('Рассчитанная Fs: %.0f Гц\n', calculated_Fs);
        fprintf('Длительность: %.2f секунд\n', duration);
        fprintf('Количество отсчетов: %d\n', length(y));
        
        %Прореживание и воспроизведение
        downsample_factor = 10;
        y_downsampled = downsample(y, downsample_factor);
        new_Fs = Fs / downsample_factor;
        
        fprintf('Коэффициент прореживания: %d\n', downsample_factor);
        fprintf('Новая Fs: %.0f Гц\n', new_Fs);
        
        % Визуализация
        figure;
        subplot(2,2,1);
        t_original = (0:length(y)-1)/Fs;
        plot(t_original(1:2000), y(1:2000), 'b-');
        title('Оригинальный сигнал');
        xlabel('Время (с)'); ylabel('Амплитуда'); grid on;
        
        subplot(2,2,2);
        t_downsampled = (0:length(y_downsampled)-1)/new_Fs;
        plot(t_downsampled(1:200), y_downsampled(1:200), 'ro-');
        title('Прореженный сигнал');
        xlabel('Время (с)'); ylabel('Амплитуда'); grid on;
        
        %Анализ спектра
        % Спектр оригинального сигнала
        N = length(y);
        Y = fft(y);
        freq = (0:N-1)*(Fs/N);
        amp_spectrum_original = (2/N) * abs(Y(1:N/2));
        
        subplot(2,2,3);
        plot(freq(1:N/2), amp_spectrum_original, 'b-');
        title('Спектр оригинального сигнала');
        xlabel('Частота (Гц)'); ylabel('Амплитуда'); grid on;
        xlim([0 10000]);
        
        % Спектр прореженного сигнала
        N_ds = length(y_downsampled);
        Y_ds = fft(y_downsampled);
        freq_ds = (0:N_ds-1)*(new_Fs/N_ds);
        amp_spectrum_ds = (2/N_ds) * abs(Y_ds(1:N_ds/2));
        
        subplot(2,2,4);
        plot(freq_ds(1:N_ds/2), amp_spectrum_ds, 'r-');
        title('Спектр прореженного сигнала');
        xlabel('Частота (Гц)'); ylabel('Амплитуда'); grid on;
        xlim([0 10000]);
        
        % Определение ширины спектра
        threshold = 0.01 * max(amp_spectrum_original);
        spectrum_width_original = max(freq(amp_spectrum_original > threshold));
        spectrum_width_ds = max(freq_ds(amp_spectrum_ds > threshold));
        
        fprintf('Ширина спектра оригинального: %.0f Гц\n', spectrum_width_original);
        fprintf('Ширина спектра прореженного: %.0f Гц\n', spectrum_width_ds);
        
       
        choice = input('Воспроизвести звук? (1-оригинальный, 2-прореженный, 0-нет): ');
        if choice == 1
            sound(y, Fs);
        elseif choice == 2
            sound(y_downsampled, new_Fs);
        end
        
    else
        fprintf('Файл %s не найден!\n', filename);
        fprintf('Запишите аудиофайл с помощью Audacity или другого редактора\n');
    end
end