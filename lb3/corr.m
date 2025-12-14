a = [6, 2, 3, -2, -4, -4, 1, 1];
b = [3, 1, 5, 0, -3, -4, 2, 3];
c = [-1, -1, 3, -9, 2, -8, 4, -4];

fprintf('Ненормированная корреляция:\n');
fprintf('     a      b      c\n');
fprintf('a  %5d  %5d  %5d\n', ...
    unnorm_corr(a,a), unnorm_corr(a,b), unnorm_corr(a,c));
fprintf('b  %5d  %5d  %5d\n', ...
    unnorm_corr(b,a), unnorm_corr(b,b), unnorm_corr(b,c));
fprintf('c  %5d  %5d  %5d\n', ...
    unnorm_corr(c,a), unnorm_corr(c,b), unnorm_corr(c,c));

fprintf('\nНормированная корреляция:\n');
fprintf('     a      b      c\n');
fprintf('a  %5.2f  %5.2f  %5.2f\n', ...
    norm_corr(a,a), norm_corr(a,b), norm_corr(a,c));
fprintf('b  %5.2f  %5.2f  %5.2f\n', ...
    norm_corr(b,a), norm_corr(b,b), norm_corr(b,c));
fprintf('c  %5.2f  %5.2f  %5.2f\n', ...
    norm_corr(c,a), norm_corr(c,b), norm_corr(c,c));

num_in_journal = 10;
f1 = num_in_journal;
f2 = num_in_journal + 4;
f3 = num_in_journal * 2 + 1;

s1 = @(t) cos(2 * pi * f1 * t);
s2 = @(t) cos(2 * pi * f2 * t);
s3 = @(t) cos(2 * pi * f3 * t);

a_signal = @(t) 5 * s1(t) + 4 * s2(t) + s3(t);
b_signal = @(t) 3 * s1(t) + s3(t);

Fs = 1000; 
t = 0:1/Fs:10; 
a_discrete = a_signal(t);
b_discrete = b_signal(t);

unnorm_result = unnorm_corr(a_discrete, b_discrete);
norm_result = norm_corr(a_discrete, b_discrete);

fprintf('Частоты: f1 = %d Гц, f2 = %d Гц, f3 = %d Гц\n', f1, f2, f3);
fprintf('Ненормированная корреляция: %.6f\n', unnorm_result);
fprintf('Нормированная корреляция: %.6f\n', norm_result);

figure('Position', [100, 100, 1200, 400]);
subplot(1,2,1);
plot(t(1:1000), a_discrete(1:1000), 'b-', 'LineWidth', 1.5);
hold on;
plot(t(1:1000), b_discrete(1:1000), 'r--', 'LineWidth', 1.5);
xlabel('Время, с');
ylabel('Амплитуда');
title('Сигналы a(t) и b(t) (первые 1000 отсчетов)');
legend('a(t) = 5s1+4s2+s3', 'b(t) = 3s1+s3');
grid on;

subplot(1,2,2);
plot(a_discrete(1:100), b_discrete(1:100), 'ko', 'MarkerSize', 6);
xlabel('a(t)');
ylabel('b(t)');
title('Фазовый портрет (первые 100 отсчетов)');
grid on;

function corr = norm_corr(a, b)
    if length(a) ~= length(b)
        error('Vectors must be same size');
    end
    norm_coef = sqrt(sum(a.^2)) * sqrt(sum(b.^2));
    corr = sum(a .* b) / norm_coef;
end

function corr = unnorm_corr(a, b)
    if length(a) ~= length(b)
        error('Vectors must be same size');
    end
    corr = sum(a .* b);
end