%% Q2
% تنظیمات اولیه
fs = 8000;  % نرخ نمونه‌برداری
block_size = 32;  % یا 64
bit_depths = [4, 8, 12];  % تعداد بیت‌های کوانتیزاسیون
frequencies = [200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 3500];  % فرکانس‌های مختلف برای آزمایش
rmse_values = zeros(length(frequencies), length(bit_depths));

% تعریف تابع پنجره سینوسی
window = sqrt(2) * sin(((1:block_size*2) - 0.5) * pi / (2 * block_size));

for f_idx = 1:length(frequencies)
    f = frequencies(f_idx);
    t = (0:fs-1)/fs;
    x = sin(2 * pi * f * t);  % تولید سیگنال سینوسی
    
    for b_idx = 1:length(bit_depths)
        b = bit_depths(b_idx);
        
        % پردازش MDCT
        X_mdct = mdct_transform1(x, block_size, window);
        X_quant = quantize1(X_mdct, b);  % کوانتیزه کردن
        x_reconstructed = imdct_transform1(X_quant, block_size, window);
        
        % محاسبه RMSE با تنظیم اندازه سیگنال بازسازی‌شده
        min_length = min(length(x), length(x_reconstructed));
        rmse_values(f_idx, b_idx) = sqrt(mean((x(1:min_length) - x_reconstructed(1:min_length)).^2));
        
        % نمایش در حوزه زمان و فرکانس برای چند فرکانس نمونه
        if f_idx <= 3 && b_idx == 1  % نمایش برای سه فرکانس اول
            figure;
            subplot(2,1,1); plot(x); hold on; plot(x_reconstructed, '--r');
            title(['Time Domain Signal (f = ', num2str(f), ' Hz)']); legend('Original', 'Reconstructed');
            subplot(2,1,2); plot(abs(fft(x))); hold on; plot(abs(fft(x_reconstructed)), '--r');
            title('Frequency Domain Comparison'); legend('Original', 'Reconstructed');
        end
    end
end

% رسم نمودار RMSE
figure;
plot(frequencies, rmse_values, 'o-');
xlabel('Frequency (Hz)'); ylabel('RMSE');
title('RMSE vs Frequency for Different Bit Depths');
legend(arrayfun(@(b) [num2str(b), ' bits'], bit_depths, 'UniformOutput', false));
grid on;



%% Q3
% تنظیمات اولیه
fs = 8000;  % نرخ نمونه‌برداری
block_size = 32;  % می‌تواند 64 هم باشد
bit_depths = [4, 8, 12];  % تعداد بیت‌های کوانتیزاسیون
frequencies = [200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 3500];  % فرکانس‌های مختلف برای آزمایش
windows = {'rect', 'hann', 'hamm', 'black', 'kaiser'}; % انواع تابع پنجره
rmse_values = zeros(length(frequencies), length(bit_depths), length(windows));

% پردازش سیگنال برای هر فرکانس، بیت و پنجره
for w_idx = 1:length(windows)
    win_type = windows{w_idx};
    
    % ساخت پنجره
    switch win_type
        case 'rect', window = ones(1, block_size * 2);
        case 'hann', window = hann(block_size * 2)';
        case 'hamm', window = hamming(block_size * 2)';
        case 'black', window = blackman(block_size * 2)';
        case 'kaiser', window = kaiser(block_size * 2, 5)'; % بتای ۵ برای Kaiser
    end
    
    for f_idx = 1:length(frequencies)
        f = frequencies(f_idx);
        t = (0:fs-1)/fs;
        x = sin(2 * pi * f * t);  % تولید سیگنال سینوسی

        for b_idx = 1:length(bit_depths)
            b = bit_depths(b_idx);

            % پردازش MDCT
            X_mdct = mdct_transform2(x, block_size, window);
            X_quant = quantize2(X_mdct, b);  % کوانتیزه کردن
            x_reconstructed = imdct_transform2(X_quant, block_size, window);

            % محاسبه RMSE
            min_length = min(length(x), length(x_reconstructed));
            rmse_values(f_idx, b_idx, w_idx) = sqrt(mean((x(1:min_length) - x_reconstructed(1:min_length)).^2));
        end
    end
end

% رسم نمودار RMSE برای هر پنجره
figure;
for w_idx = 1:length(windows)
    subplot(ceil(length(windows)/2), 2, w_idx);
    plot(frequencies, squeeze(rmse_values(:, :, w_idx)), 'o-');
    xlabel('Frequency (Hz)'); ylabel('RMSE');
    title(['RMSE vs Frequency - ', windows{w_idx}, ' Window']);
    legend(arrayfun(@(b) [num2str(b), ' bits'], bit_depths, 'UniformOutput', false));
    grid on;
end
%% Q4
% بارگذاری فایل WAV
[file, fs] = audioread('C:\Users\clover\Downloads\project_music.wav');
[num_samples, num_channels] = size(file);
block_size = 32;
bit_depths = [4, 8, 12];
windows = {'none', 'hann'};
rmse_values = zeros(length(bit_depths), length(windows), num_channels);

% پردازش برای هر کانال
for ch = 1:num_channels
    x = file(:, ch)';
    
    for w_idx = 1:length(windows)
        win_type = windows{w_idx};
        
        % تعریف تابع پنجره
        switch win_type
            case 'none', window = ones(1, block_size * 2);
            case 'hann', window = hann(block_size * 2)';
        end
        
        for b_idx = 1:length(bit_depths)
            b = bit_depths(b_idx);
            
            % پردازش MDCT
            X_mdct = mdct_transform3(x, block_size, window);
            X_quant = quantize3(X_mdct, b);
            x_reconstructed = imdct_transform3(X_quant, block_size, window);
            
            % محاسبه RMSE
            min_length = min(length(x), length(x_reconstructed));
            rmse_values(b_idx, w_idx, ch) = sqrt(mean((x(1:min_length) - x_reconstructed(1:min_length)).^2));
        end
    end
end

% نمایش جدول RMSE
for ch = 1:num_channels
    fprintf('Channel %d:\n', ch);
    disp(array2table(squeeze(rmse_values(:, :, ch)), 'VariableNames', windows, 'RowNames', cellstr(num2str(bit_depths'))));
end

% ذخیره فایل صوتی بازسازی شده
output_file = 'output.wav';
audiowrite(output_file, x_reconstructed', fs);


%% توابع کمکی
function X = mdct_transform2(x, N, window)
    x_padded = [zeros(1, N), x, zeros(1, N)];  % پدینگ
    num_blocks = floor(length(x_padded) / N) - 1;
    X = zeros(N, num_blocks);
    for k = 1:num_blocks
        block = x_padded((k-1)*N + (1:2*N)) .* window;
        X(:, k) = dct(block, N, 'Type', 4);  % MDCT
    end
end

function x_reconstructed = imdct_transform2(X, N, window)
    num_blocks = size(X, 2);
    x_reconstructed = zeros(1, (num_blocks+1)*N);
    for k = 1:num_blocks
        block = idct(X(:, k), N, 'Type', 4) .* window(1:N)';  % اعمال پنجره
        idx = (k-1)*N + (1:N);  
        x_reconstructed(idx) = x_reconstructed(idx) + block'; % همپوشانی
    end
    x_reconstructed = x_reconstructed(N+1:end);  
end

function X_quant = quantize2(X, b)
    X_max = max(abs(X(:)));
    levels = 2^b;
    X_quant = round((X / X_max) * (levels / 2)) * (X_max / (levels / 2));
end


function X = mdct_transform1(x, N, window)
    x_padded = [zeros(1, N), x, zeros(1, N)];  % پدینگ
    num_blocks = floor(length(x_padded) / N) - 1;
    X = zeros(N, num_blocks);
    for k = 1:num_blocks
        block = x_padded((k-1)*N + (1:2*N)) .* window;
        X(:, k) = dct(block, N, 'Type', 4);  % MDCT
    end
end

function x_reconstructed = imdct_transform1(X, N, window)
    num_blocks = size(X, 2);
    x_reconstructed = zeros(1, (num_blocks+1)*N);
    for k = 1:num_blocks
        block = idct(X(:, k), N, 'Type', 4) .* window(1:N)';  % اصلاح اعمال پنجره
        idx = (k-1)*N + (1:N);  % تصحیح ایندکس‌ها
        x_reconstructed(idx) = x_reconstructed(idx) + block'; % اعمال همپوشانی
    end
    x_reconstructed = x_reconstructed(N+1:end);  % حذف پدینگ اولیه
end

function X_quant = quantize1(X, b)
    X_max = max(abs(X(:)));
    levels = 2^b;
    X_quant = round((X / X_max) * (levels / 2)) * (X_max / (levels / 2));
end
function X = mdct_transform3(x, N, window)
    x_padded = [zeros(1, N), x, zeros(1, N)];
    num_blocks = floor(length(x_padded) / N) - 1;
    X = zeros(N, num_blocks);
    for k = 1:num_blocks
        block = x_padded((k-1)*N + (1:2*N)) .* window;
        X(:, k) = dct(block, N, 'Type', 4);
    end
end

function x_reconstructed = imdct_transform3(X, N, window)
    num_blocks = size(X, 2);
    x_reconstructed = zeros(1, (num_blocks+1)*N);
    for k = 1:num_blocks
        block = idct(X(:, k), N, 'Type', 4) .* window(1:N)';
        idx = (k-1)*N + (1:N);
        x_reconstructed(idx) = x_reconstructed(idx) + block';
    end
    x_reconstructed = x_reconstructed(N+1:end);
end

function X_quant = quantize3(X, b)
    X_max = max(abs(X(:)));
    levels = 2^b;
    X_quant = round((X / X_max) * (levels / 2)) * (X_max / (levels / 2));
end