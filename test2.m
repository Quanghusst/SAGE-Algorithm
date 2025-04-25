% Thông số cơ bản
rolloff = 0.25;      % Hệ số roll-off của bộ lọc RRC
span = 6;            % Độ dài xung (tính bằng số ký tự symbol)
sps = 8;             % Số mẫu trên mỗi symbol (samples per symbol)

% Tạo bộ lọc Root Raised Cosine
rrcFilter = rcosdesign(rolloff, span, sps, 'sqrt');

% Hiển thị xung RRC
figure;
t = (-span/2:1/sps:span/2);
plot(t, rrcFilter, 'LineWidth',1.5);
title('Xung Root Raised Cosine (RRC)');
xlabel('Thời gian (symbol)');
ylabel('Biên độ');
grid on;

% Giả lập truyền đi một chuỗi symbol ngẫu nhiên
numSymbols = 100;
data = randi([0 1], numSymbols, 1)*2 - 1; % Tạo chuỗi bit ngẫu nhiên {-1,1}

% Tạo tín hiệu phát đi bằng cách lấy mẫu symbol theo RRC
signal_tx = upfirdn(data, rrcFilter, sps);

% Trục thời gian cho tín hiệu phát
figure;
time_tx = (0:length(signal_tx)-1)/sps;
plot(time_tx, signal_tx);
title('Tín hiệu truyền đi qua bộ lọc RRC');
xlabel('Thời gian (symbol)');
ylabel('Biên độ');
grid on;
