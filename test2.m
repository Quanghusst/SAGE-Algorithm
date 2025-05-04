%% Khởi tạo thông số
fc = 2.4e9;               % Tần số sóng mang (Hz)
c = 3e8;
lambda = c / fc;

antenna_pos = [0 0.5 1.0] * lambda / 2;  % Mảng anten tuyến tính
M = length(antenna_pos);                % Số phần tử anten

Tsymbol = 1e-6;          % Thời gian 1 symbol (1 µs)
sps = 108;               % Số mẫu trên mỗi symbol
T0 = Tsymbol / sps;      % Thời gian lấy mẫu (s)
span = 140;              % Số symbol mà xung trải dài

%% Sinh xung RRC (Root Raised Cosine)
beta = 0.25;  % Hệ số roll-off
p = rcosdesign(beta, span, sps, 'sqrt');  % xung RRC
N = length(p);                            % Tổng số mẫu của xung
t_vec = (0:N-1) * T0;                     % Trục thời gian

% Xem xung
figure;
plot(t_vec * 1e6, p);
xlabel('Thời gian (µs)');
ylabel('Biên độ');
title('Xung RRC với sps = 108, span = 140');

%% Mã hóa dãy tín hiệu truyền
numSymbols = 100;                % Số lượng symbol muốn truyền
dataBits = randi([0 1], 1, numSymbols);  % Sinh dãy bit ngẫu nhiên
symbols = 2*dataBits - 1;        % Điều chế BPSK

% Phát tín hiệu bằng cách nhân từng symbol với xung RRC
tx_signal = upfirdn(symbols, p, sps);  % Tín hiệu đã được lọc

% Nếu bạn cần chèn sóng mang:
t = (0:length(tx_signal)-1)*T0;
carrier = cos(2*pi*fc*t);
tx_passband = tx_signal .* carrier;

% Vẽ tín hiệu truyền
figure;
plot(t(1:1000)*1e6, tx_signal(1:1000));
xlabel('Thời gian (µs)');
ylabel('Biên độ');
title('Tín hiệu truyền baseband');

