load("pos.mat"); load("ketQuaKhoiTao"); load("IR_12.mat");
%% Cài đặt hệ thống
fc = 2e9;  % Tần số sóng mang 2 GHz
c = 3e8;
lambda = c / fc;

md.type = 'RRC'; % Xung RRC
md.Tp = 0.5e-9; % Chu kỳ xung (băng thông 2GHz)
md.beta = 0.6; %  Hệ số giảm tốc 

% Mảng anten 10 phần tử cách nhau lung tung ;
% Khoảng cách giữa các anten là [0 sqrt(sum( diff(pos_centers(:,11:20), 1, 2).^2, 1))];
antenna_pos = pos_centers(:, 11:20);
M = length(antenna_pos);
%% Tập giá trị quét
T0 = 4.6414e-12;
tau_grid = (0 : 15000-1) * T0; % Khoảng thời gian 0 đến 7e-8 giây
phi_grid = deg2rad(-90:5:90);
fd_grid = 500;  % Doppler ±500 Hz
%% Cài đặt RRC pulse
p = generatePulse(md, A(1, 1), tau_grid, 0).';

%% Khởi tạo danh sách tham số cho 25 đường truyền 
initializationZ = A(1:25, :);
theta_list = struct('tau', [], 'phi', [], 'fd', [], 'alpha', []); % Khởi tạo cấu trúc mẫu
for i = 1:size(initializationZ, 1)
    theta_list(i).tau = initializationZ(i, 1);
    theta_list(i).phi = initializationZ(i, 2);
    theta_list(i).fd = initializationZ(i, 6);
    theta_list(i).alpha = initializationZ(i, 7);
end
L = length(theta_list);

%% Tín hiệu thu r: tổng của các tín hiệu từ theta_list + nhiễu
r = IR_12.';

%% Chạy thuật toán SAGE (phiên bản step-by-step)
nCycles = 20;
max_steps = nCycles * L;   % Lặp 20 vòng cho 25 đường truyền => có 20 x 25 vòng lặp

for mu = 1:max_steps
    [theta_list, l_updated] = sage_step(r, theta_list, mu, ...
        p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc);

    fprintf('Step %2d: updated path #%d\n', mu, l_updated);
end
