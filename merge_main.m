load("pos.mat"); load("ketQuaKhoiTao"); load("IR_12.mat");
%% Cài đặt hệ thống
fc = 2e9;  % Tần số sóng mang 2 GHz
c = 3e8;
lambda = c / fc;

md.type = 'RRC'; % Xung RRC
md.Tp = 0.5e-9; % Chu kỳ xung (băng thông 2GHz)
md.beta = 0.6; %  Hệ số giảm tốc 

% Mảng anten 10 phần tử cách nhau lung tung ;
%antenna_pos =  [0 sqrt(sum( diff(pos_centers(:,11:20), 1, 2).^2, 1))];
antenna_pos = pos_centers(:, 11:20);
M = length(antenna_pos);
%% Tập giá trị quét
T0 = 4.6414e-12;
tau_grid = 0:T0:14999*T0; % Khoảng thời gian 0 đến 7e-8 giây
phi_grid = -90:5:90; % nên để radian
fd_grid = -500:50:500;  % Doppler ±500 Hz
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
max_steps = 6;  % ví dụ: 20 chu kỳ cho 3 đường truyền
for mu = 1:max_steps
    [theta_list, l_updated] = sage_step(r, theta_list, mu, ...
        p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc);

    fprintf('Step %2d: updated path #%d\n', mu, l_updated);
end
