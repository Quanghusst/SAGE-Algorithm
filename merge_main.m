load("pos.mat"); load("ketquaCoPhi.mat");
%% Cài đặt hệ thống
fc = 2e9;  % Tần số sóng mang 2 GHz
c = 3e8;
lambda = c / fc;

% Mảng anten 3 phần tử cách nhau λ/2
% antenna_pos = [0 0.5 1.0] * lambda / 2;
antenna_pos =  [0 sqrt(sum( diff(pos_centers(:,11:20), 1, 2).^2, 1))];
M = length(antenna_pos);

%% Cài đặt RRC pulse
Tsymbol = 1e-6;       % thời gian 1 symbol (1 µs)
sps = 8;              % số mẫu trên mỗi symbol
T0 = Tsymbol / sps;   % thời gian lấy mẫu (phụ thuộc vào sps)

beta = 0.6;          % Roll-off factor
span = 6;             % số symbol mỗi xung RRC trải dài
p = rcosdesign(beta, span, sps, 'sqrt');  % xung RRC
    u = generatePulse(md, A(i, 1), tau, 0); % Tạo xung RRC với trễ A(i, 1);
N = length(p);        % số mẫu xung

t_vec = (0:N-1)*T0;   % trục thời gian tương ứng với p

%% Tập giá trị quét
%tau_grid = 0:T0:50e-9;
tau_grid = 0:4.6414e-12:14999*4.6414e-12; % Khoảng thời gian 0 đến 7e-8 giây
phi_grid = -90:5:90;
fd_grid = -500:50:500;  % Doppler ±500 Hz

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
r = zeros(M, N);
for l = 1:L
    tau_l = theta_list(l).tau;
    phi_l = theta_list(l).phi;
    fd_l  = theta_list(l).fd;
    alpha_l = theta_list(l).alpha;

    % Vector hướng anten
    phi_rad = deg2rad(phi_l);
    steering = exp(-1j * 2*pi/lambda * antenna_pos * sin(phi_rad)).';

    % Doppler và shaping pulse
    doppler_term = exp(-1j * 2*pi * fd_l * t_vec);
    pulse_shifted = interp1(t_vec, p, t_vec - tau_l, 'linear', 0);

    s_l = alpha_l * (steering) * (pulse_shifted .* doppler_term);

    % Cộng vào tín hiệu thu
    r = r + s_l;
end

% Cộng thêm nhiễu Gaussian trắng
noise_power = 0.01;
r = r + sqrt(noise_power/2) * (randn(M, N) + 1j*randn(M, N));

%% Chạy thuật toán SAGE (phiên bản step-by-step)
max_steps = 6;  % ví dụ: 20 chu kỳ cho 3 đường truyền
for mu = 1:max_steps
    [theta_list, l_updated] = sage_step(r, theta_list, mu, ...
        p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc);

    fprintf('Step %2d: updated path #%d\n', mu, l_updated);
end
