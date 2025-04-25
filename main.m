%% Cài đặt hệ thống
fc = 2.4e9;  % Tần số sóng mang 2.4 GHz
c = 3e8;
lambda = c / fc;

% Mảng anten 3 phần tử cách nhau λ/2
antenna_pos = [0 0.5 1.0] * lambda / 2;
M = length(antenna_pos);

%% Cài đặt RRC pulse
Tsymbol = 1e-6;       % thời gian 1 symbol (1 µs)
sps = 8;              % số mẫu trên mỗi symbol
T0 = Tsymbol / sps;   % thời gian lấy mẫu (phụ thuộc vào sps)

beta = 0.25;          % Roll-off factor
span = 6;             % số symbol mỗi xung RRC trải dài
p = rcosdesign(beta, span, sps, 'sqrt');  % xung RRC
N = length(p);        % số mẫu xung

t_vec = (0:N-1)*T0;   % trục thời gian tương ứng với p

%% Tập giá trị quét
tau_grid = 0:T0:50e-9;
phi_grid = -90:5:90;
fd_grid = -500:50:500;  % Doppler ±500 Hz

%% Khởi tạo danh sách tham số cho 3 đường truyền (giả lập)
theta_list(1) = struct('tau', 10e-9, 'phi', 30,  'fd', 100,  'alpha', 1+1j);
theta_list(2) = struct('tau', 15e-9, 'phi', -20, 'fd', -50,  'alpha', 0.8);
theta_list(3) = struct('tau', 25e-9, 'phi', 70,  'fd', 200,  'alpha', 0.5-0.3j);
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
max_steps = 60;  % ví dụ: 20 chu kỳ cho 3 đường truyền
for mu = 1:max_steps
    [theta_list, l_updated] = sage_step(r, theta_list, mu, ...
        p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc);

    fprintf('Step %2d: updated path #%d\n', mu, l_updated);
end
