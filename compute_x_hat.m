%% Bước E: loại bỏ ảnh hưởng các đường truyền khác
function x_hat = compute_x_hat(r, l, theta_list, p, T0, antenna_pos, fc)
% r: tín hiệu thu (MxN)
% l: chỉ số đường truyền hiện tại cần xử lý
% theta_list: danh sách tham số của tất cả đường truyền (struct array)
%             mỗi phần tử gồm: .tau, .phi, .fd, .alpha
% p: shaping pulse
% T0: khoảng lấy mẫu
% antenna_pos: vị trí anten (1xM)
% fc: tần số sóng mang

c = 3e8;
lambda = c / fc;
[M, N] = size(r);
t_vec = (0:N-1) * T0;

% Khởi tạo tín hiệu tổng suy diễn từ các đường khác
s_total = zeros(M, N);

for j = 1:length(theta_list)
    if j == l
        continue; % bỏ qua đường cần xử lý
    end
    
    tau_j = theta_list(j).tau;
    phi_j = theta_list(j).phi;
    fd_j  = theta_list(j).fd;
    alpha_j = theta_list(j).alpha;
    
    % Vector hướng anten
    phi_rad = deg2rad(phi_j);
    steering = exp(-1j*2*pi/lambda * antenna_pos * sin(phi_rad)).';
    
    % Doppler và shaping pulse
    doppler_term = exp(-1j*2*pi*fd_j * t_vec);
    pulse_shifted = interp1(t_vec, p, t_vec - tau_j, 'linear', 0);
    
    % Tín hiệu tái tạo
    s_j = alpha_j * (steering) * (pulse_shifted .* doppler_term);
    
    % Cộng vào tổng
    s_total = s_total + s_j;
end

% Tính x̂_l: tín hiệu còn lại cho đường truyền l
x_hat = r - s_total;

end
