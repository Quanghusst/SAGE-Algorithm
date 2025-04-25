%% M_Step Sage
function [tau_opt, phi_opt, fd_opt, alpha_opt] = sage_mstep_grid(x_hat, ...
    p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc)

% x_hat: tín hiệu đã trừ sóng khác (E-step)
% p: hàm shaping pulse (dạng vector)
% tau_grid, phi_grid, fd_grid: các giá trị cần quét
% T0: thời gian lấy mẫu
% antenna_pos: vị trí phần tử anten (vector 1xM)
% fc: tần số sóng mang

c = 3e8;                  % vận tốc ánh sáng
lambda = c / fc;         % bước sóng
[M, N] = size(x_hat);    % M: số anten, N: số mẫu

% Chuẩn bị dạng tín hiệu mẫu
t_vec = (0:N-1)*T0;

max_z = -inf;
tau_opt = 0; phi_opt = 0; fd_opt = 0;

%% Grid search cho tau, phi, fd
for tau = tau_grid
    for phi = phi_grid
        for fd = fd_grid
            % Vector hướng anten a(phi)
            phi_rad = deg2rad(phi);
            steering = exp(-1j*2*pi/lambda * antenna_pos * sin(phi_rad)).';
            
            % Doppler compensation
            doppler_term = exp(-1j*2*pi*fd * t_vec);
            
            % Pulse shift
            pulse_shifted = interp1(t_vec, p, t_vec - tau, 'linear', 0);
            
            % Tín hiệu tổng hợp
            s = (steering) * (pulse_shifted .* doppler_term);
            
            % Tính hàm z
            z = abs(sum(sum(conj(s) .* x_hat)));

            % Lưu tham số tối ưu nếu z lớn hơn
            if z > max_z
                max_z = z;
                tau_opt = tau;
                phi_opt = phi;
                fd_opt = fd;
            end
        end
    end
end

%% Sau khi tìm được bộ (tau, phi, fd) tốt nhất → tính alpha
phi_rad = deg2rad(phi_opt);
steering = exp(-1j*2*pi/lambda * antenna_pos * sin(phi_rad)).';
doppler_term = exp(-1j*2*pi*fd_opt * t_vec);
pulse_shifted = interp1(t_vec, p, t_vec - tau_opt, 'linear', 0);
s_opt = (steering) * (pulse_shifted .* doppler_term);
alpha_opt = sum(sum(conj(s_opt) .* x_hat)) / sum(sum(abs(s_opt).^2));

end
