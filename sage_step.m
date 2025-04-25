%% Core Sage Algorithm
function [theta_list, l_updated] = sage_step(r, theta_list, mu, ...
    p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc)

% r: tín hiệu thu (MxN)
% theta_list: danh sách tham số hiện tại (1xL struct)
% mu: chỉ số bước hiện tại (iteration step)
% p: xung shaping pulse
% tau_grid, phi_grid, fd_grid: các giá trị quét
% T0: khoảng lấy mẫu
% antenna_pos: vị trí các anten
% fc: tần số sóng mang

L = length(theta_list);
l = mod(mu - 1, L) + 1;  % xác định đường truyền cần cập nhật

%% Bước E: loại bỏ ảnh hưởng các đường truyền khác
x_hat = compute_x_hat(r, l, theta_list, p, T0, antenna_pos, fc);

%% Bước M: tìm tham số tối ưu cho đường truyền l
[tau_opt, phi_opt, fd_opt, alpha_opt] = sage_mstep_grid(...
    x_hat, p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc);

%% Cập nhật tham số đường truyền l
theta_list(l).tau = tau_opt;
theta_list(l).phi = phi_opt;
theta_list(l).fd  = fd_opt;
theta_list(l).alpha = alpha_opt;

% Trả lại chỉ số đường truyền đã cập nhật
l_updated = l;
end
