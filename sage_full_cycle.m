%% Core EM Algorithm
function theta_list = sage_full_cycle(r, theta_list, ...
    p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc)

% r: tín hiệu thu (MxN)
% theta_list: danh sách các tham số hiện tại của L đường truyền
% p: shaping pulse
% tau_grid, phi_grid, fd_grid: các giá trị quét
% T0: thời gian lấy mẫu
% antenna_pos: vị trí anten
% fc: tần số sóng mang

L = length(theta_list);  % Số đường truyền

for l = 1:L
    % E-step: tính x̂_l (tín hiệu riêng phần còn lại)
    x_hat = compute_x_hat(r, l, theta_list, p, T0, antenna_pos, fc);

    % M-step: tìm tham số tối ưu cho đường truyền l
    [tau_opt, phi_opt, fd_opt, alpha_opt] = sage_mstep_grid(...
        x_hat, p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc);

    % Cập nhật lại vào danh sách theta_list
    theta_list(l).tau = tau_opt;
    theta_list(l).phi = phi_opt;
    theta_list(l).fd  = fd_opt;
    theta_list(l).alpha = alpha_opt;
end

end
