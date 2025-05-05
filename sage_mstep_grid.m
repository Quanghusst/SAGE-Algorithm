function [tau_opt, phi_opt, fd_opt, alpha_opt] = sage_mstep_grid(...
    x_hat, p, tau_grid, phi_grid, fd_grid, T0, antenna_pos, fc)

  c      = 3e8;
  lambda = c/fc;
  [M, N] = size(x_hat);
  t_vec  = (0:N-1) * T0;

  % kiểm tra kích thước p
  assert(length(p) == N, 'length(p) phải = N.');

  max_z = -inf;
  for tau = tau_grid
    pulse = interp1(t_vec, p, t_vec - tau, 'linear', 0);
    for phi = phi_grid
      phase    = antenna_pos.' * [cos(phi); sin(phi)];      % M×1
      steering = exp(1j*2*pi/lambda * phase);
      for fd = fd_grid
        doppler = exp(-1j*2*pi*fd * t_vec);                 % 1×N
        s       = steering * (pulse .* doppler);           % M×N
        z       = abs(sum(conj(s).*x_hat, 'all'));
        if z > max_z
          max_z   = z;
          tau_opt = tau;
          phi_opt = phi;
          fd_opt  = fd;
        end
      end
    end
  end

  % Tính alpha
  phase_opt    = antenna_pos.' * [cos(phi_opt); sin(phi_opt)];
  steering_opt = exp(1j*2*pi/lambda * phase_opt);
  doppler_opt  = exp(-1j*2*pi*fd_opt * t_vec);
  pulse_opt    = interp1(t_vec, p, t_vec - tau_opt, 'linear', 0);
  s_opt        = steering_opt * (pulse_opt .* doppler_opt);
  alpha_opt    = sum(conj(s_opt).*x_hat, 'all') / sum(abs(s_opt).^2, 'all');
end
