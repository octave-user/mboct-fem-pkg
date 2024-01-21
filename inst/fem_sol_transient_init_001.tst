## fem_sol_transient_init.m:01
%!test
%! enable_plot = true;
%! if (enable_plot)
%!   close all;
%! endif
%! for D=[0, 20, 40, 100]
%!   M = 10;
%!   K = 100;
%!   R = 0;
%!   U0 = 1;
%!   dU0_dt = -1;
%!   options.solver = "mldivide";
%!   options.delta = 0.5;
%!   lambda = -D / (2 * M) + [1, -1] * sqrt((D / (2 * M))^2 - K / M);
%!   if (all(isreal(lambda)))
%!     T = 4 / abs(lambda(1));
%!   else
%!     T = 4 * 2 * pi / imag(lambda(1));
%!   endif
%!   dt = T / 10000;
%!   t = 0:dt:T;
%!   U = dU_dt = dU_dt2 = zeros(1, numel(t));
%!   U(:, 1) = U0;
%!   dU_dt(:, 1) = dU0_dt;
%!   dU_dt2(:, 1) = M \ (R - K * U0 - D * dU0_dt);
%!   sol_dat = fem_sol_transient_init(M, D, K, dt, options);
%!   for i=2:columns(U)
%!     [U(:, i), dU_dt(:, i), dU_dt2(:, i)] = fem_sol_transient_step(U(:, i - 1), dU_dt(:, i - 1), dU_dt2(:, i - 1), R, sol_dat);
%!   endfor
%!   A1 = (dU0_dt - lambda(2) * U0) / (lambda(1) - lambda(2));
%!   A2 = U0 - A1;
%!   Uref = real(A1 * exp(lambda(1) * t) + A2 * exp(lambda(2) * t));
%!   dUref_dt = real(lambda(1) * A1 * exp(lambda(1) * t) + lambda(2) * A2 * exp(lambda(2) * t));
%!   dUref_dt2 = real(lambda(1)^2 * A1 * exp(lambda(1) * t) + lambda(2)^2 * A2 * exp(lambda(2) * t));
%!   if (enable_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(t, U(1, :), "-x;U(t);1");
%!     plot(t, Uref, "-;Uref(t);0");
%!     xlabel("t [s]");
%!     ylabel("U [m]");
%!     grid on;
%!     grid minor on;
%!     title("displacement versus time");
%!     figure("visible", "off");
%!     hold on;
%!     plot(t, dU_dt(1, :), "-x;dU(t)/dt;1");
%!     plot(t, dUref_dt, "-;dUref(t)/dt;0");
%!     xlabel("t [s]");
%!     ylabel("dU/dt [m/s]");
%!     grid on;
%!     grid minor on;
%!     title("velocity versus time");
%!     figure("visible", "off");
%!     hold on;
%!     plot(t, dU_dt2(1, :), "-x;d^2U(t)/dt^2;1");
%!     plot(t, dUref_dt2, "-;d^2Uref(t)/dt^2;0");
%!     xlabel("t [s]");
%!     ylabel("d^2U/dt^2 [m/s^2]");
%!     grid on;
%!     grid minor on;
%!     title("acceleration versus time");
%!   endif
%!   tol = eps^0.3;
%!   assert_simple(U, Uref, tol * max(abs(Uref)));
%!   assert_simple(dU_dt, dUref_dt, tol * max(abs(dUref_dt)));
%!   assert_simple(dU_dt2, dUref_dt2, tol * max(abs(dUref_dt2)));
%! endfor
