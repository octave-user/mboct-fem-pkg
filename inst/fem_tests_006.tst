## fem_tests.m:06
%!test
%! #############################################################
%! ## TEST 6: Test case for RBE3 element
%! ## model of two brick elements connected to one rbe3 element
%! #############################################################
%! close all;
%! b_plot = false;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! l = 80e-3;
%! beta = atan(c/b);
%! r = sqrt((c/2)^2 + (b/2)^2);
%! tol = eps^0.5;
%! scale_def = 10e-3;
%! Fx = 40000;
%! Fy = 50;
%! Fz = 100;
%! Mx = 50;
%! My = 30;
%! Mz = 40;
%! Phi9 = beta + pi / 2;
%! Phi10 = 3/2 * pi - beta;
%! Phi11 = pi / 2 - beta;
%! Phi12 = 3/2 * pi + beta;
%! Ft = Mx / (4 * r);
%! Fx9  = Fx / 4 - (Mz + Fy * l) / (2 * b) + (My - Fz * l) / (2 * c);
%! Fx10 = Fx / 4 + (Mz + Fy * l) / (2 * b) + (My - Fz * l) / (2 * c);
%! Fx11 = Fx / 4 - (Mz + Fy * l) / (2 * b) + (Fz * l - My) / (2 * c);
%! Fx12 = Fx / 4 + (Mz + Fy * l) / (2 * b) + (Fz * l - My) / (2 * c);
%! Fy9 = Fy / 4 + Ft * cos(Phi9);
%! Fy10 = Fy / 4 + Ft * cos(Phi10);
%! Fy11 = Fy / 4 + Ft * cos(Phi11);
%! Fy12 = Fy / 4 + Ft * cos(Phi12);
%! Fz9 = Fz / 4 + Ft * sin(Phi9);
%! Fz10 = Fz / 4 + Ft * sin(Phi10);
%! Fz11 = Fz / 4 + Ft * sin(Phi11);
%! Fz12 = Fz / 4 + Ft * sin(Phi12);
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  # 1
%!       0,  0.5 * b,  0.5 * c;  # 2
%!       0, -0.5 * b,  0.5 * c;  # 3
%!       0.5 * a, -0.5 * b,  0.5 * c;  # 4
%!       0.5 * a,  0.5 * b, -0.5 * c;  # 5
%!       0,  0.5 * b, -0.5 * c;  # 6
%!       0, -0.5 * b, -0.5 * c;  # 7
%!       0.5 * a, -0.5 * b, -0.5 * c,  # 8
%!       a,  0.5 * b,  0.5 * c;  # 9
%!       a, -0.5 * b,  0.5 * c;  # 10
%!       a,  0.5 * b, -0.5 * c;  # 11
%!       a, -0.5 * b, -0.5 * c,  # 12
%!       a + l,        0,        0]; # 13
%! Phi1 = [0, 30, 60] * pi / 180;
%! Phi2 = [0, 20, 90] * pi / 180;
%! Phi3 = [0, 70, 10] * pi / 180;
%! for j=1:length(Phi1)
%!   Rref = euler123_to_rotation_matrix([Phi1(j); Phi2(j); Phi3(j)]);
%!   Tref = [Rref, zeros(3, 3);
%!           zeros(3, 3), Rref];
%!   if (b_plot)
%!     figure("visible", "off");
%!   endif
%!   for i=1:2
%!     if (i==2)
%!       idx_node = 1:13;
%!     else
%!       idx_node = 1:12;
%!     endif
%!     data(i, j).mesh.nodes = [X(idx_node, :), zeros(length(idx_node), 3)] * Tref.';
%!     data(i, j).mesh.elements.iso8 = int32([1:8;
%!                                            9, 1, 4, 10, 11, 5, 8, 12]);
%!     data(i, j).mesh.materials.iso8 = int32([1; 1]);
%!     if (i==2)
%!       data(i, j).mesh.elements.rbe3.nodes = int32([13, 9, 10, 11, 12]);
%!       data(i, j).mesh.elements.rbe3.weight = ones(1, 4);
%!     endif
%!     E = 210000e6;
%!     nu = 0.3;
%!     data(i, j).mesh.material_data.rho = 7850;
%!     data(i, j).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!     data(i, j).load_case.locked_dof = false(rows(data(i, j).mesh.nodes), 6);
%!     data(i, j).load_case.locked_dof([2, 3, 6, 7], 1:6) = true;
%!     if i==2
%!       data(i, j).load_case.loaded_nodes = int32([13]);
%!       data(i, j).load_case.loads = [Fx, Fy, Fz, Mx, My, Mz] * Tref.';
%!     else
%!       data(i, j).load_case.loaded_nodes = int32([9; 10; 11; 12]);
%!       data(i, j).load_case.loads = [Fx9,  Fy9,  Fz9,  0, 0, 0;
%!                                     Fx10, Fy10, Fz10, 0, 0, 0;
%!                                     Fx11, Fy11, Fz11, 0, 0, 0;
%!                                     Fx12, Fy12, Fz12, 0, 0, 0] * Tref.';
%!     endif
%!     data(i, j).dof_map = fem_ass_dof_map(data(i, j).mesh, data(i, j).load_case);
%!     [data(i, j).mat_ass.K, ...
%!      data(i, j).mat_ass.R] = fem_ass_matrix(data(i, j).mesh, ...
%!                                             data(i, j).dof_map, ...
%!                                             [FEM_MAT_STIFFNESS, ...
%!                                              FEM_VEC_LOAD_CONSISTENT], ...
%!                                             data(i, j).load_case);
%!     assert_simple(rank(data(i, j).mat_ass.K), columns(data(i, j).mat_ass.K));
%!     data(i, j).U = full(data(i, j).mat_ass.K \ data(i, j).mat_ass.R);
%!     data(i, j).sol_stat.def = fem_post_def_nodal(data(i, j).mesh, data(i, j).dof_map, data(i, j).U);
%!     if (b_plot)
%!       fem_post_sol_plot(data(i, j).mesh);
%!       fem_post_sol_plot(data(i, j).mesh, data(i, j).sol_stat, scale_def / max(norm(data(1,1).sol_stat.def(1:12, 1:3), "rows")));
%!     endif
%!   endfor
%!   assert_simple(data(2, j).sol_stat.def(1:12, :) * Tref, data(1, 1).sol_stat.def, tol * max(norm(data(1,1).sol_stat.def, "rows")));
%!   if (b_plot)
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!   endif
%! endfor
%! return;
%! ndof = data(2,1).dof_map.ndof;
%! edof = data(2,1).dof_map.edof.rbe3;
%! K = full(data(2,1).mat_ass.K);
%! sdof = 1:ndof(12, 3);
%! mdof = ndof(13, :);
%! U1 = data(2,1).U(sdof);
%! K11 = K(sdof, sdof);
%! R1 = K11 * U1;
%! UM = data(2,1).U(mdof);
%! RM = data(2,1).load_case.loads.';
%! B1 = K(sdof, edof(1:6));
%! B2 = K(sdof, edof(7:12));
%! B3 = K(edof(1:6), sdof).';
%! B4 = K(edof(7:12), sdof).';
%! xi = data(2,1).mesh.nodes(data(2,1).mesh.elements.rbe3.nodes, 1:3).';
%! xi = xi(:, 2:end) - xi(:, 1);
%! S = zeros(columns(xi) * 6, 6);
%! for i=1:columns(xi)
%!   S((i - 1) * 6 + (1:6), :) = eye(6);
%!   S((i - 1) * 6 + (1:3), 4:6) = -skew(xi(:, i));
%! endfor
%! Lc = max(norm(xi, "cols"));
%! W = ones(1, columns(xi)*6);
%! W([4:6:end,5:6:end,6:6:end]) = Lc^2;
%! W = diag(W);
%! X = inv(S.' * W * S);
%! B = W * S * X;
%! lambda1 = data(2,1).U(edof(1:6));
%! lambda2 = data(2,1).U(edof(7:12));
%! R = B1 * lambda1 + B2 * lambda2;
%! assert_simple(-lambda1 - lambda2, RM);
%! assert_simple(R, -R1, tol * norm(R1));
%! assert_simple(B3.' * U1, UM, tol * norm(UM));
%! assert_simple(B4.' * U1, UM, tol * norm(UM));
%! assert_simple(B1, B2);
%! assert_simple(B1, B3);
%! assert_simple(B1, B4);
