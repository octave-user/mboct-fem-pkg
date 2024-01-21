## fem_tests.m:53
%!test
%! ## TEST 53
%! ## Robert Gasch, Klaus Knothe
%! ## Strukturdynamik Band 1
%! ## Diskrete Systeme
%! ## 1987
%! ## page 102, figure 2.5a
%!
%! a = 1800e-3;
%! b = 4000e-3;
%! c = 2200e-3;
%! m1 = 2.2e-1;
%! m2 = 4.1e-1;
%! J1 = 0.1e-3;
%! J2 = 0.5e-3;
%! d1 = 5e-3;
%! d2 = 8e-3;
%! d3 = 7e-3;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 0;
%! X = [0,         0,   0;
%!      a,         0,   0;
%!      a + b,     0,   0;
%!      a + b + c, 0,   0];
%!
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! mesh.material_data.rho = rho;
%! mesh.materials.beam2 = ones(3, 1, "int32");
%! mesh.elements.beam2(1).nodes = int32([1, 2]);
%! mesh.elements.beam2(2).nodes = int32([2, 3]);
%! mesh.elements.beam2(3).nodes = int32([3, 4]);
%! mesh.elements.beam2(1).section.A = d1^2 * pi / 4;
%! mesh.elements.beam2(1).section.Ay = 9/10 * d1^2 * pi / 4;
%! mesh.elements.beam2(1).section.Az = 9/10 * d1^2 * pi / 4;
%! mesh.elements.beam2(1).section.Iy = d1^4 * pi / 64;
%! mesh.elements.beam2(1).section.Iz = d1^4 * pi / 64;
%! mesh.elements.beam2(1).section.It = d1^4 * pi / 32;
%! mesh.elements.beam2(1).e2 = [0; 1; 0];
%! mesh.elements.beam2(2).section.A = d2^2 * pi / 4;
%! mesh.elements.beam2(2).section.Ay = 9/10 * d2^2 * pi / 4;
%! mesh.elements.beam2(2).section.Az = 9/10 * d2^2 * pi / 4;
%! mesh.elements.beam2(2).section.Iy = d2^4 * pi / 64;
%! mesh.elements.beam2(2).section.Iz = d2^4 * pi / 64;
%! mesh.elements.beam2(2).section.It = d2^4 * pi / 32;
%! mesh.elements.beam2(2).e2 = [0; 1; 0];
%! mesh.elements.beam2(3).section.A = d3^2 * pi / 4;
%! mesh.elements.beam2(3).section.Ay = 9/10 * d3^2 * pi / 4;
%! mesh.elements.beam2(3).section.Az = 9/10 * d3^2 * pi / 4;
%! mesh.elements.beam2(3).section.Iy = d3^4 * pi / 64;
%! mesh.elements.beam2(3).section.Iz = d3^4 * pi / 64;
%! mesh.elements.beam2(3).section.It = d3^4 * pi / 32;
%! mesh.elements.beam2(3).e2 = [0; 1; 0];
%! mesh.elements.bodies(1).nodes = int32(2);
%! mesh.elements.bodies(2).nodes = int32(3);
%! mesh.elements.bodies(1).m = m1;
%! mesh.elements.bodies(1).J = diag([0, J1, 0]);
%! mesh.elements.bodies(1).lcg = zeros(3, 1);
%! mesh.elements.bodies(2).m = m2;
%! mesh.elements.bodies(2).J = diag([0, J2, 0]);
%! mesh.elements.bodies(2).lcg = zeros(3, 1);
%! load_case_dof.locked_dof = false(size(mesh.nodes));
%! load_case_dof.locked_dof(:, [1,2,4,6]) = true;
%! load_case_dof.locked_dof(1, 1:3) = true;
%! load_case_dof.locked_dof(4, 2:3) = true;
%! gz = [-9.81, -1.625];
%! load_case = struct("g", cell(size(gz)));
%! for i=1:numel(load_case)
%!   load_case(i).g = [0; 0; gz(i)];
%! endfor
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.dm] = fem_ass_matrix(mesh, ...
%!                               dof_map, ...
%!                               [FEM_MAT_MASS, ...
%!                                FEM_MAT_STIFFNESS, ...
%!                                FEM_VEC_LOAD_CONSISTENT, ...
%!                                FEM_SCA_TOT_MASS], ...
%!                               load_case);
%! [sol] = fem_sol_static(mesh, dof_map, mat_ass);
%! Ia = d1^4 * pi / 64;
%! Ib = d2^4 * pi / 64;
%! Ic = d3^4 * pi / 64;
%! Mref = diag([m1, J1, m2, J2]);
%! Kref = [3 * E * Ia / a^3 + 12 * E * Ib / b^3, 3 * E * Ia / a^2 - 6 * E * Ib / b^2, -12 * E * Ib / b^3, -6 * E * Ib / b^2;
%!         3 * E * Ia / a^2 - 6 * E * Ib / b^2, 3 * E * Ia / a + 4 * E * Ib / b, 6 * E * Ib / b^2, 2 * E * Ib / b;
%!         -12 * E * Ib / b^3, 6 * E * Ib / b^2, 3 * E * Ic / c^3 + 12 * E * Ib / b^3, -3 * E * Ic / c^2 + 6 * E * Ib / b^2;
%!         -6 * E * Ib / b^2, 2 * E * Ib / b, -3 * E * Ic / c^2 + 6 * E * Ib / b^2, 3 * E * Ic / c + 4 * E * Ib / b];
%! Rref = [m1 * gz; zeros(size(gz)); m2 * gz; zeros(size(gz))];
%! Uref = Kref \ Rref;
%! U = zeros(size(Uref));
%! for i=1:numel(load_case)
%!  U(:, i) = sol.def(2:3, [3, 5], i).'(:);
%! endfor
%! tol = 1e-5;
%! assert_simple(U, Uref, tol * norm(Uref));
%! lambdaref = eig(Kref, Mref);
%! lambda = eig(mat_ass.K, mat_ass.M);
%! lambdaref = sort(lambdaref);
%! lambda = sort(lambda)(1:4);
%! tolf = 1e-5;
%! tolm = eps;
%! assert_simple(lambda, lambdaref, tolf * max(lambdaref));
%! assert_simple(mat_ass.dm, m1 + m2, tolm * (m1 + m2));
