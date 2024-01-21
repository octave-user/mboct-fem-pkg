## fem_tests.m:54
%!test
%! ## TEST 54
%! a = 0.2;
%! b = 0.4;
%! c = 0.3;
%! rho = 7850;
%! m = rho * a * b * c;
%! Jx = m * (b^2 + c^2) / 12;
%! Jy = m * (a^2 + c^2) / 12;
%! Jz = m * (a^2 + b^2) / 12;
%! m1_2 = rho * (a / 2) * b * c;
%! Jx1_2 = m1_2 * (b^2 + c^2) / 12;
%! Jy1_2 = m1_2 * ((a / 2)^2 + c^2) / 12;
%! Jz1_2 = m1_2 * ((a / 2)^2 + b^2) / 12;
%! J = diag([Jx, Jy, Jz]);
%! J1_2 = diag([Jx1_2, Jy1_2, Jz1_2]);
%! e1 = [1; 0.5; 0.3];
%! e2 = [0.2; 1; 0.2];
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R = [e1, e2, e3];
%! R *= diag(1 ./ norm(R, "cols"));
%! J = R * J * R.';
%! J = 0.5 * (J + J.');
%! J1_2 = R * J1_2 * R.';
%! J1_2 = 0.5 * (J1_2 + J1_2.');
%! lcg = R * [12e-2; 23e-2; 34e-2];
%! lcg1 = lcg + R * [0.25 * a; 0; 0];
%! lcg2 = lcg - R * [0.25 * a; 0; 0];
%! mesh.nodes = zeros(1, 6);
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! mesh.elements.bodies(1).nodes = int32(1);
%! mesh.elements.bodies(1).m = m / 2;
%! mesh.elements.bodies(1).J = J1_2;
%! mesh.elements.bodies(1).lcg = lcg1;
%! mesh.elements.bodies(2).nodes = int32(1);
%! mesh.elements.bodies(2).m = m / 2;
%! mesh.elements.bodies(2).J = J1_2;
%! mesh.elements.bodies(2).lcg = lcg2;
%! load_case.locked_dof = false(size(mesh.nodes));
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [M] = fem_ass_matrix(mesh, ...
%!                      dof_map, ...
%!                      [FEM_MAT_MASS, ...
%!                       FEM_MAT_MASS_SYM, ...
%!                       FEM_MAT_MASS_SYM_L], ...
%!                      load_case);
%! Mref = [      m * eye(3), -m * skew(lcg);
%!         -m * skew(lcg).', J - m * skew(lcg) * skew(lcg)];
%! tol = eps;
%! assert_simple(issymmetric(M, tol));
%! assert_simple(isdefinite(M, tol));
%! assert_simple(full(M), Mref, eps * norm(Mref));
