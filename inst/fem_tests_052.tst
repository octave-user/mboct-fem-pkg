## fem_tests.m:52
%!test
%! ## TEST 52
%! m = 1.5;
%! J = diag([1e-3, 2e-3, 3e-3]);
%! e1 = [1; 0.5; 0.3];
%! e2 = [0.2; 1; 0.2];
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R = [e1, e2, e3];
%! R *= diag(1 ./ norm(R, "cols"));
%! J = R * J * R.';
%! J = 0.5 * (J + J.');
%! lcg = [2e-2; 3e-2; 4e-2];
%! mesh.nodes = zeros(1, 6);
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! mesh.elements.bodies.nodes = int32(1);
%! mesh.elements.bodies.m = m;
%! mesh.elements.bodies.J = J;
%! mesh.elements.bodies.lcg = lcg;
%! load_case.locked_dof = false(size(mesh.nodes));
%! load_case.g = [5; 4; -3];
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [M, MU, ML, R, dm] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_MASS, ...
%!                                      FEM_MAT_MASS_SYM, ...
%!                                      FEM_MAT_MASS_SYM_L, ...
%!                                      FEM_VEC_LOAD_CONSISTENT, ...
%!                                      FEM_SCA_TOT_MASS], ...
%!                                     load_case);
%! Mref = [      m * eye(3), -m * skew(lcg);
%!         -m * skew(lcg).', J - m * skew(lcg) * skew(lcg)];
%! tol = eps;
%! assert_simple(full(R), [m * load_case.g; cross(lcg, load_case.g) * m], tol * norm(m * load_case.g));
%! assert_simple(issymmetric(M, tol));
%! assert_simple(isdefinite(M, tol));
%! assert_simple(full(M), Mref, tol * norm(Mref));
%! assert_simple(full(MU + ML - diag(diag(ML))), Mref, tol * norm(Mref));
%! assert_simple(full(MU + ML - diag(diag(MU))), Mref, tol * norm(Mref));
%! assert_simple(full(ML - MU.'), zeros(6, 6), tol * norm(Mref));
%! assert_simple(dm, m, tol * m);
