%!test
%! try
%! e1 = [0.3; 0.5; 0.2];
%! e2 = [0.7; 0.8; 0.9];
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R = [e1, e2, e3];
%! R *= diag(1 ./ norm(R, "cols"));
%! a = 0.75;
%! b = 0.45;
%! X = [ 0.5 * a,   0.5 * b,      0;
%!      -0.5 * a,   0.5 * b,      0;
%!      -0.5 * a,  -0.5 * b,      0;
%!       0.5 * a,  -0.5 * b,      0;
%!      -0.5 * a,  -0.5 * b - b,  0;
%!       0.5 * a,  -0.5 * b - b,  0];
%! mesh.nodes = [X * R.', zeros(rows(X), 3)];
%! sigma = 1;
%! rho = 1.2041;
%! c = 343.25;
%! omega = 2 * pi * [100, 2000, 10000];
%! Uz = 0.3e-6;
%! dA = 2 * a * b;
%! for i=1:rows(mesh.nodes)
%!   mesh.elements.joints(i).nodes = int32(i);
%!   mesh.elements.joints(i).C = R.' * [eye(3), zeros(3, 3)];
%!   for j=1:numel(omega)
%!     load_case(j).joints(i).U = [0; 0; Uz];
%!   endfor
%! endfor
%! mesh.elements.structural_boundary.iso4 = int32([1, 2, 3, 4;
%!                                                 4, 3, 5, 6]);
%! mesh.materials.structural_boundary.iso4 = ones(2, 1, "int32");
%! mesh.material_data(1).c = c;
%! mesh.material_data(1).rho = rho;
%! load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%! load_case_dof.dof_in_use = false(rows(mesh.nodes), 7);
%! load_case_dof.dof_in_use(:, 1:3) = true;
%! load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                       FEM_VEC_LOAD_FLUID_STRUCT], ...
%!                                      load_case);
%! opt_sol.solver = "umfpack";
%! mat_ass.R = complex(mat_ass.R);
%! U = fem_sol_factor(mat_ass.K, opt_sol) \ mat_ass.R;
%! [sol.def, sol.Phi] = fem_post_def_nodal(mesh, dof_map, U);
%! sol.def = complex(sol.def);
%! sol.omega = omega;
%! [sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                           dof_map, ...
%!                                           [FEM_SCA_EFFECTIVE_RADIATED_POWER_C], ...
%!                                           load_case, ...
%!                                           sol);
%! Pref = 0.5 * sigma * rho * c * abs(1j * omega * Uz).^2 * dA;
%! tol = eps^0.9;
%! assert_simple(sum(sol.acoustic_intensity.P.iso4, 1), Pref, tol * norm(Pref));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%! e1 = [0.3; 0.5; 0.2];
%! e2 = [0.7; 0.8; 0.9];
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R = [e1, e2, e3];
%! R *= diag(1 ./ norm(R, "cols"));
%! a = 0.75;
%! b = 0.45;
%! X = [ -0.5 * a,   0.5 * b,      0;
%!              0,   0.5 * b,      0;
%!        0.5 * a,   0.5 * b,      0;
%!       -0.5 * a,         0,      0;
%!              0,         0,      0;
%!        0.5 * a,         0,      0;
%!       -0.5 * a,  -0.5 * b,      0;
%!              0,  -0.5 * b,      0;
%!        0.5 * a,  -0.5 * b,      0];
%! mesh.nodes = [X * R.', zeros(rows(X), 3)];
%! sigma = 1;
%! rho = 1.2041;
%! c = 343.25;
%! omega = 2 * pi * [100, 2000, 10000];
%! Uz = 0.3e-6;
%! dA = a * b;
%! for i=1:rows(mesh.nodes)
%!   mesh.elements.joints(i).nodes = int32(i);
%!   mesh.elements.joints(i).C = R.' * [eye(3), zeros(3, 3)];
%!   for j=1:numel(omega)
%!     load_case(j).joints(i).U = [0; 0; Uz];
%!   endfor
%! endfor
%! mesh.elements.structural_boundary.tria6h = int32([1, 7, 9, 4, 8, 5;
%!                                                   3, 1, 9, 2, 5, 6]);
%! mesh.materials.structural_boundary.tria6h = ones(2, 1, "int32");
%! mesh.elements.structural_boundary.tria6 = int32([1, 7, 9, 4, 8, 5;
%!                                                  3, 1, 9, 2, 5, 6]);
%! mesh.materials.structural_boundary.tria6 = ones(2, 1, "int32");
%! mesh.material_data(1).c = c;
%! mesh.material_data(1).rho = rho;
%! load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%! load_case_dof.dof_in_use = false(rows(mesh.nodes), 7);
%! load_case_dof.dof_in_use(:, 1:3) = true;
%! load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                       FEM_VEC_LOAD_FLUID_STRUCT], ...
%!                                      load_case);
%! opt_sol.solver = "umfpack";
%! mat_ass.R = complex(mat_ass.R);
%! U = fem_sol_factor(mat_ass.K, opt_sol) \ mat_ass.R;
%! [sol.def, sol.Phi] = fem_post_def_nodal(mesh, dof_map, U);
%! sol.def = complex(sol.def);
%! sol.omega = omega;
%! [sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                           dof_map, ...
%!                                           [FEM_SCA_EFFECTIVE_RADIATED_POWER_C], ...
%!                                           load_case, ...
%!                                           sol);
%! Pref = 0.5 * sigma * rho * c * abs(1j * omega * Uz).^2 * dA;
%! tol = eps^0.9;
%! assert_simple(sum(sol.acoustic_intensity.P.tria6h, 1), Pref, tol * norm(Pref));
%! assert_simple(sum(sol.acoustic_intensity.P.tria6, 1), Pref, tol * norm(Pref));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
