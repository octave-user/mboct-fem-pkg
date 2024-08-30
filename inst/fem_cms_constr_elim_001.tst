## fem_cms_constr_elim.m:01
%!test
%! try
%! ## Build a simple mesh made of a single hexahedron
%! X = [ 1,  1,  1;
%!      -1,  1,  1;
%!      -1, -1,  1;
%!       1, -1,  1;
%!       1,  1, -1;
%!      -1,  1, -1;
%!      -1, -1, -1;
%!       1, -1, -1];
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32(1:8);
%! mesh.materials.iso8 = int32(1);
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! ## Clamp all bottom nodes to the ground
%! load_case.locked_dof(5:8, :) = true;
%! ## Apply several loads at the top nodes
%! load_case.loaded_nodes = int32(1:4).';
%! load_case.loads = [1e10, 0, 0;
%!                    0, 2e10, 0;
%!                    0, 0, -3e10;
%!                    0, 1e10, 0];
%! ## Impose a constraint to the top nodes in a way,
%! ## that the displacement for all four nodes is identical in x, y and z-direction
%! mesh.elements.joints.nodes = int32(1:4);
%! mesh.elements.joints.C = [  eye(3, 6),  -eye(3, 6), zeros(3, 6), zeros(3, 6);
%!                           zeros(3, 6),   eye(3, 6),  -eye(3, 6), zeros(3, 6);
%!                           zeros(3, 6), zeros(3, 6),   eye(3, 6),  -eye(3, 6)];
%! ## Build the degree of freedom table   
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! ## Assemble global matrices
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_MASS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! ## Eliminate redundant equations
%! [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%! ## Solve the reduced set of equations
%! Ured = Tred * fem_sol_linsolve(Kred, Rred);
%! sol_red.def = zeros(size(mesh.nodes));
%! ## Compute nodal displacements
%! for i=1:columns(dof_map.ndof)
%!   idx = find(dof_map.ndof(:, i) > 0);
%!   sol_red.def(idx, i) = Ured(dof_map.ndof(idx, i));
%! endfor
%! tol = eps^0.8;
%! ## Solve the full set of equations
%! [sol_ref] = fem_sol_static(mesh, dof_map, mat_ass);
%! ## Check if we got the same result with reduced set of equations and full set of equations
%! assert_simple(sol_red.def, sol_ref.def, tol * norm(sol_ref.def));
%! ## Check if constraint we imposed for the top nodes is satisfied
%! assert_simple(sol_red.def(1:4, :), repmat(sol_red.def(1, :), 4, 1), tol * norm(sol_red.def(1, :)));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
