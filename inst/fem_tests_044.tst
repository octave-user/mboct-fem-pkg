## fem_tests.m:44
%!test
%! ## TEST 44
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   a = 50e-3;
%!   b = 20e-3;
%!   c = 15e-3;
%!   x = a * [ 0,  1,  0, 0, 1, 0, 1/2, 1/2,   0, 1/2, 1/2,   0, 0, 1, 0 ];
%!   y = b * [ 0,  0,  1, 0, 0, 1,   0, 1/2, 1/2,   0, 1/2, 1/2, 0, 0, 1 ];
%!   z = 0.5 * c * [-1, -1, -1, 1, 1, 1,  -1,  -1,  -1,   1,   1,   1, 0, 0, 0 ];
%!   for j=1:100
%!     e1 = rand(3, 1);
%!     e2 = rand(3, 1);
%!     e3 = cross(e1, e2);
%!     e2 = cross(e3, e1);
%!     e1 /= norm(e1);
%!     e2 /= norm(e2);
%!     e3 /= norm(e3);
%!     R = [e1, e2, e3];
%!     mesh.nodes = [(R * [x; y; z]).', zeros(numel(x), 3)];
%!     mesh.elements.penta15 = int32(1:rows(mesh.nodes));
%!     mesh.material_data.E = 210000e6;
%!     mesh.material_data.nu = 0.3;
%!     mesh.material_data.rho = 7850;
%!     mesh.materials.penta15 = int32(1);
%!     load_case.locked_dof = false(size(mesh.nodes));
%!     dof_map = fem_ass_dof_map(mesh, load_case);
%!     [mat_ass.K, ...
%!      mat_ass.M] = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS, ...
%!                                   FEM_MAT_MASS, ...
%!                                   FEM_MAT_MASS_LUMPED, ...
%!                                   FEM_VEC_LOAD_CONSISTENT, ...
%!                                   FEM_SCA_TOT_MASS], ...
%!                                   load_case);
%!     assert_simple(isdefinite(mat_ass.K) == 0);
%!     assert_simple(isdefinite(mat_ass.M) == 1);
%!     assert_simple(rank(mat_ass.K), columns(mat_ass.K) - 6);
%!     tol_eig = sqrt(eps);
%!     N = 20;
%!     shift = 1e-2 * max(max(abs(mat_ass.K))) / max(max(abs(mat_ass.M)));
%!     sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, shift);
%!     sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     FEM_SCA_STRESS_VMIS, ...
%!                                     load_case, ...
%!                                     sol_eig);
%!     assert_simple(max(sol_eig.f(1:6)) < eps^0.3 * min(sol_eig.f(7:end)));
%!     assert_simple(max(max(abs(sol_eig.stress.vmis.penta15(:,:,1:6)))) < eps^0.5 * min(min(abs(sol_eig.stress.vmis.penta15(:, :, 7:end)))));
%!     sol_eig_t = sol_eig;
%!     for i=1:size(sol_eig.def, 3)
%!       sol_eig_t.def(:, 1:3, i) = (R.' * sol_eig_t.def(:, 1:3, i).').';
%!       sol_eig_t.def(:, :, i) /= max(max(abs(sol_eig_t.def(:, :, i))));
%!     endfor
%!     if (j == 1)
%!       U_ref_eig = sol_eig_t.def;
%!       f_ref = sol_eig.f;
%!       vmis_eig = sol_eig.stress.vmis.penta15;
%!     else
%!       assert_simple(sol_eig.f(7:end), f_ref(7:end), sqrt(eps) * max(f_ref));
%!       assert_simple(sol_eig.stress.vmis.penta15, vmis_eig, tol_eig * max(max(max(abs(vmis_eig)))));
%!       for i=7:size(sol_eig.def, 3)
%!         try
%!           assert_simple(sol_eig_t.def(:, :, i), U_ref_eig(:, :, i), tol_eig);
%!         catch
%!           assert_simple(sol_eig_t.def(:, :, i), -U_ref_eig(:, :, i), tol_eig);
%!         end_try_catch
%!       endfor
%!     endif
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
