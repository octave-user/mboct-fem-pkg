## fem_tests.m:43
%!test
%! try
%! ## TEST 43
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
%!     load_case.locked_dof([1, 3, 4, 6, 9, 12, 13, 15], 1:3) = true;
%!     load_case.loads = [[-1 * e3.'; -2 * e3.'; -1 * e3.'], zeros(3, 3)];
%!     load_case.loaded_nodes = int32([2; 14; 5]);
%!     dof_map = fem_ass_dof_map(mesh, load_case);
%!     [mat_ass.K, ...
%!      mat_ass.M, ...
%!      mat_ass.Mdiag, ...
%!      mat_ass.R, ...
%!      mat_ass.mtot] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_STIFFNESS, ...
%!                                      FEM_MAT_MASS, ...
%!                                      FEM_MAT_MASS_LUMPED, ...
%!                                      FEM_VEC_LOAD_CONSISTENT, ...
%!                                      FEM_SCA_TOT_MASS], ...
%!                                     load_case);
%!     assert_simple(isdefinite(mat_ass.K));
%!     assert_simple(isdefinite(mat_ass.M));
%!     assert_simple(isdefinite(mat_ass.Mdiag));
%!     m_a = 0.5 * a * b * c * mesh.material_data.rho;
%!     tol_m = eps^0.8;
%!     tol_stat = sqrt(eps);
%!     tol_eig = sqrt(eps);
%!     assert_simple(mat_ass.mtot, m_a, tol_m * m_a);
%!
%!     sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!     sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      FEM_SCA_STRESS_VMIS, ...
%!                                      load_case, ...
%!                                      sol_stat);
%!     sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, 3);
%!     sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     FEM_SCA_STRESS_VMIS, ...
%!                                     load_case, ...
%!                                     sol_eig);
%!     for i=1:size(sol_eig.def, 3)
%!       sol_eig.def(:, 1:3, i) = (R.' * sol_eig.def(:, 1:3, i).').';
%!       sol_eig.def(:, :, i) /= max(max(abs(sol_eig.def(:, :, i))));
%!     endfor
%!     if (j == 1)
%!       U_ref_stat = R.' * sol_stat.def(:, 1:3).';
%!       U_ref_eig = sol_eig.def;
%!       f_ref = sol_eig.f;
%!       vmis_stat = sol_stat.stress.vmis.penta15;
%!       vmis_eig = sol_eig.stress.vmis.penta15;
%!     else
%!       assert_simple(R.' * sol_stat.def(:, 1:3).', U_ref_stat, tol_stat * norm(U_ref_stat));
%!       assert_simple(sol_eig.f, f_ref, sqrt(eps) * max(f_ref));
%!       assert_simple(sol_stat.stress.vmis.penta15, vmis_stat, tol_stat * max(max(max(abs(vmis_stat)))));
%!       assert_simple(sol_eig.stress.vmis.penta15, vmis_eig, tol_eig * max(max(max(abs(vmis_eig)))));
%!       for i=1:size(sol_eig.def, 3)
%!         Phi1 = sol_eig.def(:, :, i)(:);
%!         Phi2 = U_ref_eig(:, :, i)(:);
%!         assert_simple(min([norm(Phi1 - Phi2), norm(Phi1 + Phi2)]) < tol_eig);
%!       endfor
%!     endif
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
