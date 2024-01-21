## fem_pre_mesh_extrude_surf.m:01
%!test
%! ## TEST 1
%! rndstate = rand("state");
%! rand("seed", 0);
%! unwind_protect
%!   for i=1:50
%!     mesh.nodes = [0.0, 0.0, -1.0, 0.0, 0.0, 0.0;
%!                   1.0, 0.0, -1.0, 0.0, 0.0, 0.0;
%!                   0.0, 1.0, -1.0, 0.0, 0.0, 0.0;
%!                   0.5, 0.0, -1.0, 0.0, 0.0, 0.0;
%!                   0.5, 0.5, -1.0, 0.0, 0.0, 0.0;
%!                   0.0, 0.5, -1.0, 0.0, 0.0, 0.0];
%!     e1 = rand(3, 1);
%!     e2 = rand(3, 1);
%!     e3 = cross(e1, e2);
%!     e2 = cross(e3, e1);
%!     R = [e1, e2, e3];
%!     R /= diag(norm(R, "cols"));
%!     mesh.nodes(:, 1:3) *= R.';
%!     mesh.elements.tria6h = int32(1:6);
%!     mesh.groups.tria6h.id = int32(1);
%!     mesh.groups.tria6h.name = "surface";
%!     mesh.groups.tria6h.elements = int32(1);
%!     mesh.group.tria6h.nodes = unique(mesh.elements.tria6h);
%!     h = 2;
%!     N = 10;
%!     mesh.material_data.E = 210000e6;
%!     mesh.material_data.nu = 0.3;
%!     mesh.material_data.rho = 7850;
%!     [mesh.nodes, mesh.elements.penta15] = fem_pre_mesh_extrude_surf(mesh, "tria6h", int32(1), repmat(h / N, 1, N));
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!     dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!     [mat_ass.K, ...
%!      mat_ass.M, ...
%!      mat_ass.mtot, ...
%!      mat_ass.mat_info, ...
%!      mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_MAT_STIFFNESS, ...
%!                                           FEM_MAT_MASS, ...
%!                                           FEM_SCA_TOT_MASS]);
%!     N = 50;
%!     rho = 1e-2 * max(max(abs(mat_ass.K))) / max(max(abs(mat_ass.M)));
%!     sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, rho);
%!     mref = mesh.material_data.rho;
%!     rtol = sqrt(eps);
%!     assert_simple(mat_ass.mtot, mref, rtol * mref);
%!     assert_simple(all(sol_eig.f(1:6) < 1e-4 * max(sol_eig.f)));
%!     assert_simple(all(sol_eig.f(7:end) > 1e-1 * max(sol_eig.f)));
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", rndstate);
%! end_unwind_protect
