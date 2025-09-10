## fem_cms_create.m:01
%!test
%! try
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1e-1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! a = 150e-3 / SI_unit_m;
%! b = 20e-3 / SI_unit_m;
%! c = 45e-3 / SI_unit_m;
%! d = 10e-3 / SI_unit_m;
%! e = 10e-3 / SI_unit_m;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  ##  1
%!       0,  0.5 * b,  0.5 * c;  ##  2
%!       0, -0.5 * b,  0.5 * c;  ##  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ##  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ##  5
%!       0,  0.5 * b, -0.5 * c;  ##  6
%!       0, -0.5 * b, -0.5 * c;  ##  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ##  8
%!       a,  0.5 * b,  0.5 * c;  ##  9
%!       a, -0.5 * b,  0.5 * c;  ## 10
%!       a,  0.5 * b, -0.5 * c;  ## 11
%!       a, -0.5 * b, -0.5 * c,  ## 12
%!       a + d,        0,        0;  ## 13
%!       -e,        0,        0]; ## 14
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8; 9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3(1).weight = ones(1, 4);
%! mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%! mesh.elements.rbe3(2).weight = ones(1, 4);
%! E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%! nu = 0.3;
%! mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! cms_opt.nodes.modal.number = int32(14);
%! cms_opt.nodes.interfaces.number = int32(13);
%! cms_opt.tol = 1e-3;
%! sol = {"pardiso", "mumps", "umfpack", "chol", "lu", "mldivide"};
%! alg = {"shift-invert", "diag-shift-invert", "unsymmetric", "eliminate"};
%! scaling = {"none", "max K", "max M", "max K,M", "norm K", "norm M", "norm K,M", "diag K", "diag M", "lambda", "Tred", "mean M,K", "mean K,M"};
%! use_static_modes = [true, false];
%! tol = 1e-6;
%! for stat_modes=use_static_modes
%!   cms_opt.static_modes = stat_modes;
%!   for iter=[10];
%!     for modes=int32([0, 4, 8, 10])
%!       lambda_ref = [];
%!       Phi_ref = [];
%!       for iscal=1:numel(scaling)
%!      cms_opt.scaling=scaling{iscal};
%!      for isol=1:numel(sol)
%!        cms_opt.solver = sol{isol};
%!        for ialg=1:numel(alg)
%!          for invariants=[true, false]
%!            for verbose=[false]
%!              for threads=int32([1, 2])
%!                cms_opt.verbose = verbose;
%!                cms_opt.modes.number = modes;
%!                cms_opt.number_of_threads = threads;
%!                cms_opt.threshold_elem = int32(1);
%!                cms_opt.algorithm = alg{ialg};
%!                cms_opt.invariants = invariants;
%!                cms_opt.refine_max_iter = iter;
%!                cms_opt.pre_scaling = true;
%!                [mesh_cms, ...
%!                 mat_ass_cms, ...
%!                 dof_map_cms, ...
%!                 sol_eig_cms, ...
%!                 cms_opt] = fem_cms_create(mesh, load_case, cms_opt);
%!                [Phi, lambda] = eig(mat_ass_cms.Kred, mat_ass_cms.Mred);
%!                [lambda, idx] = sort(diag(lambda));
%!                Phi = mat_ass_cms.Tred * Phi(:, idx);
%!                Phi *= diag(1 ./ max(abs(Phi), [], 1));
%!                if (numel(lambda_ref))
%!                  assert_simple(lambda, lambda_ref, tol * max(abs(lambda)));
%!                  for j=1:columns(Phi)
%!                    f = min(max(abs(Phi(:, j) + Phi_ref(:, j))), max(abs(Phi(:, j) - Phi_ref(:, j))));
%!                    assert_simple(f < tol);
%!                  endfor
%!                else
%!                  lambda_ref = lambda;
%!                  Phi_ref = Phi;
%!                endif
%!              endfor
%!            endfor
%!          endfor
%!        endfor
%!      endfor
%!       endfor
%!     endfor
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
