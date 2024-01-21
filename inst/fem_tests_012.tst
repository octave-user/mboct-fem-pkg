## fem_tests.m:12
%!test
%! ## TEST 12
%! tol = eps^0.9;
%! mesh.nodes = [2, 3, 4;
%!               6, 3, 2;
%!               2, 5, 1;
%!               4, 3, 6];
%! mesh.nodes = [mesh.nodes;
%!               0.5 * (mesh.nodes(1, :) + mesh.nodes(2, :));
%!               0.5 * (mesh.nodes(2, :) + mesh.nodes(3, :));
%!               0.5 * (mesh.nodes(1, :) + mesh.nodes(3, :));
%!               0.5 * (mesh.nodes(1, :) + mesh.nodes(4, :));
%!               0.5 * (mesh.nodes(2, :) + mesh.nodes(4, :));
%!               0.5 * (mesh.nodes(3, :) + mesh.nodes(4, :))];
%! mesh.nodes = [mesh.nodes, zeros(rows(mesh.nodes), 3)];
%! mesh.elements.tet10 = int32(1:10);
%! mesh.materials.tet10 = int32(1);
%! E = 480;
%! nu = 1/3;
%! mesh.material_data.rho = 1;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.loaded_nodes = zeros(1, 0, "int32");
%! load_case.loads = zeros(0, 3);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! iperm = dof_map.ndof(:, 1:3).';
%! iperm = reshape(iperm, numel(iperm), 1);
%! [K, ...
%!  M, ...
%!  Mlumped, ...
%!  dm] = fem_ass_matrix(mesh, ...
%!                       dof_map, ...
%!                       [FEM_MAT_STIFFNESS, ...
%!                        FEM_MAT_MASS, ...
%!                        FEM_MAT_MASS_LUMPED, ...
%!                        FEM_SCA_TOT_MASS]);
%! Vref = 4;
%! assert_simple(dm, Vref, tol * Vref);
%! assert_simple(isdefinite(K), false);
%! assert_simple(isdefinite(M), true);
%! assert_simple(isdefinite(Mlumped), true);
%! Kref = [447	324	72	1	-6	-12	54	48	0	94	66	36	-152	-90	12	55	42	-12	-311	-252	-24	-431	-306	-132	95	60	24	148	114	36
%!         324	1032	162	24	-104	-42	24	216	12	60	232	84	-180	-32	72	48	112	-30	-180	-992	-90	-288	-1040	-306	84	128	42	84	448	96
%!         72	162	339	0	-30	-35	0	24	54	24	60	94	-24	36	-8	0	-6	19	-24	-126	-275	-96	-234	-395	24	30	59	24	84	148
%!         1	24	0	87	-54	-36	18	-24	0	10	-18	-12	-32	-54	12	-83	90	12	19	0	0	11	6	-12	-59	72	48	28	-42	-12
%!         -6	-104	-30	-54	132	54	-12	72	12	0	76	36	36	268	72	54	-260	-54	-18	-32	-18	-6	-28	6	18	-272	-126	-12	148	48
%!         -12	-42	-35	-36	54	87	0	24	18	0	36	46	48	108	76	12	-90	-83	-12	-18	-17	-12	-6	11	12	-126	-167	0	60	64
%!         54	24	0	18	-12	0	108	0	0	-36	-12	0	72	12	0	-90	36	0	-198	-72	0	18	12	0	-18	-24	0	72	36	0
%!         48	216	24	-24	72	24	0	432	0	-24	-144	-48	24	288	48	72	-360	-72	-144	-792	-72	24	72	-24	-48	-72	-24	72	288	144
%!         0	12	54	0	12	18	0	0	108	0	-24	-36	0	24	72	0	-36	-90	0	-36	-198	0	-12	18	0	-12	-18	0	72	72
%!         94	60	24	10	0	0	-36	-24	0	204	108	72	104	60	24	-26	-24	0	58	36	24	-350	-216	-96	-98	-36	-24	40	36	-24
%!         66	232	60	-18	76	36	-12	-144	-24	108	492	216	48	308	96	-30	-68	12	54	88	36	-234	-860	-252	18	-392	-180	0	268	0
%!         36	84	94	-12	36	46	0	-48	-36	72	216	312	24	120	140	-12	-12	10	36	36	58	-132	-324	-386	12	-180	-242	-24	72	4
%!         -152	-180	-24	-32	36	48	72	24	0	104	48	24	1416	648	144	-392	-336	0	232	336	96	136	216	48	-680	-504	-240	-704	-288	-96
%!         -90	-32	36	-54	268	108	12	288	24	60	308	120	648	3936	864	-312	-1424	96	456	352	192	216	256	-288	-648	-1568	-576	-288	-2384	-576
%!         12	72	-8	12	72	76	0	48	72	24	96	140	144	864	1416	-48	-96	-248	144	96	232	48	-144	-152	-240	-432	-680	-96	-576	-848
%!         55	48	0	-83	54	12	-90	72	0	-26	-30	-12	-392	-312	-48	376	0	-96	-152	-192	0	-116	-72	48	292	144	0	136	288	96
%!         42	112	-6	90	-260	-90	36	-360	-36	-24	-68	-12	-336	-1424	-96	0	928	0	-96	256	96	-72	-176	72	216	736	216	144	256	-144
%!         -12	-30	19	12	-54	-83	0	-72	-90	0	12	10	0	96	-248	-96	0	376	96	192	136	48	72	-116	-48	72	148	0	-288	-152
%!         -311	-180	-24	19	-18	-12	-198	-144	0	58	54	36	232	456	144	-152	-96	96	1048	576	192	292	168	-48	-308	-144	-96	-680	-672	-288
%!         -252	-992	-126	0	-32	-18	-72	-792	-36	36	88	36	336	352	96	-192	256	192	576	2176	288	192	736	168	-144	-224	-72	-480	-1568	-528
%!         -24	-90	-275	0	-18	-17	0	-72	-198	24	36	58	96	192	232	0	96	136	192	288	760	0	120	148	-96	-72	-164	-192	-480	-680
%!         -431	-288	-96	11	-6	-12	18	24	0	-350	-234	-132	136	216	48	-116	-72	48	292	192	0	984	648	144	-152	-72	48	-392	-408	-48
%!         -306	-1040	-234	6	-28	-6	12	72	-12	-216	-860	-324	216	256	-144	-72	-176	72	168	736	120	648	2208	432	-216	256	144	-240	-1424	-48
%!         -132	-306	-395	-12	6	11	0	-24	18	-96	-252	-386	48	-288	-152	48	72	-116	-48	168	148	144	432	984	48	144	136	0	48	-248
%!         95	84	24	-59	18	12	-18	-48	0	-98	18	12	-680	-648	-240	292	216	-48	-308	-144	-96	-152	-216	48	696	216	144	232	504	144
%!         60	128	30	72	-272	-126	-24	-72	-12	-36	-392	-180	-504	-1568	-432	144	736	72	-144	-224	-72	-72	256	144	216	1056	432	288	352	144
%!         24	42	59	48	-126	-167	0	-24	-18	-24	-180	-242	-240	-576	-680	0	216	148	-96	-72	-164	48	144	136	144	432	696	96	144	232
%!         148	84	24	28	-12	0	72	72	0	40	0	-24	-704	-288	-96	136	144	0	-680	-480	-192	-392	-240	0	232	288	96	1120	432	192
%!         114	448	84	-42	148	60	36	288	72	36	268	72	-288	-2384	-576	288	256	-288	-672	-1568	-480	-408	-1424	48	504	352	144	432	3616	864
%!         36	96	148	-12	48	64	0	144	72	-24	0	4	-96	-576	-848	96	-144	-152	-288	-528	-680	-48	-48	-248	144	144	232	192	864	1408];
%! assert_simple(full(K(iperm,iperm)), Kref, tol * norm(Kref));
%! assert_simple(eig(K), eig(Kref), tol * max(eig(Kref)));
%! assert_simple(sum(full(diag(Mlumped))) / 3, dm, tol * dm);
