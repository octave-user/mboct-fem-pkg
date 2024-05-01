## fem_pre_mesh_import.m:104
%!test
%! ### TEST 104
%! ### Code_Aster TLL100 V4.21.100 12/12/2011
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     dx = 2e-3;
%!     L = 100e-3;
%!     W = dx;
%!     H = dx;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L=%g;\n", L);
%!     fprintf(fd, "W=%g;\n", W);
%!     fprintf(fd, "H=%g;\n", H);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {-L,0,0,dx};\n");
%!     fputs(fd, "Point(2) = {-L,W,0,dx};\n");
%!     fputs(fd, "Point(3) = {-L,W,H,dx};\n");
%!     fputs(fd, "Point(4) = {-L,0,H,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {2 * L,0,0} {\n");
%!     fputs(fd, "  Surface{6};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"convection\",1) = {6};\n");
%!     fputs(fd, "Physical Surface(\"source\",2) = {tmp[0]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"tria6", "tet10"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   lambda = 1;
%!   rho = 1000;
%!   cp = 1;
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh.material_data.cp = cp;
%!   thetae = 100;
%!   theta0 = 100;
%!   h = 100;
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_THERMAL;
%!   mesh.elements.convection.tria6.nodes = mesh.elements.tria6([mesh.groups.tria6.elements], :);
%!   mesh.elements.convection.tria6.h = repmat(h, size(mesh.elements.convection.tria6.nodes));
%!   load_case.convection.tria6.theta = repmat(thetae, size(mesh.elements.convection.tria6.nodes));
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   e1 = [1; 0.5; 0.3];
%!   e2 = [0.1; 1; 0.2];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   e1 /= norm(e1);
%!   e2 /= norm(e2);
%!   e3 /= norm(e3);
%!   R = [e1, e2, e3];
%!   [x, idx_x] = sort(mesh.nodes(:, 1));
%!   mesh.nodes = [mesh.nodes(:, 1:3) * R.', mesh.nodes(:, 4:6) * R.'];
%!   [mat_ass.C, ...
%!    mat_ass.Kk, ...
%!    mat_ass.Qc] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_HEAT_CAPACITY, ...
%!                                  FEM_MAT_THERMAL_COND, ...
%!                                  FEM_VEC_LOAD_THERMAL], ...
%!                                 load_case);
%!   dt = rho * cp * dx^2 / lambda;
%!   alpha = 0.5;
%!   sol.t = 0:dt:5;
%!   sol.theta = zeros(dof_map.totdof, numel(sol.t));
%!   sol.theta(:, 1) = theta0;
%!   A = (1 / dt) * mat_ass.C + alpha * mat_ass.Kk;
%!   opts.number_of_threads = mbdyn_solver_num_threads_default();
%!   Afact = fem_sol_factor(A, opts);
%!   f = interp1([0, dt, sol.t(end)],[1, 0, 0], sol.t, "linear");
%!   for i=2:numel(sol.t)
%!     Qci = mat_ass.Qc * (f(i) * alpha + f(i - 1) * (1 - alpha));
%!     sol.theta(:, i) = Afact \ (mat_ass.C * (sol.theta(:, i - 1)) / dt - mat_ass.Kk * (sol.theta(:, i - 1) * (1 - alpha)) + Qci);
%!   endfor
%!   zetai = [1.428870011214117, 4.305801413119235, 7.228109771627203,10.2002625882959, 13.21418568384292, 16.25936122550391, 19.3270342916027,22.41084832844337, 25.50638298897743, 28.61058193655589, 31.72131067126047,34.83705394555897, 37.95671619253474, 41.07949058633233, 44.2047724244735,47.33210101880358, 50.46112006932919, 53.59155018094449, 56.72316948239282,59.85579973950923, 62.98929625447501, 66.12354041468544, 69.2584341235045,72.39389558637052, 75.52985608591556, 78.66625748768846, 81.80305029184685,84.94019209721015, 88.0776463800104, 91.21538151498968, 94.3533699848317,97.49158773721112, 100.6300136583878, 103.7686291395413, 106.9074177174041,110.0463647748095, 113.1854572898642, 116.3246836248031, 119.4640333473853,122.6034970790071, 125.7430663658736, 128.8827335673996, 132.0224917613492,135.1623346613892, 138.3022565459636, 141.4422521963534, 144.5823168427469,147.7224461170831, 150.8626360117056, 154.0028828430108, 157.1431832194004,160.2835340131002, 163.4239323343534, 166.5643755105532, 169.7048610649942,172.845386699872, 175.9859502802961, 179.1265498201947, 182.2671834695183,185.4078495028307, 188.5485463090021, 191.6892723819093, 194.8300263119383,197.9708067786496, 201.1116125437999, 204.2524424443358, 207.3932953878596,210.5341703467098, 213.6750663534141, 216.8159824962144, 219.9569179153474,223.0978717991966, 226.2388433810609, 229.3798319360387, 232.5208367783268,235.6618572582787, 238.8028927612871, 241.943942702743, 245.0850065294672,248.2260837157502, 251.3671737619164, 254.5082761929693, 257.649390556998,260.7905164237957, 263.9316533836523, 267.0728010459829, 270.2139590383147,273.3551270053282, 276.4963046076781, 279.6374915213833, 282.778687436733,285.9198920576064, 289.0611051007226, 292.2023262949574, 295.3435553808254,298.4847921092646, 301.6260362421629, 304.7672875509309, 307.908545816245,311.0498108278112, 314.1910823834387, 317.3323602891387, 320.4736443587389,323.6149344124281, 326.7562302783298, 329.8975317905115, 333.0388387894455,336.1801511215153, 339.3214686387974, 342.4627911987849, 345.6041186642261,348.7454509027197, 351.886787786553, 355.0281291925481, 358.1694750018927,361.3108250998642, 364.4521793757087, 367.593537722451, 370.7349000366847,373.8762662185854, 377.0176361714929, 380.1590098020788, 383.3003870200034,386.4417677378807, 389.5831518711236, 392.7245393378147, 395.8659300588318,399.0073239571865, 402.1487209585936, 405.2901209909211, 408.4315239843747,411.5729298711667, 414.7143385856815, 417.8557500639707, 420.9971642444068,424.1385810669848, 427.280000473465, 430.4214224073116, 433.5628468136392,436.7042736390664, 439.8457028320225, 442.9871343418989, 446.1285681200393,449.2700041186825, 452.4114422917058, 455.5528825939937, 458.6943249821001,461.8357694134989, 464.9772158465583, 468.1186642411863, 471.2601145581708,474.4015667592944, 477.5430208073148, 480.6844766661983, 483.8259343007499,486.9673936764074, 490.1088547598619, 493.2503175184598, 496.3917819204395,499.5332479347844, 502.6747155312872, 505.8161846805574, 508.9576553537711,512.0991275229361, 515.2406011606936, 518.3820762403626, 521.5235527358955,524.6650306218695, 527.8065098735448, 530.9479904665941, 534.089472377555,537.230955582856, 540.372440060392, 543.5139257879564, 546.655412744179,549.7969009076752, 552.9383902580064, 556.0798807751443, 559.2213724394408,562.3628652316395, 565.5043591329607, 568.6458541250273, 571.7873501898963,574.9288473099032, 578.0703454679418, 581.2118446471584, 584.3533448310872,587.4948460036311, 590.6363481490084, 593.777851251788, 596.9193552968669,600.0608602694423, 603.2023661550344, 606.343872939461, 609.4853806088179,612.6268891495113, 615.7683985482095, 618.9099087918319, 622.0514198676036,625.1929317629782, 628.3344444656705, 631.4759579636318, 634.6174722450585,637.7589872983816, 640.9005031122405, 644.0420196755233, 647.1835369774376,650.3250550070022, 653.4665737538566, 656.608093207812, 659.7496133586814,662.8911341965854, 666.0326557117976, 669.1741778947833, 672.3157007361975,675.4572242268097, 678.5987483576628, 681.7402731198679, 684.8817985047631,688.0233245038157, 691.1648511086532, 694.3063783110547, 697.4479061029703,700.5894344764122, 703.7309634236445, 706.8724929370167, 710.0140230089969,713.1555536322219, 716.2970847994138, 719.4386165035222, 722.5801487374497,725.7216814943702, 728.8632147675165, 732.004748550242, 735.1462828360155,738.2878176184463, 741.4293528911456, 744.5708886479713, 747.7124248829122,750.8539615897272, 753.9954987626505, 757.137036395925, 760.2785744838565,763.4201130207428, 766.5616520011193, 769.7031914195529, 772.8447312707053,775.9862715493073, 779.1278122502242, 782.2693533682929, 785.410894898573,788.5524368361317, 791.6939791761165, 794.835521913755, 797.9770650443427,801.1186085633075, 804.2601524660121, 807.4016967480407, 810.5432414049633,813.684786432449, 816.8263318262129, 819.9678775820473, 823.1094236958442,826.2509701634526, 829.3925169808963, 832.5340641442106, 835.675611649493,838.8171594928957, 841.9587076706794, 845.1002561790286, 848.2418050144058,851.3833541729987, 854.5249036513569, 857.6664534459915, 860.8080035534068,863.9495539702256, 867.0911046930156, 870.232655718518, 873.3742070434552,876.5157586645979, 879.6573105787785, 882.798862782908, 885.9404152738143,889.0819680485198, 892.2235211040246, 895.3650744373655, 898.5066280456265,901.6481819259338, 904.789736075486, 907.9312904914243, 911.0728451710367,914.2144001116022, 917.3559553104354, 920.4975107648907, 923.6390664723595,926.7806224302977, 929.9221786361142, 933.0637350873423, 936.2052917815032,939.3468487161582, 942.4884058889019, 945.6299632973918, 948.7715209392748,951.9130788121515, 955.0546369138292, 958.1961952420654, 961.3377537946167,964.4793125692906, 967.6208715639405, 970.762430776409, 973.9039902045986,977.0455498465103, 980.1871096998352, 983.3286697628015, 986.4702300333498,989.6117905094368, 992.7533511891697, 995.8949120706039, 999.0364731518365];
%!   Ai = 4 * sin(zetai) ./ (2 * zetai + sin(2 * zetai));
%!   theta_ref = zeros(numel(x), numel(sol.t));
%!   for i=1:numel(zetai)
%!     theta_ref += Ai(i) * exp(-zetai(i)^2 * lambda / (rho * cp * L^2) * sol.t) .* cos(zetai(i) * x / L);
%!   endfor
%!   theta_ref = theta0 * theta_ref;
%!   if (do_plot)
%!     for i=[1:10,11:1000:numel(sol.t)]
%!       figure("visible", "off");
%!       hold("on");
%!       plot(x, theta_ref(:, i), "-;reference;k");
%!       plot(x, sol.theta(idx_x, i), "-;solution;r");
%!       xlabel("x [m]");
%!       ylabel("theta [degC]");
%!       grid on;
%!       grid minor on;
%!       title(sprintf("t=%gs", sol.t(i)));
%!     endfor
%!   endif
%!   tol = 3e-3;
%!   assert_simple(max(max(abs(sol.theta(idx_x, 100:end) - theta_ref(:, 100:end)))) < tol * abs(theta0));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
