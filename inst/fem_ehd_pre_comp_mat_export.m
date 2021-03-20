## Copyright (C) 2016(-2020) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} fem_ehd_pre_comp_mat_export(@var{comp_mat}, @var{options}, @var{output_file})
## @deftypefnx {} fem_ehd_pre_comp_mat_export(@var{comp_mat}, @var{options}, @var{output_filename})
## Print a compliance matrix to a file which can be loaded by MBDyn's module-hydrodynamic_plain_bearing2.
##
## @var{comp_mat} @dots{} Return value from fem_ehd_pre_comp_mat_struct or fem_ehd_pre_comp_mat_unstruct.
##
## @var{options}.matrix_type @dots{} One of "modal", "nodal", "nodal substruct", "nodal substruct total".
##
## @var{output_file} @dots{} Open file descriptor where the compliance matrices are written to.
##
## @var{output_filename} @dots{} Output filename where the compliance matrices are written to.
##
## @seealso{fem_ehd_pre_comp_mat_struct, fem_ehd_pre_comp_mat_unstruct}
## @end deftypefn

function fem_ehd_pre_comp_mat_export(comp_mat, options, varargin)
  if (nargin < 1 || nargin > 3 || nargout > 0)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isfield(options, "matrix_type"))
    options.matrix_type = "default";
  endif

  switch (options.matrix_type)
    case "default"
      if (isfield(comp_mat, "C"))
        options.matrix_type = "nodal";
      else
        options.matrix_type = "modal";
      endif
  endswitch

  owns_fd = false;
  fd = -1;
  
  unwind_protect
    if (nargin < 3)
      fd = stdout;
    elseif (ischar(varargin{1}))
      owns_fd = true;
      [fd, msg] = fopen(varargin{1}, "w");
      if (fd == -1)
        error("failed to open file \"%s\": %s", varargin{1}, msg);
      endif
    else
      fd = varargin{1};
    endif

    file_format = 12;

    d = comp_mat.bearing_dimensions.bearing_diameter;
    w = comp_mat.bearing_dimensions.bearing_width;

    x = comp_mat.bearing_surf.grid_x;
    z = comp_mat.bearing_surf.grid_z;
    z -= (max(z) + min(z)) / 2;

    N0 = numel(x) * numel(z);
    N1 = (numel(x) - 1) * numel(z);

    if (numel(x) < 3)
      error("numel(comp_mat.bearing_surf.grid_x must be at least three")
    endif

    if (abs((x(end) - x(1)) / (d * pi) - 1) >  sqrt(eps))
      error("invalid range for comp_mat.bearing_surf.grid_x");
    endif

    ic = zeros(numel(x), numel(z), "int32");

    for i=1:numel(x)
      for j=1:numel(z)
        ic(i, j) = (i - 1) * length(z) + j;
      endfor
    endfor

    ix = [(1:numel(x) + 1) + 2, 1:2];
    iz = 1:numel(z);
    x = [x([end - (2:-1:1)]) - d * pi, x(1:end - 1), x(1:2) + d * pi];
    ic = ic([1:end - 1, 1:2, end - 2:end - 1], :);

    fprintf(fd, "file format:\t%d\n", file_format);
    fprintf(fd, "\nbearing diameter:\t%.16g\n", d);
    fprintf(fd, "\nbearing width:\t%.16g\n", w);
    fprintf(fd, "\ncircumferential grid:\t%d\n", numel(ix));

    for i=1:numel(x)
      fprintf(fd, "%d\t%.16g\n", i, x(i));
    endfor

    fprintf(fd, "\naxial grid:\t%d\n", numel(z));

    for i=1:numel(z)
      fprintf(fd, "%d\t%.16g\n", i, z(i));
    endfor

    fprintf(fd, "\nnodes:\t%d\t%d\n", numel(ic), N1);
    idx = int32(0);

    for i=1:rows(ic)
      for j=1:columns(ic)
        fprintf(fd, "%d\t%d\t%d\t%d\n", ++idx, ic(i, j), ix(i), iz(j));
      endfor
    endfor

    fprintf(fd, "\nreference pressure:\t%e\n", comp_mat.reference_pressure);

    warning("error", "Octave:singular-matrix", "local");
    warning("error", "Octave:nearly-singular-matrix", "local");

    switch (options.matrix_type)
      case "nodal"
        if (~isfield(comp_mat, "C"))
          error("missing field \"C\"");
        endif

        if (rows(comp_mat.C) ~= N0 || columns(comp_mat.C) ~= N0)
          error("invalid size for comp_mat.C");
        endif

        C = comp_mat.C(1:end - numel(z), 1:end - numel(z));

        print_matrix(fd, C, "compliance matrix");

        fprintf(fd, "\n\n## cond(C)=%e\n", cond(C));
      case "nodal substruct"
        if (~isfield(comp_mat, "C"))
          error("missing field \"C\"");
        endif

        if (~isfield(comp_mat, "D"))
          error("missing field \"D\"");
        endif

        if (~isfield(comp_mat, "E"))
          error("missing field \"E\"");
        endif

        if (rows(comp_mat.C) ~= N0 || columns(comp_mat.C) ~= N0)
          error("invalid size for comp_mat.C");
        endif

        if (rows(comp_mat.D) ~= N0 )
          error("invalid size for comp_mat.D");
        endif

        if (columns(comp_mat.E) ~= N0)
          error("invalid size for comp_mat.E");
        endif

        C = comp_mat.C(1:end - numel(z), 1:end - numel(z));

        D = comp_mat.D(1:end - numel(z), :);
        E = comp_mat.E(:, 1:end - numel(z));

        print_matrix(fd, C, "compliance matrix substruct");

        fprintf(fd, "\n\n## cond(C)=%e\n", cond(C));

        print_matrix(fd, D, "substruct contrib matrix");
        print_matrix(fd, E, "substruct residual matrix");
      case {"nodal substruct total", "modal substruct total"}
        if (~isfield(comp_mat, "D"))
          error("missing field \"D\"");
        endif

        if (~isfield(comp_mat, "E"))
          error("missing field \"E\"");
        endif

        if (rows(comp_mat.D) ~= N0 )
          error("invalid size for comp_mat.D");
        endif

        if (columns(comp_mat.E) ~= N0)
          error("invalid size for comp_mat.E");
        endif

        D = comp_mat.D(1:end - numel(z), :);
        E = comp_mat.E(:, 1:end - numel(z));
	
        print_matrix(fd, D, "substruct total contrib matrix");
        print_matrix(fd, E, "substruct total residual matrix");
      case "modal"
        if (~isfield(comp_mat, "KPhi"))
          error("missing field \"KPhi\"");
        endif

        if (~isfield(comp_mat, "RPhi"))
          error("missing field \"RPhi\"");
        endif

        if (~isfield(comp_mat, "Phin"))
          error("missing field \"Phin\"");
        endif

        if (rows(comp_mat.Phin) ~= N0)
          error("invalid size for comp_mat.Phin");
        endif

        if (columns(comp_mat.RPhi) ~= N0)
          error("invalid size for comp_mat.RPhi");
        endif

        Phin = comp_mat.Phin(1:end - numel(z), :);
        KPhi = comp_mat.KPhi;
        RPhi = comp_mat.RPhi(:, 1:end - numel(z));

        print_matrix(fd, Phin, "mode shapes");
        print_matrix(fd, KPhi \ RPhi, "modal load");

        fprintf(fd, "\n\n## cond(KPhi)=%e\n", cond(KPhi));
      otherwise
	error("unknown option matrix_type=\"%s\"", options.matrix_type);
    endswitch

  unwind_protect_cleanup
    if (owns_fd && fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

function print_matrix(fd, C, C_name)
  fprintf(fd, "\n%s:\t%d\t%d\n", C_name, rows(C), columns(C));

  for i=1:rows(C)
    fprintf(fd, "%.16e\t", C(i, :));
    fprintf(fd, "\n");
  endfor
endfunction

%!demo
%! ############################################################################################################################
%! ## EHD TEST CASE according
%! ## Freund, Norman Owen, A thermo-elasto-hydrodynamic study of journal bearings, Doctor of Philosophy thesis, University of
%! ## Wollongong. Dept. of Mechanical Engineering, University of Wollongong, 1995. http://ro.uow.edu.au/theses/1568
%! ############################################################################################################################
%!
%! close all;
%! ## Figure 46b, page 175
%! ref_data_w = [  0, 195;
%!                10, 122;
%!                20,  65;
%!                30,  28;
%!                40,  18;
%!                50,  20;
%!                60,  30;
%!                70,  40;
%!                80,  48;
%!                90,  50;
%!               100,  60;
%!               110,  45;
%!               120,  30;
%!               130,   0];
%! ref_p = [  0,     0;
%!           10,   100;
%!           20,   200;
%!           30,   300;
%!           40,   375;
%!           50,   425;
%!           60,   525;
%!           70,   650;
%!           80,   750;
%!           90,   875;
%!          100,   975;
%!          110,  1100;
%!          120,  1175;
%!          130,  1175;
%!          140,  1050;
%!          150,   925;
%!          160,   700;
%!          170,   450;
%!          180,   200;
%!          190,    50;
%!          200,     0];
%! output_file = "";
%! have_mesh_size_binary = false;
%! unwind_protect
%!   output_file = tempname();
%!   if (ispc())
%!     output_file(output_file == "\\") = "/";
%!   endif
%!   [status, output] = shell("which fem_pre_mesh_size", true);
%!   if (status == 0)
%!     have_mesh_size_binary = true;
%!   endif
%!   if (~have_mesh_size_binary)
%!     error("fem_pre_mesh_size was not installed\nrun ./configure && make install inside the src directory");
%!   endif
%!   SI_unit_meter = 1;
%!   SI_unit_kilogram = 1;
%!   SI_unit_second = 1;
%!   SI_unit_kelvin = 1;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   ## Table 6, page 154
%!   ## Table 7, page 156
%!   param.E = 200000e6 / SI_unit_pascal;
%!   param.nu = 0.3;
%!   param.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.Di = 200e-3 / SI_unit_meter;       %! inner diameter shell (named D according to Freund)
%!   param.Do = 299e-3 / SI_unit_meter;       %! outer diameter bearing
%!   param.Wo = 200e-3 / SI_unit_meter;       %! total widht of bearing (named L according to Freund)
%!   param.d1 = 1e-3 / SI_unit_meter;         %! diaphragm thickness (1mm-5mm)
%!   param.d2 = 64e-3 / SI_unit_meter;        %! diaphragm half free edge width
%!   param.d3 = 132 * pi / 180 * param.Di / 2;# diaphragm length
%!   param.d4 = 7.5e-3 / SI_unit_meter;       %! diaphragm half end width
%!   param.d5 = 0e-3 / SI_unit_meter;         %! diaphragm straight extent
%!   param.cr = 150e-6 / SI_unit_meter;       %! bearing radial clearance
%!   param.etal = 0.00929 / (SI_unit_pascal * SI_unit_second); %! dynamic viscosity lubricant
%!   param.fact_etav = 1;                     %! etav/etal
%!   param.rhol = 844 / (SI_unit_kilogram / SI_unit_meter^3); %! density lubricant
%!   param.betal = 2.4e9 / SI_unit_pascal;    %! bulk modulus
%!   param.pcl = 0 / SI_unit_pascal;          %! cavitation pressure
%!   param.pin = 0 / SI_unit_pascal;        %! oil feed pressure
%!   param.pside = 0 / SI_unit_pascal;      %! side pressure
%!   param.Tl = 40 / SI_unit_kelvin;          %! liquid temperature
%!   param.d7 = 1e-3 / SI_unit_meter;         %! circumferential length of oil supply slot
%!   param.U1 = 28 / (SI_unit_meter / SI_unit_second); %! fluid velocity at the journal surface (named U2 according to Freund)
%!   param.U2 = 0 / (SI_unit_meter / SI_unit_second); %! fluid velocity at the shell surface (named U1 according to Freund)
%!   param.delta = 180 * pi / 180;            %! position of eccentricity
%!   param.epsilon = 0.3;                     %! relative eccentricity
%!   param.diaph_sec_n = 72;                  %! number of cross sections for diaphragm
%!   param.h = 15e-3 / SI_unit_meter;         %! mesh size for hydraulic mesh
%!   param.h1 = 4e-3 / SI_unit_meter;         %! mesh size in the area of the diaphragm
%!   param.h2 = 20e-3 / SI_unit_meter;        %! mesh size outside the are of the diaphragm
%!   param.ht = 10e-3 / SI_unit_meter;        %! mesh transition region
%!   param.pref = 1e6 / SI_unit_pascal;       %! reference pressure
%!   param.damp_alpha = 0 / (1 / SI_unit_second);  %! mass damping factor
%!   param.damp_beta = 0 / SI_unit_second;         %! stiffness damping factor
%!   n = 2;
%!   k = 72;
%!   param.t1 = n * param.Di * pi / param.U1;
%!   param.dt = param.t1 / (n * k);
%!   if (~have_mesh_size_binary)
%!     opt_mesh.mesh.element_size = param.h;
%!   endif
%!   opt_mesh.mesh.jacobian_range = [0.5, 1.5];
%!   opt_mesh.verbose = false;
%!   opt_mesh.output_file = [output_file, "_msh"];
%!   options.geo_tol = sqrt(eps);
%!   options_mbdyn.mbdyn_command = "mbdyn";

%!   group_defs(end + 1).id = 1;
%!   group_defs(end).name = "node_id_shell_support";
%!   group_defs(end).R = eye(3);
%!   group_defs(end).X0 = zeros(3, 1);
%!   group_defs(end).Xi = zeros(3, 1);
%!   group_defs(end).type = "cylinder";
%!   group_defs(end).geometry.rmin = 0.5 * param.Do;
%!   group_defs(end).geometry.rmax = 0.5 * param.Do;
%!   group_defs(end).geometry.zmin = -0.5 * param.Wo;
%!   group_defs(end).geometry.zmax = 0.5 * param.Wo;
%!   group_defs(end).compliance_matrix.matrix_type = "none";

%!   group_defs(end + 1).id = 2;
%!   group_defs(end).name = "node_id_shell_bearing";
%!   group_defs(end).R = eye(3);
%!   group_defs(end).X0 = zeros(3, 1);
%!   group_defs(end).Xi = zeros(3, 1);
%!   group_defs(end).type = "cylinder";
%!   group_defs(end).geometry.rmin = 0.5 * param.Di;
%!   group_defs(end).geometry.rmax = 0.5 * param.Di;
%!   group_defs(end).geometry.zmin = -0.5 * param.Wo;
%!   group_defs(end).geometry.zmax = 0.5 * param.Wo;
%!   group_defs(end).compliance_matrix.matrix_type = "nodal substruct";
%!   group_defs(end).compliance_matrix.bearing_type = "shell";
%!   group_defs(end).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs(end).compliance_matrix.reference_pressure = param.pref;
%!   group_defs(end).compliance_matrix.mesh_size = param.h;
%!   group_defs(end).bearing = "elem_id_bearing";
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([output_file, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [output_file, ".geo"]);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.OCCUnionUnify = 0;\n");
%!     fputs(fd, "Wd_0 = 2 * d2;\n");
%!     fputs(fd, "Wd_1 = 2 * d4;\n");
%!     fputs(fd, "td = d1;\n");
%!     fputs(fd, "dxo = d7;\n");
%!     fputs(fd, "Phid_0 = 0 / (0.5 * Di);\n");
%!     fputs(fd, "Phid_1 = d5 / (0.5 * Di);\n");
%!     fputs(fd, "Phid_2 = d3 / (0.5 * Di);\n");
%!     fputs(fd, "hd = (Do - Di) / 2. - td;\n");
%!     fputs(fd, "Phio = dxo / (0.5 * Di);\n");
%!     fputs(fd, "diaph_sec_m = 4;\n");
%!     fputs(fd, "diaph_sec_p[] = {};\n");
%!     fputs(fd, "For i In {0:diaph_sec_n - 1}\n");
%!     fputs(fd, "    j = i * diaph_sec_m;\n");
%!     fputs(fd, "    alpha = i / (diaph_sec_n - 1);\n");
%!     fputs(fd, "    Phii = (Phid_0 + 0.5 * Phio + Phid_2 * alpha);   \n");
%!     fputs(fd, "    If (Phid_2 * alpha >= Phid_1)\n");
%!     fputs(fd, "        Wdi = Wd_0 + (Wd_1 - Wd_0) * (Phid_2 * alpha - Phid_1) / (Phid_2 - Phid_1);\n");
%!     fputs(fd, "    Else\n");
%!     fputs(fd, "        Wdi = Wd_0;\n");
%!     fputs(fd, "    EndIf\n");
%!     fputs(fd, "    diaph_sec_p[j + 0] = newp;\n");
%!     fputs(fd, "    Point(diaph_sec_p[j + 0]) = {(0.5 * Di + td) * Cos(Phii), (0.5 * Di + td) * Sin(Phii), -0.5 * Wdi};\n");
%!     fputs(fd, "    diaph_sec_p[j + 1] = newp;\n");
%!     fputs(fd, "    Point(diaph_sec_p[j + 1]) = {(0.5 * Do + hd) * Cos(Phii), (0.5 * Do + hd) * Sin(Phii), -0.5 * Wdi};\n");
%!     fputs(fd, "    diaph_sec_p[j + 2] = newp;\n");
%!     fputs(fd, "    Point(diaph_sec_p[j + 2]) = {(0.5 * Do + hd) * Cos(Phii), (0.5 * Do + hd) * Sin(Phii), 0.5 * Wdi};\n");
%!     fputs(fd, "    diaph_sec_p[j + 3] = newp;\n");
%!     fputs(fd, "    Point(diaph_sec_p[j + 3]) = {(0.5 * Di + td) * Cos(Phii), (0.5 * Di + td) * Sin(Phii), 0.5 * Wdi};\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "bearing_sec_p[] = {};\n");
%!     fputs(fd, "bearing_sec_p[0] = newp;\n");
%!     fputs(fd, "Point(bearing_sec_p[0]) = {0.5 * Di, 0., -0.5 * Wo};\n");
%!     fputs(fd, "bearing_sec_p[1] = newp;\n");
%!     fputs(fd, "Point(bearing_sec_p[1]) = {0.5 * Do, 0., -0.5 * Wo};\n");
%!     fputs(fd, "bearing_sec_p[2] = newp;\n");
%!     fputs(fd, "Point(bearing_sec_p[2]) = {0.5 * Do, 0., 0.5 * Wo};\n");
%!     fputs(fd, "bearing_sec_p[3] = newp;\n");
%!     fputs(fd, "Point(bearing_sec_p[3]) = {0.5 * Di, 0., 0.5 * Wo};\n");
%!     fputs(fd, "Phi = Phid_0 - 0.5 * Phio;\n");
%!     fputs(fd, "ossl_sec_p[] = {};\n");
%!     fputs(fd, "ossl_sec_p[0] = newp;\n");
%!     fputs(fd, "Point(ossl_sec_p[0]) = {(0.5 * Di - td) * Cos(Phi), (0.5 * Di - td) * Sin(Phi), -0.5 * Wd_0};\n");
%!     fputs(fd, "ossl_sec_p[1] = newp;\n");
%!     fputs(fd, "Point(ossl_sec_p[1]) = {(0.5 * Do + hd) * Cos(Phi), (0.5 * Do + hd) * Sin(Phi), -0.5 * Wd_0};\n");
%!     fputs(fd, "ossl_sec_p[2] = newp;\n");
%!     fputs(fd, "Point(ossl_sec_p[2]) = {(0.5 * Do + hd) * Cos(Phi), (0.5 * Do + hd) * Sin(Phi), 0.5 * Wd_0};\n");
%!     fputs(fd, "ossl_sec_p[3] = newp;\n");
%!     fputs(fd, "Point(ossl_sec_p[3]) = {(0.5 * Di - td) * Cos(Phi), (0.5 * Di - td) * Sin(Phi), 0.5 * Wd_0};\n");
%!     fputs(fd, "diaph_sec_l[] = {};\n");
%!     fputs(fd, "For i In {0:diaph_sec_n - 1}\n");
%!     fputs(fd, "    For j In {0:diaph_sec_m - 1}\n");
%!     fputs(fd, "        k0 = i * diaph_sec_m + j;\n");
%!     fputs(fd, "        If (j == diaph_sec_m - 1)\n");
%!     fputs(fd, "           k1 = i * diaph_sec_m;\n");
%!     fputs(fd, "        Else\n");
%!     fputs(fd, "           k1 = k0 + 1;\n");
%!     fputs(fd, "        EndIf\n");
%!     fputs(fd, "        diaph_sec_l[k0] = newreg;\n");
%!     fputs(fd, "        Line(diaph_sec_l[k0]) = {diaph_sec_p[k0], diaph_sec_p[k1]};\n");
%!     fputs(fd, "    EndFor\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "bearing_sec_l[] = {};\n");
%!     fputs(fd, "For i In {0:#bearing_sec_p[] - 1} \n");
%!     fputs(fd, "    bearing_sec_l[i] = newreg;\n");
%!     fputs(fd, "    If (i == #bearing_sec_p[] - 1)\n");
%!     fputs(fd, "       j = 0;\n");
%!     fputs(fd, "    Else\n");
%!     fputs(fd, "       j = i + 1;\n");
%!     fputs(fd, "    EndIf\n");
%!     fputs(fd, "    Line(bearing_sec_l[i]) = {bearing_sec_p[i], bearing_sec_p[j]};\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "ossl_sec_l[] = {};\n");
%!     fputs(fd, "For i In {0:#ossl_sec_p[] - 1} \n");
%!     fputs(fd, "    ossl_sec_l[i] = newreg;\n");
%!     fputs(fd, "    If (i == #ossl_sec_p[] - 1)\n");
%!     fputs(fd, "       j = 0;\n");
%!     fputs(fd, "    Else\n");
%!     fputs(fd, "       j = i + 1;\n");
%!     fputs(fd, "    EndIf\n");
%!     fputs(fd, "    Line(ossl_sec_l[i]) = {ossl_sec_p[i], ossl_sec_p[j]};\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "diaph_sec_ll[] = {};\n");
%!     fputs(fd, "For i In {0:diaph_sec_n - 1}\n");
%!     fputs(fd, "    diaph_sec_ll[i] = newreg;\n");
%!     fputs(fd, "    Line Loop(diaph_sec_ll[i]) = {diaph_sec_l[i * diaph_sec_m],\n");
%!     fputs(fd, "                                  diaph_sec_l[i * diaph_sec_m + 1],\n");
%!     fputs(fd, "                                  diaph_sec_l[i * diaph_sec_m + 2],\n");
%!     fputs(fd, "                                  diaph_sec_l[i * diaph_sec_m + 3]};\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "bearing_sec_ll = newreg;\n");
%!     fputs(fd, "Line Loop(bearing_sec_ll) = {bearing_sec_l[]};\n");
%!     fputs(fd, "ossl_sec_ll = newreg;\n");
%!     fputs(fd, "Line Loop(ossl_sec_ll) = {ossl_sec_l[]};\n");
%!     fputs(fd, "bearing_sec_s = newreg;\n");
%!     fputs(fd, "Plane Surface(bearing_sec_s) = {bearing_sec_ll};\n");
%!     fputs(fd, "ossl_sec_s = newreg;\n");
%!     fputs(fd, "Plane Surface(ossl_sec_s) = {ossl_sec_ll};\n");
%!     fputs(fd, "diaph_v = newv;\n");
%!     fputs(fd, "ThruSections(diaph_v) = {diaph_sec_ll[]};\n");
%!     fputs(fd, "For i In {0:#diaph_sec_ll[] - 1}\n");
%!     fputs(fd, "    Recursive Delete { Curve{diaph_sec_ll[i]}; }\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "For i In {0:#diaph_sec_l[] - 1}\n");
%!     fputs(fd, "    Recursive Delete { Curve{diaph_sec_l[i]}; }\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "For i In {0:#diaph_sec_p[] - 1}\n");
%!     fputs(fd, "    Recursive Delete { Curve{diaph_sec_p[i]}; }\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "bearing_rot_v1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{bearing_sec_s}; };\n");
%!     fputs(fd, "bearing_rot_v2[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{bearing_rot_v1[0]}; };\n");
%!     fputs(fd, "bearing_rot_v3[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{bearing_rot_v2[0]}; };\n");
%!     fputs(fd, "bearing_rot_v4[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{bearing_rot_v3[0]}; };\n");
%!     fputs(fd, "ossl_rot_v[] = Extrude {{0, 0, 1}, {0, 0, 0}, Phio} { Surface{ossl_sec_s}; };\n");
%!     fputs(fd, "bearing_v5 = newv;\n");
%!     fputs(fd, "BooleanUnion(bearing_v5) = {Volume{bearing_rot_v1[1]}; Delete; }{ Volume{bearing_rot_v2[1], bearing_rot_v3[1], bearing_rot_v4[1]}; Delete; };\n");
%!     fputs(fd, "bearing_v = newv;\n");
%!     fputs(fd, "BooleanDifference(bearing_v) = {Volume{bearing_v5}; Delete; }{ Volume{diaph_v, ossl_rot_v[1]}; Delete; };\n");
%!     fputs(fd, "bearing_bnd[] = Unique(Abs(Boundary{Volume{bearing_v};}));\n");
%!     fputs(fd, "Physical Volume(\"bearing\", 1) = {bearing_v};\n");
%!     fputs(fd, "For i In {0:#bearing_bnd[] - 1}\n");
%!     fputs(fd, "    Physical Surface(i) = {bearing_bnd[i]};\n");
%!     fputs(fd, "EndFor\n");
%!     if (have_mesh_size_binary)
%!       fputs(fd, "iCurrField = 0;\n");
%!       fputs(fd, "iCurrField++;\n");
%!       fputs(fd, "Field[iCurrField] = ExternalProcess;\n");
%!       fputs(fd, "Field[iCurrField].CommandLine = Sprintf(\"fem_pre_mesh_size diaphragm Wo=%g Di=%g Do=%g d1=%g d2=%g d3=%g d4=%g d7=%g h1=%g h2=%g ht=%g\", Wo, Di, Do, 1.5 * d1, d2, d3, d4, d7, h1, h2, ht);\n");
%!       fputs(fd, "Background Field = iCurrField;\n");
%!     endif
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect

%!   mesh = fem_pre_mesh_unstruct_create([output_file, ".geo"], param, opt_mesh);
%!   mesh.groups.tria6 = fem_pre_mesh_groups_create(mesh, group_defs, options.geo_tol).tria6;
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(find([mesh.groups.tet10.id == 1])).elements) = 1;
%!   mesh.material_data.rho = param.rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(param.E, param.nu);

%!   cms_opt.invariants = true;
%!   cms_opt.refine_max_iter = int32(0);
%!   cms_opt.number_of_threads = int32(4);
%!   cms_opt.verbose = false;
%!   cms_opt.modes.number = 0;
%!   cms_opt.element.name = "elem_id_diaphragm_cms";

%!   node_set = int32(rows(mesh.nodes) + (1:numel(group_defs)));
%!   node_names = {group_defs.name};
%!   cms_opt.nodes.modal.number = node_set(1);
%!   cms_opt.nodes.modal.name = node_names{1};

%!   for j=1:numel(node_set) - 1
%!     cms_opt.nodes.interfaces(j).number = node_set(j + 1);
%!     cms_opt.nodes.interfaces(j).name = node_names{j + 1};
%!   endfor

%!   idx_grp_itf = find([group_defs.id] > 1);
%!   idx_grp_modal = find([group_defs.id] == 1);

%!   mesh.nodes([cms_opt.nodes.interfaces.number], 1:3) = [group_defs(idx_grp_itf).Xi].';
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = group_defs(idx_grp_modal).Xi.';

%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [group_defs(idx_grp_itf).id], node_set(idx_grp_itf));

%!   idx_grp_outer = find([mesh.groups.tria6.id] == 1);

%!   load_case.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.tria6(idx_grp_outer).nodes, 1:3) = true; %! Norman Owen Freund 1995, page 8
%!   load_case.locked_dof(cms_opt.nodes.modal.number, 1:6) = true;

%!   bearing_surf = repmat(struct("group_idx", [], "X0", [], "R", [], "options", [], "name", [], "bearing", []), 1, numel(group_defs));
%!   num_comp_mat = int32(0);

%!   for j=1:numel(group_defs)
%!     switch group_defs(j).compliance_matrix.matrix_type
%!       case "none"
%!       otherwise
%!         ++num_comp_mat;
%!         bearing_surf(num_comp_mat).group_idx = find([mesh.groups.tria6.id] == group_defs(j).id);
%!         bearing_surf(num_comp_mat).name = group_defs(j).name;
%!         bearing_surf(num_comp_mat).bearing = group_defs(j).bearing;
%!         bearing_surf(num_comp_mat).X0 = group_defs(j).X0;
%!         bearing_surf(num_comp_mat).R = group_defs(j).R;
%!         bearing_surf(num_comp_mat).options = group_defs(j).compliance_matrix;
%!         bearing_surf(num_comp_mat).master_node_no = node_set(j);

%!         switch group_defs(j).type
%!           case "cylinder"
%!             bearing_surf(num_comp_mat).r = mean([group_defs(j).geometry.rmax, group_defs(j).geometry.rmin]);
%!             bearing_surf(num_comp_mat).w = group_defs(j).geometry.zmax - group_defs(j).geometry.zmin;
%!           otherwise
%!             error("bearing geometry type \"%s\" not implemented", group_defs(j).type);
%!         endswitch

%!         bearing_surf(num_comp_mat).nodes = mesh.groups.tria6(find([[mesh.groups.tria6].id] == group_defs(j).id)).nodes;
%!     endswitch
%!   endfor

%!   bearing_surf = bearing_surf(1:num_comp_mat);

%!   [load_case_pressure, bearing_surf] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf);

%!   load_case = fem_pre_load_case_merge(load_case, load_case_pressure);

%!   [mesh, ...
%!    mat_ass, ...
%!    dof_map, ...
%!    sol_eig, ...
%!    cms_opt] = fem_cms_create(mesh, load_case, cms_opt);

%!   mat_ass.Dred = param.damp_alpha * mat_ass.Mred + param.damp_beta * mat_ass.Kred;

%!   fem_cms_export([output_file, "_cms"], mesh, dof_map, mat_ass, cms_opt);

%!   comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, ...
%!                                            mat_ass, ...
%!                                            dof_map, ...
%!                                            cms_opt, ...
%!                                            bearing_surf);

%!   for j=1:numel(comp_mat)
%!     comp_mat_file = [output_file, "_", bearing_surf(j).bearing, ".dat"];
%!     fem_ehd_pre_comp_mat_export(comp_mat(j), bearing_surf(j).options, comp_mat_file);
%!   endfor

%!   unwind_protect
%!     fd = -1;

%!     [fd, msg] = fopen([output_file, "_shell.set"], "w");

%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", [output_file, "_shell.set"], msg);
%!     endif

%!     for j=1:numel(comp_mat)
%!       fprintf(fd, "set: number_of_nodes_x = %d;\n", numel(comp_mat(j).bearing_surf.grid_x) + 1);
%!       fprintf(fd, "set: number_of_nodes_z = %d;\n", numel(comp_mat(j).bearing_surf.grid_z));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect

%!   mbdyn_pre_write_param_file([output_file, ".set"], param);

%!   options_mbdyn.output_file = [output_file, "_mbd"];
%!   options_mbdyn.f_run_mbdyn2easyanim = false;
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([output_file, "_mbd.mbdyn"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [output_file, "_mbd.mbdyn"]);
%!     endif
%!     fputs(fd, "set: integer ref_id_assembly = 1001;\n");
%!     fputs(fd, "set: integer ref_id_shell_bearing = 1002;\n");
%!     fputs(fd, "set: integer ref_id_shell_support = 1003;\n");
%!     fputs(fd, "set: integer ref_id_journal_bearing = 1004;\n");
%!     fputs(fd, "set: integer ref_id_journal_support = 1005;\n");
%!     fputs(fd, "set: integer node_id_shell_support = 2001;\n");
%!     fputs(fd, "set: integer node_id_shell_bearing = 2002;\n");
%!     fputs(fd, "set: integer node_id_journal_bearing = 2004;\n");
%!     fputs(fd, "set: integer joint_id_shell_support = 3001;\n");
%!     fputs(fd, "set: integer joint_id_journal_support = 3002;\n");
%!     fputs(fd, "set: integer elem_id_diaphragm_cms = 3003;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 3005;\n");
%!     fputs(fd, "set: integer number_of_nodes_x;\n");
%!     fputs(fd, "set: integer number_of_nodes_z;\n");
%!     fprintf(fd, "include: \"%s.set\";\n", output_file);
%!     fprintf(fd, "include: \"%s_shell.set\";\n", output_file);
%!     fputs(fd, "set: real omega1 = U1 / (0.5 * Di - cr);\n");
%!     fputs(fd, "set: real omega2 = U2 / (0.5 * Di);\n");
%!     fputs(fd, "set: real omega = omega1 - omega2;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        linear solver: umfpack, cc, scale, row max column max, always, max iterations, 2500;\n");
%!     fputs(fd, "        nonlinear solver: line search, default solver options, moderate nonlinear, divergence check, no;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, minmax;\n");
%!     fputs(fd, "        max iterations: 100;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-3;\n");
%!     fputs(fd, "        derivatives max iterations: 5;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fputs(fd, "        output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       print: dof stats, to file;\n");
%!     fputs(fd, "       print: dof description, to file;\n");
%!     fputs(fd, "       print: equation description, to file;\n");
%!     fputs(fd, "       structural nodes: 3;\n");
%!     fputs(fd, "       joints: 3;\n");
%!     fputs(fd, "       loadable elements: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_assembly,\n");
%!     fputs(fd, "        position, reference, global, null,\n");
%!     fputs(fd, "        orientation, reference, global, eye,\n");
%!     fputs(fd, "        velocity, reference, global, null,\n");
%!     fputs(fd, "        angular velocity, reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_shell_bearing,\n");
%!     fputs(fd, "        position, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_assembly, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_assembly,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                omega2;\n");
%!     fputs(fd, "reference: ref_id_shell_support,\n");
%!     fputs(fd, "        position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_shell_bearing, null;\n");
%!     fputs(fd, "reference: ref_id_journal_bearing,\n");
%!     fputs(fd, "        position, reference, ref_id_assembly,\n");
%!     fputs(fd, "                  epsilon * cr * cos(delta),\n");
%!     fputs(fd, "                  epsilon * cr * sin(delta),\n");
%!     fputs(fd, "                  0.,\n");
%!     fputs(fd, "        orientation, reference, ref_id_assembly, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_assembly,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                omega1;\n");
%!     fputs(fd, "reference: ref_id_journal_support,\n");
%!     fputs(fd, "        position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_journal_bearing, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_shell_support, modal,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_support, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_support, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_shell_support, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_shell_support, null;\n");
%!     fputs(fd, "        structural: node_id_shell_bearing, static,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_shell_bearing, null;\n");
%!     fputs(fd, "        structural: node_id_journal_bearing, static,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_journal_bearing, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_shell_support, total pin joint,\n");
%!     fputs(fd, "                node_id_shell_support,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                mult, time, omega2;\n");
%!     fputs(fd, "        joint: joint_id_journal_support, total pin joint,\n");
%!     fputs(fd, "                node_id_journal_bearing,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                mult, time, omega1;\n");
%!     fprintf(fd, "        include: \"%s_cms.elm\";\n", output_file);
%!     fputs(fd, "        user defined: elem_id_bearing, hydrodynamic plain bearing2,\n");
%!     fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!     fputs(fd, "                density, rhol,\n");
%!     fputs(fd, "                betal,\n");
%!     fputs(fd, "                pressure, pcl,\n");
%!     fputs(fd, "                viscosity, etal,\n");
%!     fputs(fd, "                temperature, Tl,\n");
%!     fputs(fd, "            viscosity vapor, factor, fact_etav,\n");
%!     fputs(fd, "                mesh, linear finite difference,\n");
%!     fputs(fd, "                geometry, cylindrical,\n");
%!     fputs(fd, "                        mesh position, at bearing,\n");
%!     fputs(fd, "                        bearing width, Wo,\n");
%!     fputs(fd, "                        shaft diameter, Di - 2 * cr,\n");
%!     fputs(fd, "                        bearing diameter, Di,\n");
%!     fputs(fd, "                shaft node, node_id_journal_bearing,\n");
%!     fputs(fd, "                offset, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                bearing node, node_id_shell_bearing,\n");
%!     fputs(fd, "                offset, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                        number of nodes z, number_of_nodes_z,\n");
%!     fputs(fd, "                        number of nodes Phi, number_of_nodes_x,\n");
%!     fputs(fd, "                        boundary conditions,\n");
%!     fputs(fd, "                                        pressure, pside,\n");
%!     fputs(fd, "                                        pressure, pside,\n");
%!     fputs(fd, "                        lubrication grooves, 1,\n");
%!     fputs(fd, "                                at bearing,\n");
%!     fputs(fd, "                                        pressure, pin,\n");
%!     fputs(fd, "                                        position, 0., 0.,\n");
%!     fputs(fd, "                                        rectangle, width, d7, height, Wo,\n");
%!     fputs(fd, "                compliance model,\n");
%!     fprintf(fd, "                        matrix, 1, from file, \"%s_elem_id_bearing.dat\",\n", output_file);
%!     fputs(fd, "                        E1, E,\n");
%!     fputs(fd, "                        nu1, nu,\n");
%!     fputs(fd, "                        modal element, elem_id_diaphragm_cms,\n");
%!     fputs(fd, "                pressure dof scale, pref,\n");
%!     fputs(fd, "                reynolds equation scale, dt / (rhol * cr),\n");
%!     fputs(fd, "                elasticity equation scale, dt / cr,\n");
%!     fputs(fd, "                output pressure, yes,\n");
%!     fputs(fd, "                output stress, yes,\n");
%!     fputs(fd, "                output density, yes,\n");
%!     fputs(fd, "                output friction loss, yes,\n");
%!     fputs(fd, "                output clearance, yes,\n");
%!     fputs(fd, "                output reaction force, yes,\n");
%!     fputs(fd, "                output mesh, yes,\n");
%!     fputs(fd, "                output total deformation, yes,\n");
%!     fputs(fd, "                output, yes;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   mbdyn_solver_run([output_file, "_mbd.mbdyn"], options_mbdyn);

%!   res.log_dat = mbdyn_post_load_log(options_mbdyn.output_file);

%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbdyn.output_file);

%!   [res.elem_id, res.q, res.qdot, res.qddot] = mbdyn_post_load_output_mod(options_mbdyn.output_file, numel(res.t));

%!   res.bearings = mbdyn_post_ehd_load_output(options_mbdyn.output_file, res.log_dat);

%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);

%!   res.elem_idx_cms = find(res.elem_id == res.log_dat.vars.elem_id_diaphragm_cms);
%!   res.elem_idx_bearing = find(res.log_dat.vars.elem_id_bearing == [res.log_dat.bearings.label]);

%!   cms_data.mesh = mesh;
%!   cms_data.dof_map = dof_map;
%!   cms_data.cms_opt = cms_opt;
%!   res.sol_dyn = fem_post_cms_sol_import(options_mbdyn.output_file, cms_data);

%!   figure("visible", "off");
%!   plot(res.t * SI_unit_second, 1e-3 * max(res.bearings.columns.p_n, [], 2) * SI_unit_pascal, "-;max(p(t));1");
%!   xlabel("t [s]");
%!   ylabel("p [kPa]");
%!   grid on;
%!   grid minor on;
%!   title("convergence history of peak pressure");

%!   figure("visible", "off");
%!   plot(res.t * SI_unit_second, 1e6 * max(res.bearings.columns.wtot_n, [], 2) * SI_unit_pascal, "-;max(w(t));1");
%!   xlabel("t [s]");
%!   ylabel("w [um]");
%!   grid on;
%!   grid minor on;
%!   title("convergence history of peak deformation");

%!   figure("visible", "off");
%!   Phi = res.bearings.xi(1, :) / (0.5 * param.Di);
%!   p = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.p(:, :, end), res.bearings.xi(1, :), 0);
%!   w = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.wtot(:, :, end), res.bearings.xi(1, :), 0);
%!   ax = plotyy(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, 180 / pi * Phi, 1e6 * w * SI_unit_meter);
%!   xlabel("Phi [deg]");
%!   ylabel(ax(1), "p [kPa]");
%!   ylabel(ax(2), "w [um]");
%!   grid on;
%!   grid minor on;
%!   for i=1:2
%!     xlim(ax(i), [0, 360]);
%!   endfor
%!   xticks(0:30:360);
%!   title("midplane pressure and deformation");

%!   figure("visible","off");
%!   hold on;
%!   set(plot(180 / pi * Phi, 1e6 * w * SI_unit_meter, "-;w(Phi) [um];1"), "linewidth", 5);
%!   set(plot(ref_data_w(:, 1), ref_data_w(:, 2), "-;reference w(Phi) [um];0"), "linewidth", 3);
%!   xlabel("Phi [deg]");
%!   ylabel("w [um]");
%!   grid on;
%!   grid minor on;
%!   xlim([0,360]);
%!   xticks(0:30:360);
%!   title("midplane deformation");

%!   figure("visible","off");
%!   hold on;
%!   set(plot(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, "-;p(Phi) [kPa];1"), "linewidth", 5);
%!   set(plot(ref_p(:, 1), ref_p(:, 2), "-;reference p(Phi) [kPa];0"), "linewidth", 3);
%!   xlabel("Phi [deg]");
%!   ylabel("p* [kPa]");
%!   grid on;
%!   grid minor on;
%!   xlim([0,360]);
%!   xticks(0:30:360);
%!   title("midplane pressure");

%!   figure("visible","off");
%!   set(plot(180 / pi * Phi, p * param.cr^2 / (0.5 * param.Di * param.U1 * param.etal), "-;p*(Phi) [1];1"), "linewidth", 5);
%!   xlabel("Phi [deg]");
%!   ylabel("p* [1]");
%!   grid on;
%!   grid minor on;
%!   xlim([0,360]);
%!   xticks(0:30:360);
%!   title("midplane dimensionless pressure");

%!   figure("visible", "off");
%!   hold on;
%!   h = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.h(:, :, end), res.bearings.xi(1, :), 0);
%!   set(plot(180 / pi * Phi, 1e6 * h * SI_unit_meter, "-;h(Phi) [um];1"), "linewidth", 5);
%!   set(plot(180 / pi * Phi, 1e6 * w * SI_unit_meter, "-;w(Phi) [um];3"), "linewidth", 5);
%!   xlabel("Phi [deg]");
%!   ylabel(ax(2), "h, w [um]");
%!   grid on;
%!   grid minor on;
%!   xlim([0, 360]);
%!   title("midplane clearance and deformation");

%!   figure("visible","off");
%!   contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, 1e-3 * res.bearings.columns.p(:, :, end) * SI_unit_pascal);
%!   colormap jet;
%!   colorbar;
%!   xlabel("Phi [deg]");
%!   ylabel("z [mm]");
%!   grid on;
%!   grid minor on;
%!   title("pressure distribution p [kPa]");

%!   figure("visible","off");
%!   contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, 1e6 * res.bearings.columns.wtot(:, :, end) * SI_unit_meter);
%!   colormap jet;
%!   colorbar;
%!   xlabel("Phi [deg]");
%!   ylabel("z [mm]");
%!   grid on;
%!   grid minor on;
%!   title("radial deformation [um]");

%!   figure("visible","off");
%!   contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, res.bearings.columns.rho(:, :, end) / param.rhol);
%!   colormap jet;
%!   colorbar;
%!   xlabel("Phi [deg]");
%!   ylabel("z [mm]");
%!   grid on;
%!   grid minor on;
%!   title("fractional film content [1]");

%!   figure_list();
%!   p_int = interp1(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, ref_p(:, 1), "linear");
%!   w_int = interp1(180 / pi * Phi, 1e6 * w * SI_unit_meter, ref_data_w(:, 1), "linear");
%!   assert(p_int, ref_p(:, 2), 0.15 * max(abs(ref_p(:, 2))));
%!   assert(mean(abs(p_int - ref_p(:, 2))) < 0.05 * max(abs(ref_p(:, 2))));
%!   ## Don't check the first point because of issues related to the interpolation of the compliance matrix near to the slot
%!   assert(w_int(2:end), ref_data_w(2:end, 2), 0.07 * max(abs(ref_data_w(:, 2))));
%!   assert(mean(abs(w_int(2:end) - ref_data_w(2:end, 2))) < 0.04 * max(abs(ref_data_w(:, 2))));
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
