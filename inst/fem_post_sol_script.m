## Copyright (C) 2020(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} fem_post_sol_script(@var{post_pro_script_file}, @var{options})
## Create a Gmsh script named <@var{post_pro_script_file}> for post-processing within Gmsh.
## This function is called from fem_post_sol_external.
##
## @var{options}.scale_def @dots{} scale factor for deformations
##
## @var{options}.animation_delay @dots{} delay between time steps
##
## @var{options}.skin_only @dots{} create a skin only view in order to speed up animation
##
## @var{options}.show_time @dots{} display time steps (show_time == 1) or natural frequencies (show_time == 6)
##
## @var{options}.show_element @dots{} display element outlines
##
## @var{options}.stress_output @dots{} cell array of stress field names available in the result set
##
## @var{options}.mesh_file @dots{} mesh file name to be loaded into Gmsh
##
## @var{options}.result_filenames @dots{} cell array of result filenames to be loaded into Gmsh
##
## @var{options}.rotation_angle @dots{} Euler angles for default viewpoint
##
## @var{options}.scalar_min @dots{} lower limit for color-bar
##
## @var{options}.scalar_max @dots{} upper limit for color-bar
##
## @var{options}.print_to_file @dots{} generate JPEG files for each time step
##
## @var{options}.print_and_exit @dots{} exit Gmsh as soon as JPEG files have been created
##
## @end deftypefn

function fem_post_sol_script(post_pro_script_file, options)
  if (nargin ~= 2 || nargout > 0)
    print_usage();
  endif

  if (~isfield(options, "scale_def"))
    options.scale_def = 1;
  endif

  if (~isfield(options, "animation_delay"))
    options.animation_delay = 0.1;
  endif

  if (~isfield(options, "skin_only"))
    options.skin_only = true;
  endif

  if (~isfield(options, "show_time"))
    options.show_time = 1;
  endif

  if (~isfield(options, "stress_output"))
    options.stress_output = {};
  endif

  if (~isfield(options, "show_element"))
    options.show_element = true;
  endif

  if (~isfield(options, "intervals_type"))
    options.intervals_type = 3;
  endif

  fd = -1;
  
  unwind_protect
    [fd, msg] = fopen(post_pro_script_file, "w");

    if (fd == -1)
      error("failed to open file \"%s\"", post_pro_script_file);
    endif
    
    fprintf(fd, "Merge \"%s\";\n", gmsh_relative_path(options.mesh_filename));
    fprintf(fd, "Mesh.SurfaceEdges = 0;\n");
    fprintf(fd, "Mesh.VolumeFaces = 0;\n");
    fprintf(fd, "Mesh.SurfaceFaces = 0;\n");
    fprintf(fd, "Mesh.VolumeEdges = 0;\n");
    
    for j=1:numel(options.result_filenames)
      fprintf(fd, "Merge \"%s\";\n", gmsh_relative_path(options.result_filenames{j}));
    endfor

    iview = int32(0);
    iview_print = iview;
    
    fprintf(fd, "View[%d].Type = 1;\n", iview);
    fprintf(fd, "View[%d].VectorType = 5;\n", iview);
    fprintf(fd, "View[%d].Visible = %d;\n", iview, ~options.skin_only);
    fprintf(fd, "View[%d].DisplacementFactor = %g;\n", iview, options.scale_def);
    fprintf(fd, "PostProcessing.AnimationDelay = %g;\n", options.animation_delay);
    fprintf(fd, "View[%d].ShowTime = %d;\n", iview, options.show_time);
    fprintf(fd, "View[%d].ShowElement = %d;\n", iview, options.show_element);
    fprintf(fd, "View[%d].IntervalsType = %d;\n", iview, options.intervals_type);
    fprintf(fd, "View[%d].NbIso = 20;\n", iview);
    fprintf(fd, "View[%d].DrawSkinOnly=1;\n", iview);

    stress_type = {"tau", "taum", "vmis"};

    if (isfield(options, "stress_output"))
      for j=1:numel(stress_type)
        switch (stress_type{j})
          case options.stress_output
            ++iview;
            fprintf(fd, "View[%d].Type = 1;\n", iview);
            fprintf(fd, "View[%d].Visible = 0;\n", iview);
            fprintf(fd, "View[%d].ShowTime = %d;\n", iview, options.show_time);
            fprintf(fd, "View[%d].ShowElement = %d;\n", iview, options.show_element);
            fprintf(fd, "View[%d].IntervalsType = %d;\n", iview, options.intervals_type);
            fprintf(fd, "View[%d].NbIso = 20;\n", iview);
            fprintf(fd, "View[%d].DrawSkinOnly=1;\n", iview);
        endswitch
      endfor
    endif

    if (isfield(options, "rotation_angle"))
      fputs(fd, "General.Trackball = 0;\n");
      
      for j=1:3        
        fprintf(fd, "General.Rotation%s = %g;\n", {"X", "Y", "Z"}{j}, 180 / pi * options.rotation_angle(j));
      endfor
    endif
    
    if (options.skin_only)
      iview = int32(0);
      
      fprintf(fd, "Plugin(Skin).Visible = 0;\n");
      fprintf(fd, "Plugin(Skin).FromMesh = 0;\n");
      fprintf(fd, "Plugin(Skin).View = %d;\n", iview);
      fprintf(fd, "Plugin(Skin).Run;\n");

      inumstress = int32(0);

      if (isfield(options, "stress_output"))
        for j=1:numel(stress_type)
          switch (stress_type{j})
            case options.stress_output
              ++iview;
              ++inumstress;
              fprintf(fd, "Plugin(Skin).Visible = 0;\n");
              fprintf(fd, "Plugin(Skin).FromMesh = 0;\n");
              fprintf(fd, "Plugin(Skin).View = %d;\n", iview);
              fprintf(fd, "Plugin(Skin).Run;\n");
          endswitch
        endfor
      endif
      
      iview_print = ++iview;      
      fprintf(fd, "View[%d].Type = 1;\n", iview);
      fprintf(fd, "View[%d].VectorType = 5;\n", iview);
      fprintf(fd, "View[%d].Visible = 1;\n", iview);
      fprintf(fd, "View[%d].DisplacementFactor = %g;\n", iview, options.scale_def);
      fprintf(fd, "View[%d].ShowTime = %d;\n", iview, options.show_time);
      fprintf(fd, "View[%d].ShowElement = %d;\n", iview, options.show_element);
      fprintf(fd, "View[%d].IntervalsType = %d;\n", iview, options.intervals_type);
      fprintf(fd, "View[%d].NbIso = 20;\n", iview);
      fprintf(fd, "View[%d].DrawSkinOnly=1;\n", iview);
      fprintf(fd, "View[%d].TensorType = 1;\n", iview);
      fprintf(fd, "View[%d].ExternalView = %d;\n", iview, iview + inumstress);
      
      if (isfield(options, "scalar_max"))
        if (~isfield(options, "scalar_min"))
          options.scalar_min = 0;
        endif
        
        fprintf(fd, "View[%d].RangeType = 2;\n", iview);
        fprintf(fd, "View[%d].CustomMin = 0;\n", iview, options.scalar_min);
        fprintf(fd, "View[%d].CustomMax = %g;\n", iview, options.scalar_max);
        fprintf(fd, "View[%d].SaturateValues = 1;\n", iview);
      endif

      if (isfield(options, "stress_output"))
        for j=1:numel(stress_type)
          switch (stress_type{j})
            case options.stress_output
              ++iview;
              fprintf(fd, "View[%d].Visible = 0;\n", iview);
          endswitch
        endfor
      endif
    endif

    if (isfield(options, "print_to_file"))
      for j=1:options.num_time_steps
        fprintf(fd, "View[%d].TimeStep = %d;\n", iview_print, j - 1);
        output_file = make_absolute_filename(options.print_to_file);
        output_file(find(output_file == "\\")) = "/";
        fprintf(fd, "Print \"%s_%03d.jpg\";\n", output_file, j);
      endfor
      
      if (isfield(options, "print_and_exit") && options.print_and_exit)
        fprintf(fd, "Exit;\n");
      endif      
    endif
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

function fname = gmsh_relative_path(fname)
  [fdir, fname, fext] = fileparts(fname);
  fname = [fname, fext];
endfunction
