## fem_get_software_version.m:01
%!test
%! version = fem_get_software_version("gmsh");
%! assert(size(version), [1, 3]);
%! assert(version(1) >= 4);
