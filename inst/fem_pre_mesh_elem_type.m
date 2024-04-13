function eltype_out = fem_pre_mesh_elem_type()

  persistent eltype = [];

  if (isempty(eltype))
    empty_cell = cell(1, 31);

    struct("dim", empty_cell, "id", empty_cell, "name", empty_cell, "norder", empty_cell, "nordernop", empty_cell, "promote", empty_cell);

    eltype(1).dim = 3;
    eltype(1).id = 11;
    eltype(1).name = "tet10";
    eltype(1).norder = [1,2,3,4,5,6,7,8,10,9];
    eltype(1).nordernonp = [];
    eltype(1).promote = -1;

    eltype(2).dim = 2;
    eltype(2).id = 9;
    eltype(2).name = "tria6";
    eltype(2).norder = [1,2,3,4,5,6];
    eltype(2).nordernonp = [];
    eltype(2).promote = -1;

    eltype(3).dim = 3;
    eltype(3).id = 4;
    eltype(3).name = "tet4";
    eltype(3).norder = [4,4,4,4,1,2,3,3];
    eltype(3).nordernonp = [];
    eltype(3).promote = 6;

    eltype(4).dim = 3;
    eltype(4).id = 6;
    eltype(4).name = "penta6";
    eltype(4).norder = [6,4,4,5,3,1,1,2];
    eltype(4).nordernonp = [];
    eltype(4).promote = 6;

    eltype(5).dim = 2;
    eltype(5).id = 2;
    eltype(5).name = "tria3";
    eltype(5).norder = [1,2,3,3];
    eltype(5).nordernonp = [1,2,3];
    eltype(5).promote = 7;

    eltype(6).dim = 3;
    eltype(6).id = 5;
    eltype(6).name = "iso8";
    eltype(6).norder = [5,6,7,8,1,2,3,4];
    eltype(6).nordernonp = [];
    eltype(6).promote = -1;

    eltype(7).dim = 2;
    eltype(7).id = 3;
    eltype(7).name = "iso4";
    eltype(7).norder = [1,2,3,4];
    eltype(7).nordernonp = [];
    eltype(7).promote = -1;

    eltype(8).dim = 3;
    eltype(8).id = 17;
    eltype(8).name = "iso20";
    eltype(8).norder = [5,6,7,8,1,2,3,4,17,19,20,18,9,12,14,10,11,13,15,16];
    eltype(8).nordernonp = [];
    eltype(8).promote = -1;

    eltype(9).dim = 2;
    eltype(9).id = 16;
    eltype(9).name = "quad8";
    eltype(9).norder = [1,2,3,4,5,6,7,8];
    eltype(9).nordernonp = [];
    eltype(9).promote = -1;

    eltype(10).dim = 3;
    eltype(10).id = 18;
    eltype(10).name = "penta15";
    eltype(10).norder = [1,2,3,4,5,6,7,10,8,13,15,14,9,11,12];
    eltype(10).nordernonp = [];
    eltype(10).promote = -1;

    eltype(11).dim = 3;
    eltype(11).id = 11;
    eltype(11).name = "tet10h";
    eltype(11).norder = [1,2,3,4,5,6,7,8,10,9];
    eltype(11).nordernonp = [];
    eltype(11).promote = -1;

    eltype(12).dim = 3;
    eltype(12).id = 9;
    eltype(12).name = "tria6h";
    eltype(12).norder = [1,2,3,4,5,6];
    eltype(12).nordernonp = [];
    eltype(12).promote = -1;

    eltype(13).dim = 3;
    eltype(13).id = 7;
    eltype(13).name = "pyra5";
    eltype(13).norder = [5,5,5,5,1,2,3,4];
    eltype(13).nordernonp = [1,2,3,4,5];
    eltype(13).promote = 6;

    eltype(14).dim = 3;
    eltype(14).id = 29;
    eltype(14).name = "tet20";
    eltype(14).norder = [1,5,6,2,7,8,3,9,10,17,16,20,14,19,12,18,15,13,11,4];
    eltype(14).nordernonp = [];
    eltype(14).promote = -1;

    eltype(15).dim = 2;
    eltype(15).id = 21;
    eltype(15).name = "tria10";
    eltype(15).norder = [1,2,3,4,5,6,7,8,9,10];
    eltype(15).nordernonp = [];
    eltype(15).promote = -1;

    eltype(16).dim = 3;
    eltype(16).id = 17;
    eltype(16).name = "iso20r";
    eltype(16).norder = [1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16];
    eltype(16).nordernonp = [];
    eltype(16).promote = -1;

    eltype(17).dim = 2;
    eltype(17).id = 16;
    eltype(17).name = "quad8r";
    eltype(17).norder = [3,4,1,2,7,8,5,6];
    eltype(17).nordernonp = [];
    eltype(17).promote = -1;

    eltype(18).dim = 3;
    eltype(18).id = 12;
    eltype(18).name = "iso27";
    eltype(18).norder = [1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18,21,22,24,25,23,26,27];
    eltype(18).nordernonp = [];
    eltype(18).promote = -1;

    eltype(19).dim = 2;
    eltype(19).id = 10;
    eltype(19).name = "quad9";
    eltype(19).norder = [1,2,3,4,5,6,7,8,9];
    eltype(19).nordernonp = [];
    eltype(19).promote = -1;

    eltype(20).dim = 3;
    eltype(20).id = 12;
    eltype(20).name = "iso27upc";
    eltype(20).norder = [1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18,21,22,24,25,23,26,27];
    eltype(20).nordernonp = [];
    eltype(20).promote = -1;

    eltype(21).dim = 3;
    eltype(21).id = 17;
    eltype(21).name = "iso20upc";
    eltype(21).norder = [5,6,7,8,1,2,3,4,17,19,20,18,9,12,14,10,11,13,15,16];
    eltype(21).nordernonp = [];
    eltype(21).promote = -1;

    eltype(22).dim = 3;
    eltype(22).id = 5;
    eltype(22).name = "iso8upc";
    eltype(22).norder = [5,6,7,8,1,2,3,4];
    eltype(22).nordernonp = [];
    eltype(22).promote = -1;

    eltype(23).dim = 3;
    eltype(23).id = 17;
    eltype(23).name = "iso20upcr";
    eltype(23).norder = [1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16];
    eltype(23).nordernonp = [];
    eltype(23).promote = -1;

    eltype(24).dim = 3;
    eltype(24).id = 18;
    eltype(24).name = "penta15upc";
    eltype(24).norder = [1,2,3,4,5,6,7,10,8,13,15,14,9,11,12];
    eltype(24).nordernonp = [];
    eltype(24).promote = -1;

    eltype(25).dim = 3;
    eltype(25).id = 6;
    eltype(25).name = "penta6upc";
    eltype(25).norder = [6,4,4,5,3,1,1,2];
    eltype(25).nordernonp = [];
    eltype(25).promote = 22;

    eltype(26).dim = 1;
    eltype(26).id = 1;
    eltype(26).name = "line2";
    eltype(26).norder = [1,2];
    eltype(26).nordernonp = [];
    eltype(26).promote = -1;

    eltype(27).dim = 1;
    eltype(27).id = 8;
    eltype(27).name = "line3";
    eltype(27).norder = [1,2,3];
    eltype(27).nordernonp = [];
    eltype(27).promote = -1;

    eltype(28).dim = 0;
    eltype(28).id = 15;
    eltype(28).name = "point1";
    eltype(28).norder = 1;
    eltype(28).nordernonp = [];
    eltype(28).promote = -1;

    eltype(29).dim = 1;
    eltype(29).id = 26;
    eltype(29).name = "line4";
    eltype(29).norder = [1,2,3,4];
    eltype(29).nordernonp = [];
    eltype(29).promote = -1;

    eltype(30).dim = 3;
    eltype(30).id = 11;
    eltype(30).name = "tet10upc";
    eltype(30).norder = [1,2,3,4,5,6,7,8,10,9];
    eltype(30).nordernonp = [];
    eltype(30).promote = -1;

    eltype(32).dim = 3;
    eltype(32).id = 11;
    eltype(32).name = "tet10f";
    eltype(32).norder = [1,2,3,4,5,6,7,8,10,9];
    eltype(32).nordernonp = [];
    eltype(32).promote = -1;

    eltype(33).dim = 3;
    eltype(33).id = 4;
    eltype(33).name = "tet4f";
    eltype(33).norder = [4,4,4,4,1,2,3,3];
    eltype(33).nordernonp = [];
    eltype(33).promote = 35;

    eltype(34).dim = 3;
    eltype(34).id = 6;
    eltype(34).name = "penta6f";
    eltype(34).norder = [6,4,4,5,3,1,1,2];
    eltype(34).nordernonp = [];
    eltype(34).promote = 35;

    eltype(35).dim = 3;
    eltype(35).id = 5;
    eltype(35).name = "iso8f";
    eltype(35).norder = [5,6,7,8,1,2,3,4];
    eltype(35).nordernonp = [];
    eltype(35).promote = -1;

    eltype(36).dim = 3;
    eltype(36).id = 17;
    eltype(36).name = "iso20f";
    eltype(36).norder = [5,6,7,8,1,2,3,4,17,19,20,18,9,12,14,10,11,13,15,16];
    eltype(36).nordernonp = [];
    eltype(36).promote = -1;

    eltype(37).dim = 3;
    eltype(37).id = 18;
    eltype(37).name = "penta15f";
    eltype(37).norder = [1,2,3,4,5,6,7,10,8,13,15,14,9,11,12];
    eltype(37).nordernonp = [];
    eltype(37).promote = -1;

    eltype(38).dim = 3;
    eltype(38).id = 11;
    eltype(38).name = "tet10hf";
    eltype(38).norder = [1,2,3,4,5,6,7,8,10,9];
    eltype(38).nordernonp = [];
    eltype(38).promote = -1;

    eltype(39).dim = 3;
    eltype(39).id = 7;
    eltype(39).name = "pyra5f";
    eltype(39).norder = [5,5,5,5,1,2,3,4];
    eltype(39).nordernonp = [1,2,3,4,5];
    eltype(39).promote = 35;

    eltype(40).dim = 3;
    eltype(40).id = 17;
    eltype(40).name = "iso20fr";
    eltype(40).norder = [1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16];
    eltype(40).nordernonp = [];
    eltype(40).promote = -1;

    eltype(41).dim = 3;
    eltype(41).id = 12;
    eltype(41).name = "iso27f";
    eltype(41).norder = [1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18,21,22,24,25,23,26,27];
    eltype(41).nordernonp = [];
    eltype(41).promote = -1;
  endif

  eltype_out = eltype;
endfunction

%!test
%! eltype = fem_pre_mesh_elem_type();
%! assert_simple(numel(eltype), 31);
