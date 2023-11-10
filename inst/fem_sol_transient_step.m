## Copyright (C) 2021(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{U_1}, @var{dU_dt_1}, @var{dU_dt2_1}, @var{res}] = fem_sol_transient_init(@var{U_0}, @var{dU_dt_0}, @var{dU_dt2_0}, @var{R_1}, @var{data})
## Perform a time step
##
## @var{U_0} @dots{} Displacement at time step @var{t}
##
## @var{dU_dt_0} @dots{} Velocity at time step @var{t}
##
## @var{dU_dt2_0} @dots{} Acceleration at time step @var{t}
##
## @var{R_1} @dots{} Load vector at time step @var{t} + @var{data}.dt
##
## @var{U_1} @dots{} Displacement at time @var{t} + @var{data}.dt
##
## @var{dU_dt_1} @dots{} Velocity at time @var{t} + @var{data}.dt
##
## @var{dU_dt2_1} @dots{} Acceleration at time @var{t} + @var{data}.dt
##
## @var{res} @dots{} residual error of linear solver
##
## @end deftypefn

function [U_1, dU_dt_1, dU_dt2_1, res] = fem_sol_transient_step(U_0, dU_dt_0, dU_dt2_0, R_1, data)
  if (nargin ~= 5 || nargout < 3 || nargout > 4)
    print_usage();
  endif
  
  R_eff = R_1 + data.M * (data.a0 * U_0 + data.a2 * dU_dt_0 + data.a3 * dU_dt2_0) + data.D * (data.a1 * U_0 + data.a4 * dU_dt_0 + data.a5 * dU_dt2_0);
  
  U_1 = data.Keff_fact \ R_eff;

  if (nargout >= 4)
    Keff_U1 = data.Keff * U_1;
    res = norm(Keff_U1  - R_eff) / norm(Keff_U1 + R_eff);
  endif
  
  dU_dt2_1 = data.a0 * (U_1 - U_0) - data.a2 * dU_dt_0 - data.a3 * dU_dt2_0;
  dU_dt_1 = dU_dt_0 + data.a6 * dU_dt2_0 + data.a7 * dU_dt2_1;
endfunction
