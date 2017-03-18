function r = i4_sobol_generate ( m, n, skip )

%*****************************************************************************80
%
%% I4_SOBOL_GENERATE generates a Sobol dataset.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 December 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points to generate.
%
%    Input, integer SKIP, the number of initial points to skip.
%
%    Output, real R(M,N), the points.
%
  for j = 1 : n
    seed = skip + j - 1;
    [ r(1:m,j), seed ]  = i4_sobol ( m, seed );
  end

  return
end
