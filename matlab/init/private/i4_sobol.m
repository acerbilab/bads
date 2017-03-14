function [ quasi, seed ] = i4_sobol ( dim_num, seed )

%*****************************************************************************80
%
%% I4_SOBOL generates a new quasirandom Sobol vector with each call.
%
%  Discussion:
%
%    The routine adapts the ideas of Antonov and Saleev.
%
%    Thanks to Francis Dalaudier for pointing out that the range of allowed
%    values of DIM_NUM should start at 1, not 2!  17 February 2009.
%
%    This function was modified to use PERSISTENT variables rather than
%    GLOBAL variables, 13 December 2009.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 March 2012
%
%  Author:
%
%    Original FORTRAN77 version by Bennett Fox.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Antonov, Saleev,
%    USSR Computational Mathematics and Mathematical Physics,
%    Volume 19, 1980, pages 252 - 256.
%
%    Paul Bratley, Bennett Fox,
%    Algorithm 659:
%    Implementing Sobol's Quasirandom Sequence Generator,
%    ACM Transactions on Mathematical Software,
%    Volume 14, Number 1, pages 88-100, 1988.
%
%    Bennett Fox,
%    Algorithm 647:
%    Implementation and Relative Efficiency of Quasirandom 
%    Sequence Generators,
%    ACM Transactions on Mathematical Software,
%    Volume 12, Number 4, pages 362-376, 1986.
%
%    Ilya Sobol,
%    USSR Computational Mathematics and Mathematical Physics,
%    Volume 16, pages 236-242, 1977.
%
%    Ilya Sobol, Levitan, 
%    The Production of Points Uniformly Distributed in a Multidimensional 
%    Cube (in Russian),
%    Preprint IPM Akad. Nauk SSSR, 
%    Number 40, Moscow 1976.
%
%  Parameters:
%
%    Input, integer DIM_NUM, the number of spatial dimensions.
%    DIM_NUM must satisfy 1 <= DIM_NUM <= 40.
%
%    Input/output, integer SEED, the "seed" for the sequence.
%    This is essentially the index in the sequence of the quasirandom
%    value to be generated.  On output, SEED has been set to the
%    appropriate next value, usually simply SEED+1.
%    If SEED is less than 0 on input, it is treated as though it were 0.
%    An input value of 0 requests the first (0-th) element of the sequence.
%
%    Output, real QUASI(DIM_NUM), the next quasirandom vector.
%
  persistent atmost;
  persistent dim_max;
  persistent dim_num_save;
  persistent initialized;
  persistent lastq;
  persistent log_max;
  persistent maxcol;
  persistent poly;
  persistent recipd;
  persistent seed_save;
  persistent v;

  if ( isempty ( initialized ) )
    initialized = 0;
    dim_num_save = -1;
  end

  if ( ~initialized | dim_num ~= dim_num_save )

    initialized = 1;

    dim_max = 40;
    dim_num_save = -1;
    log_max = 30;
    seed_save = -1;
%
%  Initialize (part of) V.
%
    v(1:dim_max,1:log_max) = zeros(dim_max,log_max);

    v(1:40,1) = [ ...
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]';

    v(3:40,2) = [ ...
            1, 3, 1, 3, 1, 3, 3, 1, ...
      3, 1, 3, 1, 3, 1, 1, 3, 1, 3, ...
      1, 3, 1, 3, 3, 1, 3, 1, 3, 1, ...
      3, 1, 1, 3, 1, 3, 1, 3, 1, 3 ]';

    v(4:40,3) = [ ...
               7, 5, 1, 3, 3, 7, 5, ...
      5, 7, 7, 1, 3, 3, 7, 5, 1, 1, ...
      5, 3, 3, 1, 7, 5, 1, 3, 3, 7, ...
      5, 1, 1, 5, 7, 7, 5, 1, 3, 3 ]';

    v(6:40,4) = [ ...
                     1, 7, 9,13,11, ...
      1, 3, 7, 9, 5,13,13,11, 3,15, ...
      5, 3,15, 7, 9,13, 9, 1,11, 7, ...
      5,15, 1,15,11, 5, 3, 1, 7, 9 ]';
  
    v(8:40,5) = [ ...
                           9, 3,27, ...
     15,29,21,23,19,11,25, 7,13,17, ...
      1,25,29, 3,31,11, 5,23,27,19, ...
     21, 5, 1,17,13, 7,15, 9,31, 9 ]';

    v(14:40,6) = [ ...
              37,33, 7, 5,11,39,63, ...
     27,17,15,23,29, 3,21,13,31,25, ...
      9,49,33,19,29,11,19,27,15,25 ]';

    v(20:40,7) = [ ...
                                         13, ...
     33,115, 41, 79, 17, 29,119, 75, 73,105, ...
      7, 59, 65, 21,  3,113, 61, 89, 45,107 ]';

    v(38:40,8) = [ ...
                                7, 23, 39 ]';
%
%  Set POLY.
%
    poly(1:40)= [ ...
        1,   3,   7,  11,  13,  19,  25,  37,  59,  47, ...
       61,  55,  41,  67,  97,  91, 109, 103, 115, 131, ...
      193, 137, 145, 143, 241, 157, 185, 167, 229, 171, ...
      213, 191, 253, 203, 211, 239, 247, 285, 369, 299 ];

    atmost = 2^log_max - 1;
%
%  Find the number of bits in ATMOST.
%
    maxcol = i4_bit_hi1 ( atmost );
%
%  Initialize row 1 of V.
%
    v(1,1:maxcol) = 1;

  end
%
%  Things to do only if the dimension changed.
%
  if ( dim_num ~= dim_num_save )
%
%  Check parameters.
%
    if ( dim_num < 1 | dim_max < dim_num )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'I4_SOBOL - Fatal error!\n' );
      fprintf ( 1, '  The spatial dimension DIM_NUM should satisfy:\n' );
      fprintf ( 1, '    1 <= DIM_NUM <= %d\n', dim_max );
      fprintf ( 1, '  But this input value is DIM_NUM = %d\n', dim_num );
      return
    end

    dim_num_save = dim_num;
%
%  Initialize the remaining rows of V.
%
    for i = 2 : dim_num
%
%  The bits of the integer POLY(I) gives the form of polynomial I.
%
%  Find the degree of polynomial I from binary encoding.
%
      j = poly(i);
      m = 0;

      while ( 1 )

        j = floor ( j / 2 );

        if ( j <= 0 )
          break;
        end

        m = m + 1;

      end
%
%  Expand this bit pattern to separate components of the logical array INCLUD.
%
      j = poly(i);
      for k = m : -1 : 1
        j2 = floor ( j / 2 );
        includ(k) = ( j ~= 2 * j2 );
        j = j2;
      end
%
%  Calculate the remaining elements of row I as explained
%  in Bratley and Fox, section 2.
%
      for j = m + 1 : maxcol 
        newv = v(i,j-m);
        l = 1;
        for k = 1 : m
          l = 2 * l;
          if ( includ(k) )
            newv = bitxor ( newv, l * v(i,j-k) );
          end
        end
        v(i,j) = newv;
      end
    end
%
%  Multiply columns of V by appropriate power of 2.
%
    l = 1;
    for j = maxcol-1 : -1 : 1
      l = 2 * l;
      v(1:dim_num,j) = v(1:dim_num,j) * l;
    end
%
%  RECIPD is 1/(common denominator of the elements in V).
%
    recipd = 1.0 / ( 2 * l );

    lastq(1:dim_num) = 0;

  end

  seed = floor ( seed );

  if ( seed < 0 )
    seed = 0;
  end

  if ( seed == 0 )

    l = 1;
    lastq(1:dim_num) = 0;

  elseif ( seed == seed_save + 1 )
%
%  Find the position of the right-hand zero in SEED.
%
    l = i4_bit_lo0 ( seed );

  elseif ( seed <= seed_save )

    seed_save = 0;
    l = 1;
    lastq(1:dim_num) = 0;

    for seed_temp = seed_save : seed - 1
      l = i4_bit_lo0 ( seed_temp );
      for i = 1 : dim_num
        lastq(i) = bitxor ( lastq(i), v(i,l) );
      end
    end

    l = i4_bit_lo0 ( seed );

  elseif ( seed_save + 1 < seed )

    for seed_temp = seed_save + 1 : seed - 1
      l = i4_bit_lo0 ( seed_temp );
      for i = 1 : dim_num
        lastq(i) = bitxor ( lastq(i), v(i,l) );
      end
    end

    l = i4_bit_lo0 ( seed );

  end
%
%  Check that the user is not calling too many times!
%
  if ( maxcol < l )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4_SOBOL - Fatal error!\n' );
    fprintf ( 1, '  Too many calls!\n' );
    fprintf ( 1, '  MAXCOL = %d\n', maxcol );
    fprintf ( 1, '  L =      %d\n', l );
    return
  end
%
%  Calculate the new components of QUASI.
%
  for i = 1 : dim_num
    quasi(i) = lastq(i) * recipd;
    lastq(i) = bitxor ( lastq(i), v(i,l) );
  end

  seed_save = seed;
  seed = seed + 1;

  return
end
