function bit = i4_bit_lo0 ( n )

%*****************************************************************************80
%
%% I4_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
%
%  Example:
%
%       N    Binary     BIT
%    ----    --------  ----
%       0           0     1
%       1           1     2
%       2          10     1
%       3          11     3 
%       4         100     1
%       5         101     2
%       6         110     1
%       7         111     4
%       8        1000     1
%       9        1001     2
%      10        1010     1
%      11        1011     3
%      12        1100     1
%      13        1101     2
%      14        1110     1
%      15        1111     5
%      16       10000     1
%      17       10001     2
%    1023  1111111111     1
%    1024 10000000000     1
%    1025 10000000001     1
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 February 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the integer to be measured.
%    N should be nonnegative.
%
%    Output, integer BIT, the position of the low 1 bit.
%
  bit = 0;
  i = floor ( n );

  while ( 1 )

    bit = bit + 1;
    i2 = floor ( i / 2 );

    if ( i == 2 * i2 )
      break;
    end

    i = i2;

  end

  return
end


