function [y,sd] = hetsphere(x)
%HETSPHERE Sphere function with added (known) heteroskedastic noise.

f = sum(x.^2,2);            % Noiseless function values
sd = 1 + sqrt(f);         % Standard deviation of the noise
y = f + sd.*randn(size(f)); % Noisy values
