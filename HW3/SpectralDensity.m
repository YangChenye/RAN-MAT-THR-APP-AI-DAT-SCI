function [p] = SpectralDensity(lambda,k)
%the spectral density according to the McKay Law for regular sparse matrices
p = (k.*sqrt(4.*(k-1)-lambda.^2))./(2*pi.*(k^2-lambda.^2));
end

