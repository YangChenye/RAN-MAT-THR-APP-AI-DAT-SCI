function [p] = FourierSeries(t,a,M)
%FourierSeries Fourier series is a periodic function compsed of
%harmonically related sinusoids, combined by a weighter summation
%   A strictly band-limited signal over the interval [0, 1) can be written as the weighted sum of M harmonics in terms of Fourier series
%   Input:  t[<float vector> in [0,1)]:  independent variable
%           a[<float vector> of 2M+1]:  Fourier coefficients
%           M[<int> >= 0]:  M harmonics
%   Output: p:  dependent variable
%   
p = zeros(length(t),1);
for k = -M:1:M
    p = p + a(k+M+1,1).*exp(2*pi*1i*k.*t); % a(1)->k=-M, a(M+1)->k=0, a(2M+1)->k=M
end

end

