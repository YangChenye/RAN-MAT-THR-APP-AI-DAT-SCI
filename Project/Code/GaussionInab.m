function [G] = GaussionInab(mu,sigma,a,b)
%generate a Gaussion distribution with all value in interval [a,b]
%

G = normrnd(mu,sigma);
if (G<a || G>b)
    G = GaussionInab(mu,sigma,a,b);

end

