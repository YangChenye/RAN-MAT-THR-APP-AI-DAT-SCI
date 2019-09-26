function [p] = Wigner(lamda)
% the analytical curve (Wigner's semi-circle law of the limiting distribution)
for l = 1:length(lamda)
    p(l) = (1/(2*pi))*sqrt(4-lamda(l)^2)*(abs(lamda(l))<2);
end

