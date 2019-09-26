function [p] = MarcenkoPastur(lamda,N,M)
% Marcenko-Pastur law of the limiting distribution
alpha = M./N;
for i = 1:length(lamda)
    p(i) = sqrt(4*alpha-(lamda(i)-1-alpha).^2)./(2*pi*lamda(i))*(logical((1-sqrt(alpha)).^2 <= lamda(i)) && logical(lamda(i) <= (1+sqrt(alpha)).^2));
end

end

