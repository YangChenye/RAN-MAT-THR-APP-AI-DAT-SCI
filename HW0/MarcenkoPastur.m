function [p] = MarcenkoPastur(lamda,N,M)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
alpha = M./N;
p = sqrt(4*alpha-(lamda-1-alpha).^2)./(2*pi*lamda)*(((1-sqrt(alpha)).^2 <= lamda) && (lamda <= (1+sqrt(alpha)).^2));
end

