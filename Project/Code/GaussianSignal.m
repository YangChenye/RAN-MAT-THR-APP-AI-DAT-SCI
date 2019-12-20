% Reconstruction problem: 
% Given r pairs [t,p(t)] to find fouriers coefficients

% Given r sample pairs
M = 10;
r = 26;
% Uniform distribution
t = rand(r,1);    % independent variable
% Gaussion distribution
% t = zeros(r,1);
% for i=1:1:r
%     t(i,1) = GaussionInab(0.5,0.15,0,1);
% end

% Gaussian signal
p = normpdf(t,0.5,0.15);


% To reconstruct the Fourier coefficients
F = zeros(2*M+1,r);
for k = -M:1:M
    for q = 1:1:r
        F(k+M+1,q) = (1/sqrt(r))*exp(2*pi*1i*k*t(q,1)); % F(1,1)->k=-M q=1, F(M+1,1)->k=0 q=1, F(2M+1,1)->k=M q=1
    end
end
T = F*F';
b = F*p;
a_hat = T\b;

% preconditioned system
% w = zeros(1,r);
% for q = 2:1:r-1
%     w(q) = (t(q+1)-t(q-1))/2;
% end
% w(1) = (t(2)-t(r)+1)/2;
% w(r) = (1+t(1)-t(r-1))/2;
% W = diag(w);
% T_w = F*W*F';
% b_w = F*W*p;
% a_hat = T_w\b_w;

% Plot Figure 1
figure(1)
x = linspace(0,1)';
y = normpdf(x,0.5,0.15);
plot(x,y,'r-')
hold on
y_hat = FourierSeries(x,a_hat,M)*0.2;
plot(x,y_hat,'b--')
hold on
scatter(t,p,'filled')
xlabel('t')
ylabel('p(t)')
xlim([0,1])
% ylim([0,5])
grid on
hold off
legend('True signal','Estimated signal','Samples')
print(gcf,'-dpng','-r300','Fig_Gaussian_Signal_U.png')
% close(figure(gcf))