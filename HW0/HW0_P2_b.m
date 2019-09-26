K = 1;
N = 2000;
M = 2*N;

% Creat an ensemble of K realizations of N x N Wishart matrices
A = cell(1,K);
for i = 1:1:K
    X = randn(N,M);
    A{i} = (1/N)*(X*X');
end

e = zeros(N,K);
for i = 1:1:K
    e (:, i) = eig(A{i});
end

% Plot a histogram with Normalization set to 'pdf' to produce an estimation of the probability density function.
histogram(e, 'Normalization','pdf')
ylabel('Frequency / Width');
xlabel('Eigenvalues');
title(['Averaged Histograms, K=',num2str(K),', N=',num2str(N),', M=',num2str(M)]);
hold on
% plot([-3:0.1:3], Wigner([-3:0.1:3]), '-r', 'linewidth', 2);
% legend('histograms','analytical curve');
saveas(gcf,['/Users/yangchenye/Downloads/HW0_2_b_',num2str(K),'_',num2str(N),'_',num2str(M),'.png'])
% close;