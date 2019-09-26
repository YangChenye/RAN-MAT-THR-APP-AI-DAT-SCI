clear;clc;
K = 1;
N = 4000;

mu = 0;         % mean value
var = 1/N;      % variance
sd = sqrt(var); % standard deviation 

A = cell(1, K);

% Creat an ensemble of K realizations of N x N Gaussian symmetric matrices
for i = 1:1:K                       % index to an ensemble of K realizations
    for j = 1:1:N                   % index to row of one matrix
        for k = j:1:N               % index to column of one matrix
            A{i}(j,k) = normrnd(mu,sd);
            A{i}(k,j) = A{i}(j,k);  % make the matrix symmetric
        end
    end
end

e = zeros(N,K);

for i = 1:1:K
    e (:, i) = eig(A{i});
end

% Plot a histogram with Normalization set to 'pdf' to produce an estimation of the probability density function.
histogram(e,'Normalization','pdf')
ylabel('Frequency / Width');
xlabel('Eigenvalues');
title('Averaged Histograms, K=1, N=4000');
saveas(gcf,'/Users/yangchenye/Downloads/HW0_1_b_1_4000.png')

