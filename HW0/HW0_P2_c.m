KK = [3000, 1];
NN = [100, 2000];
MM = [2*NN(1), 2*NN(2)];

for question = 1:2
    K = KK(question);
    N = NN(question);
    M = MM(question);
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
    ylabel('p(\lambda)');
    xlabel('\lambda');
    title(['Averaged Histograms, K=',num2str(K),', N=',num2str(N),', M=',num2str(M)]);
    hold on
    plot([0.1:0.1:6], MarcenkoPastur([0.1:0.1:6],N,M), '-r', 'linewidth', 2);
    legend('histograms','analytical curve');
    saveas(gcf,['/Users/yangchenye/Downloads/HW0_2_c_',num2str(K),'_',num2str(N),'_',num2str(M),'.png'])
    close;
end





