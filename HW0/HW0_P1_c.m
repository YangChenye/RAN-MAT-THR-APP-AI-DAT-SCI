% a
clear;clc;
KK1 = [10000];
NN1 = [5, 20];

for kindex = 1:length(KK1)
    for nindex = 1:length(NN1)
        
        K = KK1(kindex);
        N = NN1(nindex);

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
        histogram(e, 'Normalization','pdf')
        ylabel('p(\lambda)');
        xlabel('\lambda');
        title(['Averaged Histograms, K=',num2str(K),', N=',num2str(N)]);
        hold on
        plot([-3:0.1:3], Wigner([-3:0.1:3]), '-r', 'linewidth', 2);
        legend('histograms','analytical curve');
        saveas(gcf,['/Users/yangchenye/Downloads/HW0_1_c_',num2str(K),'_',num2str(N),'.png'])
        close;
    end
end


% b
clear;clc;
KK2 = [1];
NN2 = [10, 100, 1000, 4000];

for kindex = 1:length(KK2)
    for nindex = 1:length(NN2)
        
        K = KK2(kindex);
        N = NN2(nindex);

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
        histogram(e, 'Normalization','pdf') 
        ylabel('p(\lambda)');
        xlabel('\lambda');
        title(['Averaged Histograms, K=',num2str(K),', N=',num2str(N)]);
        hold on
        plot([-3:0.1:3], Wigner([-3:0.1:3]), '-r', 'linewidth', 2);
        legend('histograms','analytical curve');
        saveas(gcf,['/Users/yangchenye/Downloads/HW0_1_c_',num2str(K),'_',num2str(N),'.png'])
        close;
    end
end



