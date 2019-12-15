%% Figure 8
clear;clc;
parpool % open parallel pool
% Construct K realizations of matrix T. Each is formed from r sample pairs
M = [10,20,40];
beta = 0.25;
e = cell(1,length(M));
e_min = cell(1,length(M));
e_max = cell(1,length(M));
e_k = cell(1,length(M));
linetype_k = {'-k','-b','-r'};
linetype_min = {'--k','--b','--r'};
parfor index = 1:1:length(M)
    r = round((2*M(index)+1)/beta);

    K = 10000;
    T = cell(1, K);

    % Creat an ensemble of K realizations of matrices
    for index_T = 1:1:K                       % index to an ensemble of K realizations
        % the sampling points tq are i.i.d. random variables with uniform distribution U[0,1)
        % Uniform distribution
        % t = rand(r,1);    % independent variable
        % Gaussion distribution
        t = zeros(r,1);
        for i=1:1:r
            t(i,1) = GaussionInab(0.5,0.15,0,1);
        end
        
        a = zeros(2*M(index)+1,1);       % Fourier coefficients. a(1)->k=-M, a(M+1)->k=0, a(2M+1)->k=M
        for k = (M(index)+1):1:(2*M(index)+1)   % randomly generate fouriers coefficients of actual signal
            a(k,1) = rand*0.01;
        end
        for k = 1:1:M(index)             % real valued signals
            a(k,1) = conj(a(2*M(index)+2-k,1));
        end
        p = FourierSeries(t,a,M(index));% sampling
        F = zeros(2*M(index)+1,r);
        for k = -M(index):1:M(index)
            for q = 1:1:r
                F(k+M(index)+1,q) = (1/sqrt(r))*exp(2*pi*1i*k*t(q,1)); % F(1,1)->k=-M q=1, F(M+1,1)->k=0 q=1, F(2M+1,1)->k=M q=1
            end
        end
        T{index_T} = F*F';
        b = F*p;
        a_hat = T{index_T}\b;
    end

    e{index} = zeros(2*M(index)+1,K); % Because T is a (2M + 1) ¡Á (2M + 1) Hermitian Toeplitz matrix
    e_min{index} = zeros(1,K); % One matrix T has one min eigenvalue
    e_max{index} = zeros(1,K); % One matrix T has one min eigenvalue
    e_k{index} = zeros(1,K); % One matrix T has one min eigenvalue

    for index_T = 1:1:K
        e{index}(:, index_T) = eig(T{index_T});
        e_min{index}(1, index_T) = min(e{index}(:, index_T));
        e_max{index}(1, index_T) = max(e{index}(:, index_T));
        e_k{index}(1, index_T) = e_max{index}(1, index_T)/e_min{index}(1, index_T);
    end

end

% Plot a histogram with Normalization set to 'pdf' to produce an estimation of the probability density function.
% histogram(e, 'Normalization','pdf') 
figure(8)
for index = 1:1:length(M)
    % [N_cdf,edges] = histcounts(e{index},'Normalization','pdf');
    [N_cdf_min,edges_min] = histcounts(e_min{index},'Normalization','pdf');
    [N_cdf_k,edges_k] = histcounts(e_k{index},'Normalization','pdf');
    
    % edges = edges(2:end) - (edges(2)-edges(1))/2; % independent variables
    edges_min = edges_min(2:end) - (edges_min(2)-edges_min(1))/2;
    edges_k = edges_k(2:end) - (edges_k(2)-edges_k(1))/2;
    
    plot(log10(edges_k), log10(N_cdf_k), linetype_k{index});
    hold on
    plot(log10(edges_min), log10(N_cdf_min), linetype_min{index});
    hold on
end

ylabel('log_{10}f_{M,\beta}^{\kappa}(x), log_{10}f_{M,\beta}^{min}(x)');
xlabel('log_{10}x');
grid on
hold off
legend({'\kappa, M=10', '\lambda_{min}, M=10','\kappa, M=20', '\lambda_{min}, M=20','\kappa, M=40', '\lambda_{min}, M=40'},'Location','northeast')
print(figure(8),'-dpng','-r300','Fig_8.png')
% close(figure(gcf))

delete(gcp('nocreate'))  % close parallel pool
