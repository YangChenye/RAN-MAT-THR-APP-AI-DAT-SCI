%% Figure 7
clear;clc;
parpool % open parallel pool
% Construct K realizations of matrix T. Each is formed from r sample pairs
M = 40;
beta = [0.25,0.5,0.75];
e = cell(1,length(beta));
e_min = cell(1,length(beta));
linetype = {'-k','-b','-r'};
linetype_min = {'--k','--b','--r'};
parfor index = 1:1:length(beta)
    r = round((2*M+1)/beta(index));

    K = 5000;
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
        
        a = zeros(2*M+1,1);       % Fourier coefficients. a(1)->k=-M, a(M+1)->k=0, a(2M+1)->k=M
        for k = (M+1):1:(2*M+1)   % randomly generate fouriers coefficients of actual signal
            a(k,1) = rand*0.01;
        end
        for k = 1:1:M             % real valued signals
            a(k,1) = conj(a(2*M+2-k,1));
        end
        p = FourierSeries(t,a,M);% sampling
        F = zeros(2*M+1,r);
        for k = -M:1:M
            for q = 1:1:r
                F(k+M+1,q) = (1/sqrt(r))*exp(2*pi*1i*k*t(q,1)); % F(1,1)->k=-M q=1, F(M+1,1)->k=0 q=1, F(2M+1,1)->k=M q=1
            end
        end
        T{index_T} = F*F';
        b = F*p;
        a_hat = T{index_T}\b;
    end

    e{index} = zeros(2*M+1,K); % Because T is a (2M + 1) �� (2M + 1) Hermitian Toeplitz matrix
    e_min{index} = zeros(1,K); % One matrix T has one min eigenvalue

    for index_T = 1:1:K
        e{index}(:, index_T) = eig(T{index_T});
        e_min{index}(1, index_T) = min(e{index}(:, index_T))
    end

end

% Plot a histogram with Normalization set to 'pdf' to produce an estimation of the probability density function.
% histogram(e, 'Normalization','pdf') 
figure(7)
for index = 1:1:length(beta)
    [N_cdf,edges] = histcounts(e{index},'Normalization','cdf');
    [N_cdf_min,edges_min] = histcounts(e_min{index},'Normalization','cdf');
    
    edges = edges(2:end) - (edges(2)-edges(1))/2; % independent variables
    edges_min = edges_min(2:end) - (edges_min(2)-edges_min(1))/2;
    
    plot(log10(edges), log10(N_cdf), linetype{index});
    hold on
    plot(log10(edges_min), log10(N_cdf_min), linetype_min{index});
    hold on
end

ylabel('log_{10}F_{M,\beta}(x), log_{10}F_{M,\beta}^{min}(x)');
xlabel('log_{10}x');
grid on
hold off
legend({'\beta=0.25,\lambda','\beta=0.25,\lambda_m','\beta=0.5,\lambda','\beta=0.5,\lambda_m','\beta=0.75,\lambda','\beta=0.75,\lambda_m'},'Location','southeast')
print(figure(7),'-dpng','-r300','Fig_7.png')
% close(figure(gcf))

delete(gcp('nocreate'))  % close parallel pool
