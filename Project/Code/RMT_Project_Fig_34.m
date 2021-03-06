%% Figure 3
clear;clc;
parpool % open parallel pool
% Construct K realizations of matrix T. Each is formed from r sample pairs
M = [1,4,10,90];
e = cell(1,length(M));
linetype = {'*-','s-','o-','.-'};
parfor index_M = 1:1:length(M) % each loop is independent, parallel programming
    beta = 0.25;
    r = (2*M(index_M)+1)/beta;

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
        
        a = zeros(2*M(index_M)+1,1);       % Fourier coefficients. a(1)->k=-M, a(M+1)->k=0, a(2M+1)->k=M
        for k = (M(index_M)+1):1:(2*M(index_M)+1)   % randomly generate fouriers coefficients of actual signal
            a(k,1) = rand*0.01;
        end
        for k = 1:1:M(index_M)             % real valued signals
            a(k,1) = conj(a(2*M(index_M)+2-k,1));
        end
        p = FourierSeries(t,a,M(index_M));% sampling
        F = zeros(2*M(index_M)+1,r);
        for k = -M(index_M):1:M(index_M)
            for q = 1:1:r
                F(k+M(index_M)+1,q) = (1/sqrt(r))*exp(2*pi*1i*k*t(q,1)); % F(1,1)->k=-M q=1, F(M+1,1)->k=0 q=1, F(2M+1,1)->k=M q=1
            end
        end
        T{index_T} = F*F';
        b = F*p;
        a_hat = T{index_T}\b;
    end

    e{index_M} = zeros(2*M(index_M)+1,K); % Because T is a (2M + 1) �� (2M + 1) Hermitian Toeplitz matrix

    for index_T = 1:1:K
        e{index_M}(:, index_T) = eig(T{index_T});
    end

end

% Plot a histogram with Normalization set to 'pdf' to produce an estimation of the probability density function.
% histogram(e, 'Normalization','pdf') 
figure(3)
for index_M = 1:1:length(M)
    [N,edges] = histcounts(e{index_M}, 30, 'Normalization','pdf');
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N, linetype{index_M});
    hold on
end

ylabel('f_{M,\beta}(x)');
xlabel('x');
grid on
hold off
legend('r=12, M=1','r=36, M=4','r=84, M=10','r=724, M=90')
print(figure(3),'-dpng','-r300','Fig_3_c_1.png')
% close(figure(gcf))


%% Figure 4
clear;clc;
% Construct K realizations of matrix T. Each is formed from r sample pairs
M = [94,90,87,85,104];
beta = [0.15,0.25,0.35,0.45,0.55];
e = cell(1,length(M));
linetype = {'*-','s-','d-','o-','.-'};
parfor index_M = 1:1:length(M) % each loop is independent, parallel programming
    r = round((2*M(index_M)+1)/beta(index_M));

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
        
        a = zeros(2*M(index_M)+1,1);       % Fourier coefficients. a(1)->k=-M, a(M+1)->k=0, a(2M+1)->k=M
        for k = (M(index_M)+1):1:(2*M(index_M)+1)   % randomly generate fouriers coefficients of actual signal
            a(k,1) = rand*0.01;
        end
        for k = 1:1:M(index_M)             % real valued signals
            a(k,1) = conj(a(2*M(index_M)+2-k,1));
        end
        p = FourierSeries(t,a,M(index_M));% sampling
        F = zeros(2*M(index_M)+1,r);
        for k = -M(index_M):1:M(index_M)
            for q = 1:1:r
                F(k+M(index_M)+1,q) = (1/sqrt(r))*exp(2*pi*1i*k*t(q,1)); % F(1,1)->k=-M q=1, F(M+1,1)->k=0 q=1, F(2M+1,1)->k=M q=1
            end
        end
        T{index_T} = F*F';
        b = F*p;
        a_hat = T{index_T}\b;
    end

    e{index_M} = zeros(2*M(index_M)+1,K); % Because T is a (2M + 1) �� (2M + 1) Hermitian Toeplitz matrix

    for index_T = 1:1:K
        e{index_M}(:, index_T) = eig(T{index_T});
    end

end

% Plot a histogram with Normalization set to 'pdf' to produce an estimation of the probability density function.
% histogram(e, 'Normalization','pdf') 
figure(4)
for index_M = 1:1:length(M)
    [N,edges] = histcounts(e{index_M}, 30, 'Normalization','pdf');
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N, linetype{index_M});
    hold on
end

ylabel('f_{M,\beta}(x)');
xlabel('x');
grid on
hold off
legend('\beta=0.15, r=1260, M=94','\beta=0.25, r=724, M=90','\beta=0.35, r=500, M=87','\beta=0.45, r=380, M=85','\beta=0.55, r=380, M=104')
print(figure(4),'-dpng','-r300','Fig_4_c_1.png')
% close(figure(gcf))

delete(gcp('nocreate'))  % close parallel pool

