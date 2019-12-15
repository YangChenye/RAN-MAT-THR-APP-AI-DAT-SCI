%% Figure 5
clear;clc;
parpool % open parallel pool
% Construct K realizations of matrix T. Each is formed from r sample pairs
M = 200;
beta = [0.1:0.1:0.8];
e = cell(1,length(beta));
linetype = {'-k','-.k','--k'};
parfor index = 1:1:length(beta)
    r = round((2*M+1)/beta(index));

    K = 1000;
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

    e{index} = zeros(2*M+1,K); % Because T is a (2M + 1) ¡Á (2M + 1) Hermitian Toeplitz matrix

    for index_T = 1:1:K
        e{index}(:, index_T) = eig(T{index_T});
    end

end

% Plot a histogram with Normalization set to 'pdf' to produce an estimation of the probability density function.
% histogram(e, 'Normalization','pdf') 
figure(5)
for index = 1:1:length(beta)
    [N_cdf,edges] = histcounts(e{index},'Normalization','cdf');
    %N_cdf dependent variables
    edges = edges(2:end) - (edges(2)-edges(1))/2; % independent variables
    
    N_cdf_log = log10(N_cdf);
    edges_log = log10(edges);
    
    [M,I] = min(abs(N_cdf_log + 2)); % find the tangent point
    x0 = edges_log(I);
    y0 = N_cdf_log(I);
    
%     k = beta(index); % slope
%     b = y0 - k * x0; % intercept
%     ff = @(x) k * x + b; % euqation
%     yy = ff(log10(edges));
%     l1 = plot(log10(edges), yy, '-k');
    
    % use secant line as tangent line
    if ((I-1)>=1 && (I+2)<=length(edges_log))
        x1 = edges_log(I-1); 
        y1 = N_cdf_log(I-1);
        x2 = edges_log(I+1);
        y2 = N_cdf_log(I+1);
        k = (y2-y1)/(x2-x1); % slope
        b = y0 - k * x0;     % intercept
        ff = @(x) k * x + b; % euqation
    elseif ((I-1)==0)
        x2 = edges_log(I+1);
        y2 = N_cdf_log(I+1);
        k = (y2-y0)/(x2-x0); % slope
        b = y0 - k * x0;     % intercept
        ff = @(x) k * x + b; % euqation
    elseif ((I+2)>length(edges_log))
        x1 = edges_log(I-1); 
        y1 = N_cdf_log(I-1);
        k = (y0-y1)/(x0-x1); % slope
        b = y0 - k * x0;     % intercept
        ff = @(x) k * x + b; % euqation
    end
    
    yy = ff(edges_log);
    l1 = plot(edges_log, yy, '-k');
    hold on
    l2 = plot(edges_log, N_cdf_log, '--k');
    hold on
end

ylabel('log_{10}F_{M,\beta}(x)');
xlabel('log_{10}x');
grid on
hold off
legend([l1 l2],{'Tangent at -2','Simulation'})
print(figure(5),'-dpng','-r300','Fig_5.png')
% close(figure(gcf))

delete(gcp('nocreate'))  % close parallel pool
