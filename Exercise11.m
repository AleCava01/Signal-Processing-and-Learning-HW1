%% 1.1
clear;close all; clc;

% parameters (fixed)
sigma_w = 1;
L=10000 %to set higher
 
% Parameters values
N0s = [0, 5, 10, 20, 50, 100, 200];
rhos = [0, 0.5, 0.9, 0.99];
Ns = linspace(20, 200, 10);

% values to search
minN0 = 0;
minN = 0;
minRho = 0;
minMSE = 1000;

% storage for MSEs
MSEs = zeros(length(Ns), length(rhos), length(N0s));

% Monte Carlo loop
for i = 1:length(Ns)
    N = Ns(i);
    for j = 1:length(rhos)
        rho = rhos(j);
        for k = 1:length(N0s)
            N0 = N0s(k);

            % Create covariance matrix
            C = createMatrix(N,N, rho,sigma_w);

            % L loop
            A = zeros(N); % N x N matrix
            for q = 1:L
                w = sigma_w * randn(N+N0,1);
                a = sqrt(1 - rho^2)*sigma_w;
                b = [1, -rho];
                
                y = filter(a, b, w);
                last_y = y(end-N+1:end);
                A = A + last_y*transpose(last_y);
            end
            sample_cov = (1/L).*A;
            
            MSEs(i,j,k) = mean((C - sample_cov).^2, 'all');

            if (MSEs(i,j,k)<minMSE)
                minMSE = MSEs(i,j,k);
                minRho = rhos(j);
                minN = Ns(i);
                minN0 = N0s(k);
            end
        end
    end
end

minMSE
minN
minN0
minRho

% functions

function y = covarianceFunction(sigma2w, rho, k)
    y = sigma2w*rho^k;
end

function A = createMatrix(n,m, rho, sigma2w)
    A = zeros(n,m);

    for i = 1:n
        for j = 1:m
            A(i,j) = covarianceFunction(sigma2w, rho, abs(i-j));
        end
    end
end
