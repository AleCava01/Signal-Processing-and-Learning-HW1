%% 1.1
clear;close all; clc;

% parameters (fixed)
sigma_w = 1;
L=10000;
 
% Parameters values
N0s = [0, 5, 10, 20, 50, 100, 200, 400, 800];
rhos = [0, 0.5, 0.9, 0.99];
Ns = linspace(20, 200, 10);

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

        end
    end
end

MSEdB = 10*log10(MSEs); % Ns x rhos x N0s

%% plots
close all;


epsilon_dB = -38;
idx_N = 1; % N=20 samples

figure;
% squeeze extracts a 2D matrix (length(rhos) x length(N0s)); 
% transpose (') needed to have N0s on x-axis
plot(N0s(1:5), squeeze(MSEdB(idx_N, :, 1:5))', '-o', 'LineWidth', 1.5);
grid on;
xlabel('Discarded transient samples (N_0)');
ylabel('MSE [dB]');
yline(epsilon_dB, '--r', 'LineWidth', 1.5, 'DisplayName', ['Threshold \epsilon (' num2str(epsilon_dB) ' dB)']);
title(sprintf('MSE as a function of N_0 (with N = %d)', Ns(idx_N)));
legend(compose('\\rho = %.2f', rhos), 'Location', 'best');


figure;
plot(N0s, squeeze(MSEdB(idx_N, :, :))', '-o', 'LineWidth', 1.5);
grid on;
xlabel('Discarded transient samples (N_0)');
ylabel('MSE [dB]');
yline(epsilon_dB, '--r', 'LineWidth', 1.5, 'DisplayName', ['Threshold \epsilon (' num2str(epsilon_dB) ' dB)']);
title(sprintf('MSE as a function of N_0 (with N = %d)', Ns(idx_N)));
legend(compose('\\rho = %.2f', rhos), 'Location', 'best');


idx_N = 4; % N=60 samples

figure;
plot(N0s, squeeze(MSEdB(idx_N, :, :))', '-o', 'LineWidth', 1.5);
grid on;
xlabel('Discarded transient samples (N_0)');
ylabel('MSE [dB]');
title(sprintf('MSE as a function of N_0 (with N = %d)', Ns(idx_N)));
legend(compose('\\rho = %.2f', rhos), 'Location', 'best');

idx_N = 10; % N=200 samples

figure;
plot(N0s, squeeze(MSEdB(idx_N, :, :))', '-o', 'LineWidth', 1.5);
grid on;
xlabel('Discarded transient samples (N_0)');
ylabel('MSE [dB]');
title(sprintf('MSE as a function of N_0 (with N = %d)', Ns(idx_N)));
legend(compose('\\rho = %.2f', rhos), 'Location', 'best');



%% functions

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
