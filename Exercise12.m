clear; close all; clc;
rng(0);

%% Parameters
Bs          = [pi/8 pi/2 .9*pi].';
Ks          = [5 10 20 50].';
sigma_alpha = 1;
N           = 2^12;     % samples in s[n]
Nfft        = 2^12;

omega = (-Nfft/2:Nfft/2-1)' * 2*pi/Nfft;

%% Point a)
figure;
tiledlayout(length(Ks), length(Bs));

for i = 1:length(Ks)
    for j = 1:length(Bs)
        B = Bs(j);
        K = Ks(i);

        % PHASE a) Signal generation
        h = generate_filter(B, K);
        s = generate_signal(h, K, N, sigma_alpha);

        % PHASE b) PSD
        S_estimated = estimate_psd(s, N, Nfft);
        S_ideal     = sigma_alpha^2 * (abs(omega) <= B);
        S_theory    = sigma_alpha^2 * fftshift(abs(fft(h, Nfft)).^2);

        % PLOT
        nexttile;
        if(i==1 && j==1)
            plot_psd(omega, S_estimated, S_ideal, S_theory, B, K, true);
        else
            plot_psd(omega, S_estimated, S_ideal, S_theory, B, K, false);
        end
    end
end
sgtitle('PSD estimation for different B and K');

%% Point b)

Kgrid = 1:1000;
Kopt  = zeros(size(Bs));
Ematrix = zeros(length(Bs), length(Kgrid));

for i = 1:length(Bs)
    B = Bs(i);
    for j = 1:length(Kgrid)
        K = Kgrid(j);
        % generate filter
        h = generate_filter(B, K);

        % compute e
        e = metric_oob(h, B, Nfft, omega);
        Ematrix(i,j) = e;

        if(Kopt(i) == 0 && e<-40)
            Kopt(i) = K;
            fprintf('Found minimum K for B = %.4f, Kopt = %d\n', B, K);
        end
    end
end

figure;
hold on; grid on;

for i = 1:length(Bs)
    semilogx(Kgrid, Ematrix(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('B = %.2f', Bs(i)));
end
yline(-40, '--k', 'LineWidth', 1.2, 'DisplayName', 'Threshold (-40 dB)');

xlabel('K (filter half-width)');
ylabel('\epsilon(K, B) [dB]');
title('\epsilon(K, B) vs. K for different values of B');
legend('Location', 'best');
hold off;

%% Point b extra) 

Kgrid = 1:3000;
Kopt  = zeros(size(Bs));
Ematrix = zeros(length(Bs), length(Kgrid));

for i = 1:length(Bs)
    B = Bs(i);
    for j = 1:length(Kgrid)
        K = Kgrid(j);
        % generate filter
        
        h = generate_filter(B, K);
        s = generate_signal(h, K, N, sigma_alpha);

        % compute e
        e = metric_oob(s, B, Nfft, omega);
        Ematrix(i,j) = e;

        if(Kopt(i) == 0 && e<-40)
            Kopt(i) = K;
            fprintf('Found minimum K for B = %.4f, Kopt = %d\n', B, K);
        end
    end
end

figure;
hold on; grid on;

for i = 1:length(Bs)
    semilogx(Kgrid, Ematrix(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('B = %.2f', Bs(i)));
end
yline(-40, '--k', 'LineWidth', 1.2, 'DisplayName', 'Threshold (-40 dB)');

xlabel('K (filter half-width)');
ylabel('\epsilon(K, B) [dB]');
title('\epsilon(K, B) vs. K for different values of B');
legend('Location', 'best');
hold off;


%% Functions

function e = metric_oob(h, B, Nfft, omega)
    P   = fftshift(abs(fft(h, Nfft)).^2);
    inb = abs(omega) <= B;
    e   = 10*log10(sum(P(~inb)) / sum(P(inb)));
end


function h = generate_filter(B, K)
% Impulive response of the low-pass filter truncated to |k| <= K
    h = zeros(2*K+1, 1);
    for idx = 1:(2*K+1)
        k_val = idx - K - 1;   % from -K to +K
        if k_val == 0
            h(idx) = B / pi;
        else
            h(idx) = sin(B * k_val) / (pi * k_val);
        end
    end
end

function s = generate_signal(h, K, N, sigma_alpha)
% Filter the white gaussian noise removing the transition effect
    alpha  = sigma_alpha * randn(N + 2*K, 1);
    s_full = conv(h, alpha);
    s      = s_full(2*K+1 : 2*K+N);
end

function S_estimated = estimate_psd(s, N, Nfft)
% Periodogram, centered in omega = 0
    X           = fft(s, Nfft);
    S_estimated = (1/N) * abs(X).^2;
    S_estimated = fftshift(S_estimated);
end

function plot_psd(omega, S_estimated, S_ideal, S_theory, B, K, is_first)
    hold on;
    plot(omega, 10*log10(max(S_estimated, 1e-12)));
    plot(omega, 10*log10(max(S_ideal,     1e-12)), '--');
    plot(omega, 10*log10(max(S_theory,    1e-12)), 'r', 'LineWidth', 1);
    if(is_first)
        legend('Estimated PSD', 'Ideal PSD', 'Truncated filter (theory)','Location','southwest');
    end
    xlabel('\omega (rad/sample)');
    ylabel('PSD (dB)');
    title(sprintf('B = %.2f, K = %d', B, K));
    grid on;
end
