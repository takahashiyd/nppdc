%% verify_carter_mean.m -- check the bias formula used by ipdc_bias_correct
% against direct Monte Carlo of a sample coherence between complex-Gaussian
% vectors with n DOF: E|c_hat|^2 ~= C + (1-C)^2/n (Carter, Knapp & Nuttall
% 1973), and the closed-form inversion round trip.
clear; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src')); cd(root); rng(3)

Mrep = 200000;
for n = [26 56 116]
    fprintf('\nn = %d complex DOF (m = %d)\n', n, n-1);
    fprintf('%-8s %-10s %-10s %-10s\n', 'C', 'MC mean', 'Carter', 'MC se');
    for C = [0 0.1 0.3 0.534 0.8 0.95]
        % x, y complex Gaussian with coherence C: y = sqrt(C) x + sqrt(1-C) e
        x = (randn(n, Mrep) + 1i*randn(n, Mrep)) / sqrt(2);
        e = (randn(n, Mrep) + 1i*randn(n, Mrep)) / sqrt(2);
        y = sqrt(C) * x + sqrt(1 - C) * e;
        c = abs(sum(conj(x) .* y, 1)).^2 ./ (sum(abs(x).^2, 1) .* sum(abs(y).^2, 1));
        fprintf('%-8.3f %-10.5f %-10.5f %-10.5f\n', C, mean(c), C + (1 - C)^2 / n, std(c) / sqrt(Mrep));
    end
    fprintf('round trip (closed-form inversion of the MC mean):');
    for C = [0.1 0.3 0.534 0.8]
        x = (randn(n, Mrep) + 1i*randn(n, Mrep)) / sqrt(2);
        e = (randn(n, Mrep) + 1i*randn(n, Mrep)) / sqrt(2);
        y = sqrt(C) * x + sqrt(1 - C) * e;
        c = abs(sum(conj(x) .* y, 1)).^2 ./ (sum(abs(x).^2, 1) .* sum(abs(y).^2, 1));
        [~, I2j] = ipdc_bias_correct(sqrt(mean(c)), n, 1);   % dof - N + 1 = n
        fprintf('  C=%.3f -> %.4f', C, I2j);
    end
    fprintf('\n');
end
