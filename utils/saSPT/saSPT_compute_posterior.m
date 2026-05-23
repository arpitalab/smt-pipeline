function [naive_occ, post_occ, r0, Rhat, post_mean] = saSPT_compute_posterior(log_L, jumps_per_track, D, loc_error)
% SASPT_COMPUTE_POSTERIOR  Variational-Bayes posterior over (D, sigma) grid.
%
%   [naive_occ, post_occ, r0, Rhat, post_mean] = ...
%        saSPT_compute_posterior(log_L, jumps_per_track, D, loc_error)
%
%   Implements the variational EM update of Bishop PRML Ch. 10 with a
%   symmetric Dirichlet prior (alpha = 1 / (LD * LL)) over the state grid.
%   Iterates up to 100 times; prints when converged.
%
%   Inputs match saSPT_calc_likelihood outputs.
%
%   Outputs
%     naive_occ : LD x LL — prior-only occupancy
%     post_occ  : LD x LL — posterior occupancy (the saSPT result)
%     r0        : LD x LL x N — naive responsibilities per subtrack
%     Rhat      : LD x LL x N — converged responsibilities per subtrack
%     post_mean : LD x LL — Dirichlet posterior mean

    naive_assignment_probabilities = myexp(log_L);
    denom = sum(sum(naive_assignment_probabilities, 2), 1);
    r0 = naive_assignment_probabilities ./ denom;

    LD  = numel(D);
    LL  = numel(loc_error);

    numerator = r0 .* jumps_per_track(ones(LD,1), ones(LL,1), :);
    num1 = squeeze(sum(numerator, 3));
    naive_occ = num1 ./ sum(num1(:));

    L = r0;
    R = L;
    alpha_over_K = 1 / (LD * LL);
    for ii = 1:100
        n = R .* jumps_per_track(ones(LD,1), ones(LL,1), :);
        n = sum(n, 3);
        m = n + alpha_over_K * ones(size(n));
        exp_log_tau = exp(psi(m));
        Rhat = L .* exp_log_tau;
        Rhat = Rhat ./ sum(sum(Rhat, 1), 2);
        diffR = Rhat - R;
        if sum(abs(diffR(:))) / sum(R(:)) < 5e-3
            fprintf('saSPT posterior converged after %d iterations\n', ii);
            break;
        end
        R = Rhat;
    end

    denom = sum(sum(R, 2), 1);
    Rhat  = R ./ denom;
    numerator = Rhat .* jumps_per_track(ones(LD,1), ones(LL,1), :);
    num1 = squeeze(sum(numerator, 3));
    post_occ = num1 ./ sum(num1(:));
    post_mean = n / sum(n(:));
end

function L = myexp(log_L)
    log_L(isnan(log_L)) = -inf;
    log_L(log_L >= 100.0) = 100.0;
    L = exp(log_L - max(max(log_L, [], 1), [], 2));
end
