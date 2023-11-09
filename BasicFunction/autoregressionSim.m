function [Sim, mphi] = autoregressionSim(Sig)
    % input Sig = [repeat time]
    nrep = size(Sig, 1);
    n = size(Sig, 2);
    Sim = zeros(nrep, n);
    phi = nan(nrep, 1);
    for i = 1:nrep
        a = Sig(i, :);
        X = [ones(n-1, 1) a(1:end-1)'];
        y = a(2:end)';
        beta_hat = X \ y;
        phi(i) = beta_hat(2);
        
        epsilon = std(a)*randn(n, 1);
        for t = 2:n
            Sim(i, t) = phi(i) * Sim(i, t-1) + epsilon(t);
        end
    end
    mphi = median(phi, 'omitnan');
end