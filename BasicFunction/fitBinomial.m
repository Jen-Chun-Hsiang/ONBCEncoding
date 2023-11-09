function p_opt = fitBinomial(data)
    % Define the negative log-likelihood function for binomial distribution
    negLogLikelihood = @(p) -sum(log(binopdf(round(data), 1, p)));

    % Optimization settings
    options = optimset('fmincon');
    options.Display = 'off';
    
    % Use fmincon to find the p that minimizes the negative log-likelihood
    p_initial = mean(data); % Initial guess for p based on data mean
    p_opt = fmincon(negLogLikelihood, p_initial, [], [], [], [], 0, 1, [], options);
end