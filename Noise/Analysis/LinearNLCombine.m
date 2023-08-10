function [Optm, Eval] = LinearNLCombine(L, NL, Sig)
L = L(:)/max(L);
NL = NL(:)/max(NL);
Sig = Sig(:);
Len = length(L);
nDiv = ceil(Len/1200);% 60 frames * Sampling 100 Hz / 5 Hz
% CostF = @(w) 1/corrchop(Sig, L+NL*w, nDiv).^2;
CostF = @(w) 1/abs(corrchop(Sig, (1-w)*L+w*NL, nDiv));
Optm = fmincon(CostF, 0.1, [], [], [], [],...
    0,...
    1);
Eval = corrchop(Sig, L+NL*Optm, nDiv);
end