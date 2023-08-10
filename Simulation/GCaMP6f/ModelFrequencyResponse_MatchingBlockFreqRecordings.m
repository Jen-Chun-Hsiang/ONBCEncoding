%close all; clear all; clc;
Act = @(MaxAmp, Alpha, Tau_r, Tau_f, t)  MaxAmp * Alpha * (exp(- t/Tau_f) - exp(- t/Tau_r)) /(Tau_f - Tau_r);
% Act = @(MaxAmp, Alpha, Tau_r, Tau_f, t)  MaxAmp * Alpha * (exp(- t/Tau_f) - exp(- t/Tau_r));
SamplingRate = 1000; % times per second
% T = (0:TimeLength-1)/SamplingRate;
T = -1:(1/SamplingRate):3;
DecayTau = 0.1562; %0.142 110%: 0.1562
RiseTau = 0.0495; %0.045 110%: 0.0495
figure; plot(T(T>=0 & T<2 ), Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2)));
GC6f = Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2));
Freqs = [0.5 1 2 4 8 16 32];
nFreq = length(Freqs);

Colors = lines(nFreq);
%%
% Run @BlockFreqAnalysis_BipolarCellTerminal
% and save Ti=T
minV = min(abs(T-Ti'), [], 1);
tids = find(minV < 1e-6);
assert(length(tids) == length(Ti));
%%
figure;
% plot(T/SamplingRate, wave, 'k');
for i = 1:nFreq
    subplot(2, 4, i); hold on
    wave = 0.5*(sin(2*pi*Freqs(i)*T)+1);
    wave(T<0) = 0.5;
    wave(T>2) = 0.5;
%     wave = zeros(size(T));
%     wave(T>1 & T<2) = 1;
    plot(T, wave, 'k');
%     keyboard;
    Response = conv(GC6f, wave,'full')/1000;
    dResp = Response(tids);
    dT = T(tids);
%     dResp = resample(Response, 5, 132);
%     dT =  resample(T, 5, 132);
%     plot(T, Response(1:length(T)), 'Color', Colors(i, :));
    plot(dT, dResp(1:length(dT)), 'Color', Colors(i, :));
    dResp = dResp(1:length(dT));
    dResp = dResp(dT>=-0.25);
    dT = dT(dT>=-0.25);
    if i == 1
        SimT = dT;
        SimResp = nan(nFreq, length(dResp));
    end
    SimResp(i, :) = dResp;
    title(sprintf('%G (Hz)', Freqs(i)));
end
subplot(2, 4, 8);
plot(T(T>=0 & T<2 ), Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2)));
title('GCaMP6 filter (1AP)');
%%
UpFac = 5;
Fz = 37.88;
TI = linspace(SimT(1), SimT(end), length(SimT)*UpFac);
n = length(TI);  
fs = Fz*UpFac;
f = (0:n-1)*(fs/n); 
fl = f(f<Fz/2);
nf = length(fl);
SimFreqPw = nan(nFreq, nf);
for i = 1:nFreq
    YI = interp1(SimT, SimResp(i, :), TI, 'pchip');
    YI = abs(fft(YI)).^2/n;
%     YI = abs(fft(YI))/n;
    SimFreqPw(i, :) = YI(1:nf);
end
%%
figure; imagesc(log10(SimFreqPw)); colorbar
%%
fsets = {2:3, 4, 6:7, 12:13, 23:24, 44:46};
%%
f1pow = nan(nFreq-1, 1);
for i = 1:nFreq-1
    f1pow(i) = mean(SimFreqPw(i, fsets{i}));
end
f1pow = f1pow./max(f1pow);
samplefreqs = [0.5 1 2 4 8 16];
figure; plot(samplefreqs, log10(f1pow(1:6)));
%%
% % save('simf1power_voltage.mat', 'samplefreqs', 'f1pow');
% save('simf1powere.mat', 'samplefreqs', 'f1pow');