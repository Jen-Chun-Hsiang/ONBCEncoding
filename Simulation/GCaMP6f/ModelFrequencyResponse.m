close all; clear all; clc;
Act = @(MaxAmp, Alpha, Tau_r, Tau_f, t)  MaxAmp * Alpha * (exp(- t/Tau_f) - exp(- t/Tau_r)) /(Tau_f - Tau_r);
% Act = @(MaxAmp, Alpha, Tau_r, Tau_f, t)  MaxAmp * Alpha * (exp(- t/Tau_f) - exp(- t/Tau_r));
SamplingRate = 1000; % times per second
% T = (0:TimeLength-1)/SamplingRate;
T = -1:(1/SamplingRate):3;
DecayTau = 0.142;
RiseTau = 0.045;
figure; plot(T(T>=0 & T<2 ), Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2)));
GC6f = Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2));
Freqs = [0.5 1 2 4 8 16 32];
nFreq = length(Freqs);

Colors = lines(nFreq);
%% after run %ModelFrequencyResponse_MatchingBlockFreqRecordings
 load('simf1powere.mat', 'samplefreqs', 'f1pow');
%%
figure;
% plot(T/SamplingRate, wave, 'k');
subplot(2, 4, 1);
plot(T(T>=0 & T<2 ), Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2)), 'Color', [50 168 82]/255);
box off
ylabel('Amplitude (arbi.)');
xlabel('Times (s)');
title('GCaMP6 filter (1AP)');
for i = 1:nFreq-1
    subplot(2, 4, i+1); hold on
    wave = 0.5*(sin(2*pi*Freqs(i)*T)+1);
    wave(T<0) = 0.5;
    wave(T>2) = 0.5;
    
    % step function
%     wave = wave>0.5;
%     xlim([-1 0]);

    plot(T, wave, 'Color', 0.5*ones(1, 3));
    Response = conv(GC6f, wave,'full')/1000;
    plot(T, Response(1:length(T)), 'Color', [50 168 82]/255);
    title(sprintf('%G (Hz)', Freqs(i)));
    xlabel('Times (s)');
    ylabel('Amplitude');
    yticks(0:0.5:1);
    yticklabels({'0', '0.5', '1'});
end
subplot(2, 4, 8);
plot(samplefreqs, f1pow, 'k');
set(gca, 'XScale', 'log');
yticks(0:0.5:1);
yticklabels({'0', '0.5', '1'});
box off
ylabel('Amplitude (norm.)')
xlabel('Times (s)');
title('f1 power reduction');
 %%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sSupplement_GCaMP6fFrequencyResponseSimulation', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);

%%
UpFac = 5;
Fz = 200;
T = -1:(1*UpFac/SamplingRate):3;
TI = linspace(T(1), T(end), length(T)*UpFac);
n = length(TI);  
fs = Fz*UpFac;
f = (0:n-1)*(fs/n); 
fl = f(f<Fz/2);
nf = length(fl);
YI = interp1(T, randn(size(T)), TI, 'pchip');
% YI = abs(fft(YI)).^2/n;
YI = abs(fft(YI))/n;
FP = YI(1:nf);
