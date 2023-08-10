load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\DataStimulus\Ai148_AAV-Grm6Cre_121418_1_11.mat');
OUT = DG_OUT;
Grid = OUT.Grids;
%%
% 60 (main) 870 880 882 
pageid = 882;
figure; imagesc(squeeze(Grid(:, :, pageid))'); colorbar
colormap(gray);
axis off
box off

%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure6_NatMovProjection_NoiseCheckboard_pageid%d', SaveFolder, pageid);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);

%%
close all
% run @RFMapping_RF_nonlinearity_batch
% f = 37
figure;  hold on
timwin = [65 125];
x = a(t>timwin(1) & t<=timwin(2));
y = b(t>timwin(1)  & t<=timwin(2));
x = x/max(x);
y = y/max(y);
z = t(t>timwin(1)  & t<=timwin(2));
z = z - min(z);
for i = 1:5
    tids = ((i-1)*1200+1:i*1200);   
    x(tids) = x(tids)-mean(x(tids));
    y(tids) = y(tids)-mean(y(tids));
end
plot(z, x, 'b');
plot(z, y, 'k');
r = corr(x', y');
R = corrchop(a', b');
% title(sprintf('%0.3G %0.3G (%0.3G)', lag(mid)/UpSamplingFz, R, r));
xticks(0:20:60);
xticklabels({'0', '20', '40', '60'});
xlim([0 60]);
yticks(-0.5:0.5:0.5);
yticklabels({'-0.5', '0', '0.5'});
ylim([min([x(:); y(:)])-.1 max([x(:); y(:)])+0.05]);
xlabel('Times (s)');

%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure6_NatMovProjection_AlignedTraces', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);