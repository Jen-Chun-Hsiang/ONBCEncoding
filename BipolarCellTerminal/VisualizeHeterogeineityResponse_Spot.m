%%
[~, sids] = sort(Quality(Quality.^2>0.3), 'ascend');
x = dtest(sids, :);
x = x./max(x, [], 2);
figure; 
imagesc(x); colorbar
title(sprintf('%d-%d %d', ii, f, bctype));
% time
Fz = 9.47;
StmDur = 3; % in second
BasDur = 0.5;
StmDurFrm = round(StmDur*Fz);
BasDurFrm = round(BasDur*Fz);
SpnDurFrm = StmDurFrm+BasDurFrm;
t = round(linspace(-BasDurFrm/Fz, StmDurFrm/Fz, length(1:1:SpnDurFrm))*100)/100;
t = t(4:33);
t = repmat(t(:), 7, 1)';
% size
s = [20 50 100 200 400 600 800];
s = repmat(s, 30, 1);
s = s(:)';

%%
figure; imagesc(cmap); colorbar
%%
switch ii
%     case 7
%         fid = 4;
%         compids = [4 5 22 28];
    case 1
        fid = 1;
        compids = [9 4 6 28];
    case 2
        fid = 1;
        compids = [2 5 18 17];
    case 3
        fid = 1;
        compids = [4 6 25 27];
    case 4
        fid = 1;
        compids = [7 10 23 22];
    case 5
        fid = 1;
        compids = [6 12 22 28];
    case 6
        fid = 1;
        compids = [5 4 25 29];
    case 7
        fid = 2;
        compids = [8 12 27 32];
    case 8
        fid = 1;
        compids = [2 4 22 23];
    case 9
        fid = 1;
        compids = [2 3 5 15];
    case 10
        fid = 1;
        compids = [20 19 2 15];
    case 11
        fid = 1;
        compids = [23 21 5 20];
    case 12
        fid = 1;
        compids = [2 6 11 26];
    case 13
        fid = 1;
        compids = [2 5 27 20];
    case 15
        fid = 1;
        compids = [3 5 31 34];
    case 16
        fid = 1;
        compids = [15 16 2 6];
end
assert(f== fid);
nr = length(compids);
%%
nCent = nan(nROI, 2);
for i = 1:nROI
    nCent(i, :) = getROIloc(cmap == i);
end
tz = mean(PlaneIds);
nCent = [nCent tz*ones(nROI, 1)];
%%
subplot(1, 4, 4); hold on
for i = 1:nr
    text(nCent(compids(i), 1), nCent(compids(i), 2), sprintf('%d', i), 'Color', 'w');
end
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure7_SpotSelecROIs_%d', SaveFolder, ii);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
upscale = 5;
[~, sids] = sort(Quality(Quality.^2>0.3), 'ascend');
insIds = incIds(sids);
assert(all(ismember(compids, insIds)));
% close all
figure; hold on
x = dtest(sids, :);
x = x./max(x, [], 2);
imagesc(imresize(x, upscale, 'nearest')); colorbar
miny = 0;
maxy = length(sids);
for i = 1:6
    plot(30*i*ones(1, 2)*upscale, [miny+1 maxy*upscale], '--k');
end
for i = 1:nr
    loc = find(insIds == compids(i));
    if length(loc) > 1
        loc = loc(1);
    end
    rectangle('Position', [1 ((loc-1)*upscale+1) size(x, 2)*upscale, upscale],...
        'EdgeColor', 'w');
    text(-20, (loc-0.3)*upscale, num2str(i));
end
box off
axis off
title(sprintf('%d. %d-%d-%d (%s)', ii, Day, Cel,  exps(f), BCTypeLabels{find(bctype == BCTypes)}));

%%
%  SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
%  FleNam = sprintf('%sFigure7_SpotSelecROIheatmapResponses_%d', SaveFolder, ii);
%  print('-depsc','-painters','-loose', '-r300',FleNam)
%  saveas(gcf,[FleNam '.png']);
 