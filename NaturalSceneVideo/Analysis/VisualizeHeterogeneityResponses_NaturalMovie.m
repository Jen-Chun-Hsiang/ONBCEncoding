% run @NaturalMovieHeterogeneityEncodingMapping first
% when it stops, go to @GetPairedPathDistance
% first time:
%%
figure; imagesc(cmap); colorbar;
%%
% Representative for each BC types for the heterogeity responses across
% terminals
% (1). Show the align functional images on top of the structure and note the
%       ROI number
% (2). Plot the responses of each select ROI for its repeated responses to
%       the natural movies
% final selection
% BC5o:15 , BC5i:4, BC5t:5, XBC:16 , BC6:9 , BC7:12 , BC8/9:13
switch ii
    case 1
        fid = 1;
        prid = 1;
        compids = [5 6 20 17];
    case 2
        fid = 1;
        prid = 1;
        compids = [12 14 2 13];
    case 3
        fid = 2;
        prid = 1;
        compids = [3 5 17 19];
    case 4
        fid = 1;
        prid = 1;
        compids = [18 16 7];% 25 19
    case 5
        fid = 1;
        prid = 1;
        compids = [7 5 18]; % 16
    case 6
        fid = 2;
        prid = 1;
        compids = [14 3 10 11];
    case 7
        fid = 1;
        prid = 1;
        compids = [28 29 17 7];
    case 16
        fid = 1;
        prid = 1;
        compids = [12 10 5];%5 3
    case 8 % very bad
        fid = 1;
        prid = 1;
        compids = [3 5 11 13];
    case 9
        fid = 1;
        prid = 1;
        compids = [3 4 8]; % 9
    case 10 % kind of bad
        fid = 1;
        prid = 1;
        compids = [5 3 11 2];
    case 11
        fid = 1;
        prid = 1;
        compids = [5 6 17 19];
    case 12
        fid = 1;
        prid = 1;
        compids = [5 14 19];%22 19
    case 13
        fid = 1;
        prid = 1;
        compids = [4 3 22];% 18
    case 15
        fid = 1;
        prid = 1;
        compids = [35 34 11];% 7
end
assert(f== fid && p== prid);

nr = length(compids);
%%
nCent = nan(nROI, 2);
for i = 1:nROI
    nCent(i, :) = getROIloc(cmap == i);
end
tz = mean(PlaneIds);
nCent = [nCent tz*ones(nROI, 1)];
%%
% run first section of @GetPairedPathDistance
%%
subplot(1, 4, 4); hold on
for i = 1:nr
    text(nCent(compids(i), 1), nCent(compids(i), 2), sprintf('%d', i), 'Color', 'w');
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure5_NaturalMovieSelecROIs_%d', SaveFolder, ii);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%% select clip display
% dispids = [4 7 9 10];
dispids = 1:11;
ndisp = length(dispids);
figure; hold on

Fz = 100;
t = (0:size(dtest, 2)-1)/Fz;
cliplines = cumsum(ClipLens-90)/Fz;
cliplinecenter = mean([[0; cliplines(1:(end-1))] cliplines], 2);
nClip = length(cliplines);

avgtrace = mean(targetD1, 1);
avgtrace = avgtrace - min(avgtrace);
avgtrace1 = avgtrace/max(avgtrace);

avgtrace = mean(targetD2, 1);
avgtrace = avgtrace - min(avgtrace);
avgtrace2 = avgtrace/max(avgtrace);
% plot(t, avgtrace+1, 'Color', [0 180 171]/255, 'LineWidth', 2);
atlen = nan(ndisp, 3);
x = [];
xc = [];
y1 = [];
y2 = [];
clear qids
for i = 1:nr
    a = find(cTab(:, 2)== compids(i));
    qids(i) =a(1);
end
for j = 1:ndisp
    if dispids(j) == 1
        tids = 1:round(cliplines(dispids(j))*Fz);
    else
        tids = round(cliplines(dispids(j)-1)*Fz):round(cliplines(dispids(j))*Fz);
    end
    if j == 1
        atlen(j, :) = [length(tids)/Fz 0.5*length(tids)/Fz length(tids)];
        t = (0:atlen(j, 3)-1)/Fz;
    else
        atlen(j, :) = [atlen(j-1, 1)+length(tids)/Fz atlen(j-1, 1)+0.5*length(tids)/Fz atlen(j-1, 3)+length(tids)];
        t = (atlen(j-1, 3):atlen(j, 3)-1)/Fz;
    end
    text(atlen(j, 2), 1.5, sprintf('%d', dispids(j)), 'Color', 'k');
    if j < ndisp
        plot(atlen(j, 1)*ones(1, 2), [-nr*2+1; 1], '--k');
    end
    x = [x; ((0:length(t)-1)/Fz)'];
    xc = [xc; j*ones(length(t), 1)];
    y1 = [y1; targetD1(qids, tids)'];
    y2 = [y2; targetD2(qids, tids)'];
    for i = 1:nr
        fcompids = find(cTab(:, 2)== compids(i));
        plot(t, avgtrace1(tids)-(i-1), 'Color', 0.5*ones(1, 3));
        plot(t, mean(targetD1(fcompids, tids), 1)-(i-1), 'Color', [2 105 164]/255);
        plot(t, avgtrace2(tids)-(i-1)-nr, 'Color', 0.5*ones(1, 3));
        plot(t, mean(targetD2(fcompids, tids), 1)-(i-1)-nr, 'Color', [0 180 171]/255);
        if j == 1
            text(-2, -(i-1)+0.2, sprintf('%d', i), 'Color', 'k');
            text(-2, -(i-1)-nr+0.2, sprintf('%d', i), 'Color', 'k');
        end
        
    end
end
xlabel('Time (s)');
box off
axis off
sgtitle(sprintf('%d. %d-%d (%d, %d) (%s)', ii, Day, Cel,  expids(1), expids(2), BCTypeLabels{DTab(cids(1), 4)}));
%%
 SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
 FleNam = sprintf('%sFigure5_NaturalMovieSelecROIandClipsResponses_%d', SaveFolder, ii);
 print('-depsc','-painters','-loose', '-r300',FleNam)
 saveas(gcf,[FleNam '.png']);
 
 %% overlapping with average
figure;
hold on

Fz = 100;
t = (0:size(dtest, 2)-1)/Fz;
cliplines = cumsum(ClipLens-90)/Fz;
cliplinecenter = mean([[0; cliplines(1:(end-1))] cliplines], 2);
nClip = length(cliplines);

avgtrace = mean(targetD1, 1);
avgtrace = avgtrace - min(avgtrace);
avgtrace1 = avgtrace/max(avgtrace);

avgtrace = mean(targetD2, 1);
avgtrace = avgtrace - min(avgtrace);
avgtrace2 = avgtrace/max(avgtrace);
% plot(t, avgtrace+1, 'Color', [0 180 171]/255, 'LineWidth', 2);

for i = 1:nr
    fcompids = find(cTab(:, 2)== compids(i));
    plot(t, avgtrace1-(i-1)*2, 'Color', 0.5*ones(1, 3));
    plot(t, mean(targetD1(fcompids, :), 1)-(i-1)*2, 'Color', [2 105 164]/255);
    plot(t, avgtrace2-(i-1)*2-1, 'Color', 0.5*ones(1, 3));
    plot(t, mean(targetD2(fcompids, :), 1)-(i-1)*2-1, 'Color', [0 180 171]/255);
    
    text(-2, -i*2+2, sprintf('%d', i), 'Color', 'k');
    
end
% text(-2, 2, 'average', 'Color', 'k');
for i = 1:nClip
    text(cliplinecenter(i), 1.5, sprintf('%d', i), 'Color', 'k');
end
plot(cliplines(1:end-1)'.*ones(2, nClip-1), [-nr*2+1; 1].*ones(2, nClip-1), '--k');

xlabel('Time (s)');
box off
axis off
sgtitle(sprintf('%d. %d-%d (%d, %d) (%s)', ii, Day, Cel,  expids(1), expids(2), BCTypeLabels{DTab(cids(1), 4)}));
%% with average trace on the top
% figure;
% hold on
% nr = length(compids);
% Fz = 100;
% t = (0:size(dtest, 2)-1)/Fz;
% cliplines = cumsum(ClipLens-90)/Fz;
% cliplinecenter = mean([[0; cliplines(1:(end-1))] cliplines], 2);
% nClip = length(cliplines);
% 
% 
% for i = 1:nr
%     fcompids = find(cTab(:, 2)== compids(i));
%     plot(t, mean(targetD1(fcompids, :), 1)-(i-1)*2, 'Color', [2 105 164]/255);
%     plot(t, mean(targetD2(fcompids, :), 1)-(i-1)*2-1, 'Color', [0 180 171]/255);
%     text(-2, -i*2+2, sprintf('%d', i), 'Color', 'k');
% end
% % text(-2, 2, 'average', 'Color', 'k');
% for i = 1:nClip
%     text(cliplinecenter(i), 3.5, sprintf('%d', i), 'Color', 'k');
% end
% plot(cliplines(1:end-1)'.*ones(2, nClip-1), [-nr*2+1; 3].*ones(2, nClip-1), '--k');
% avgtrace = mean(targetD1, 1);
% avgtrace = avgtrace - min(avgtrace);
% avgtrace = avgtrace/max(avgtrace);
% plot(t, avgtrace+2, 'Color', [2 105 164]/255, 'LineWidth', 2);
% avgtrace = mean(targetD2, 1);
% avgtrace = avgtrace - min(avgtrace);
% avgtrace = avgtrace/max(avgtrace);
% plot(t, avgtrace+1, 'Color', [0 180 171]/255, 'LineWidth', 2);
% xlabel('Time (s)');
% box off
% axis off
% sgtitle(sprintf('%d. %d-%d (%d, %d) (%s)', ii, Day, Cel,  expids(1), expids(2), BCTypeLabels{DTab(cids(1), 4)}));