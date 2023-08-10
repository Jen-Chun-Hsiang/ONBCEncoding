
load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\ProcessedData\RFSpatialNL_04242023\Ai148_Grm6Cre_NoiseSpatialNonlinear_radius40_110818_1_5.mat')
%% decide the targeting position
figure;
plot(mean(Rd_posi(:, 3, :), 3, 'omitnan'));
%% decide the targeting ROIs
figure;
plot(squeeze(Rd_posi(221, 3, :)));
%%
nROI = size(RF_locs, 1);
nposi = size(Rd_posi, 1);
IsDisplay = 1;
%     close all
RF_locs = nan(nROI, 6, 3);
exts = -44:4:44;
exts_fine = -44:44;
j = 3;
i = 33;
v = reshape(squeeze(Rd_posi(:, j, i)), sqrt(nposi), []);
assert(~any(isnan(v(:))));
v = [v(:, 1) v v(:, end)];
v = [v(1, :); v; v(end, :)];
[x, y] = meshgrid(exts, exts);
[xx, yy] =meshgrid(exts_fine, exts_fine);
Vq = interp2(x, y,...
    v, xx, yy, 'spline');
[RF_locs(i, 4, j), maxid] = max(Vq(:));
RF_locs(i, 5, j) = xx(maxid);
RF_locs(i, 6, j) = yy(maxid);
[RF_locs(i, 1, j), maxid] = max(v(:));
RF_locs(i, 2, j) = x(maxid);
RF_locs(i, 3, j) = y(maxid);

if IsDisplay
                close all
    figure;
    subplot(1, 2, 1);
    imagesc(v); colorbar
    hold on
    plot(find(exts == RF_locs(i, 2, j)), find(exts == RF_locs(i, 3, j)), 'ok');
    title(sprintf('ROI: %d (min: %0.3G, max:%0.3G)', i, min(v(:)), max(v(:))));
    box off
    axis off
    
    subplot(1, 2, 2);
    imagesc(Vq); colorbar;
    hold on
    plot(find(exts_fine == RF_locs(i, 5, j)), find(exts_fine == RF_locs(i, 6, j)), 'ok');
    title(sprintf('ROI: %d', i));
    box off
    axis off
end

%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure5_NatMovProjection_RFNoiseMap_ROI%d', SaveFolder, i);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);

%% after determined by matching table
% and have run to %DrawROIinImage_Process
figure; 
imshow(basec); hold on
ROIids = [12 10 36]; % determined by matchingROITable
for i = 1:length(ROIids)
    [B,L] = bwboundaries(cmap == ROIids(i),'noholes');
    boundary = B{1};
    plot(boundary(:,2), boundary(:,1), 'r');
    Cent = getROIloc(cmap == ROIids(i));
    text(Cent(1), Cent(2), num2str(i), 'Color', 'y');
end
% [B,L] = bwboundaries(cmap == 10,'noholes');
% boundary = B{1};
% plot(boundary(:,2), boundary(:,1), 'r')
% Cent = getROIloc(cmap == 12);
% text(Cent(1), Cent(2), '2', 'Color', 'y');
% 
% [B,L] = bwboundaries(cmap == 36,'noholes');
% boundary = B{1};
% plot(boundary(:,2), boundary(:,1), 'r')
% Cent = getROIloc(cmap == 12);
% text(Cent(1), Cent(2), '3', 'Color', 'y');
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure5_NatMovProjection_RFNoiseMap_ROISelect', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
