% run after @DirectDistanceDataAugmentationTypeComparison
Ids = any(DepthTypeProb>0.25, 2);
A = distprob(nBC+1:end, 1:nBC);
A = A(Ids, :);
A = A(:, :);
B = DepthTypeProb(Ids, dtypes);
B = B(:, :);
C = B.*A;
C = C./max(C, [], 2);
[d, D] = max(B.*A, [], 2);
D = dtypes(D);
D(d > maxx) = 8;
[~, sids] = sort(D, 'ascend');

%%
figure; 
subplot(3, 1, 1);

imagesc(A(sids, :)'); colormap(gray); colorbar; 
title(sprintf('Min-max(%.03G ,  %.03G) size(%d, %d)', min(A(:)), max(A(:)), size(A, 1), size(A, 2)));
box off
axis off


subplot(3, 1, 2);
imagesc(B(sids, :)'); colormap(gray); colorbar; 
title(sprintf('Min-max(%.03G ,  %.03G)  size(%d, %d)', min(B(:)), max(B(:)), size(A, 1), size(A, 2)));
box off
axis off


subplot(3, 1, 3);
ids = [53 98 154];
imagesc(C(sids, :)'); colormap(gray); colorbar; 
% for i = 1:length(ids)
%     rectangle('Position', [ids(i)-0.5, 0, 1, size(C, 2)], 'EdgeColor', 'r');
% end
title(sprintf('Min-max(%.03G ,  %.03G)  size(%d, %d)', min(C(:)), max(C(:)), size(C, 1), size(C, 2)));
box off
axis off
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure4_SpotProjection_ProbabilityMap', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);

%%
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = [dazure(nType); 0*ones(1, 3)];
Dr = D(sids);
ids = find(Dr ==8);
figure; 
imagesc(Dr'); colormap(Colors);colorbar
% for i = 1:length(ids)
%     rectangle('Position', [ids(i)-0.5, 0, 1, size(C, 2)], 'EdgeColor', 'k');
% end
box off
axis off

%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure4_SpotProjection_Classification', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);

%% use @getIPLProb
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(nType);
figure; hold on
for i = 1:7
    plot(IPLx, IPLd(i, :), 'Color', Colors(i, :));
end
xlim([0.47 0.95]); 
xticks(0.5:0.2:0.9);
xticklabels({'0.5', '0.7', '0.9'});
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
xlabel('IPL depth');
ylabel('Prob (norm.)');
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure4_SpotProjection_normalizedDepthProbability', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);