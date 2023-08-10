
%%
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
%%
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
cqualthr = qualthr;
% cqualthr = 0;
close all
figure; hold on
ctyp = typ(any(Qual > cqualthr & Rept > 1, 2));
targetD = BSctr;
targetD(Qual < cqualthr) =  nan;
targetD = mean(targetD, 2, 'omitnan');
targetD = targetD(any(Qual > cqualthr, 2) & Rept > 1);
cqual = Qual;
cqual(Qual < cqualthr) =  nan;
cqual = mean(cqual, 2, 'omitnan');
cqual = cqual(any(Qual > cqualthr, 2) & Rept > 1);
mv = nan(7, 1);
se = nan(7, 1);
ql = nan(7, 1);
for i = 1:7
   cids = ctyp ==  i;
   nc = sum(cids);
   mv(i) = mean(targetD(cids));
   se(i) = std(targetD(cids))/sqrt(nc);
%    mv(i) = mean(targetD(cids)./cqual(cids));
%    se(i) = std(targetD(cids)./cqual(cids))/sqrt(nc);
   ql(i) = mean(cqual(cids));
end
b = bar(1:7, mv);
b.EdgeColor = 'w';
b.FaceColor = 'flat';
Colors = dazure(7);
e =errorbar(1:7, mv, se, 'LineStyle', 'none');
e.CapSize = 0;
e.Color = 0.5*ones(1, 3);
for i = 1:7
    b.CData(i,:) = Colors(i, :);
    plot(i+[-1 1]*0.4 , ql(i)*ones(1, 2), '--', 'Color', 0*ones(1, 3), 'LineWidth', 1);
end
xticks(1:7);
xticklabels(BCTypelabels);
yticks(0:0.3:0.6);
yticklabels({'0', '0.3', '0.6'});
ylabel('Explained variance');
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure5_LuminanceExplainedVarianceNatMovCellType', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);

%%
% figure; 
stits = {'Mean Intensity', 'Spatial nonlinearlity', 'Motion sensitivity', 'Surround strength','Transience', 'Sustain'};
% for i = 1:5
%     subplot(1, 5, i);
%    imagesc( squeeze(lagD(:, :, i))); colorbar
%    title(stits{i});
% end

%% draw feature responses (comparing the level with above lag)
% hypothesis, the stronger the feature the more precisely reflecting lags
figure;  hold on
t = (0:length(clipids)-1)/100;
for i = 1:6
    switch i
        case 1
            a = Features(findcol('RF_L'),:);
        case 2
            a = Features(findcol('RF_NL'),:);
        case 3
            a = Features(findcol('MS_VAd'),:);
        case 4
            a = Features(findcol('RF_CS'),:);
        case 5
            a = Features(findcol('TS_HP'),:);
        case 6
            a = Features(findcol('TS_LP'),:);
    end
    a = 0.9*(a-min(a))/range(a);
    plot(t, a-i+1);
%     yticks(0.5:1:6.5);
%     yticklabels(stits);
    text(1, -i+1.9, stits{i}, 'Color', Colors(i, :), 'FontSize', 12);
end
for i = 1:nclip
    ids = find(clipids' == i);
    plot(ids(end)*ones(1, 2)/100, [-5 1], '--k');
end
box off
axis off
xlim([0 t(end)]);
xlabel('Time (s)');
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure5_CovolutedVisualFeatureTrace', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
