celltype = [50 51 57 58 6 7 89];
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
%% ON
cts = 1;
celltype = [50 51 57 58 6 7 89];
peakv = [];
T = (-BasDurFrm:StmDurFrm)/Fz;
Twin = T>0.1 & T<1.5;
ii = 0;
jj = [];
for t = 1:7
    ctyp = celltype(t);
    rois = find(Typ == ctyp);
    Vars = DataROIsDfnrm;
    for k = 1:5 % for each spot size
        crois = rois(squeeze(RepReli(k, cts, rois))>0.5);
        jj = [jj; length(rois), length(crois)];
        a = squeeze(Vars(Twin, :, k, cts, crois));
        [~, sids] = sort(squeeze(RepReli(k, cts, crois)), 'descend');
        for i = 1:length(crois) % for each experiment
            b = squeeze(a(:, :, sids(i)));
            ii = ii + size(b, 2) ;
            for j = 1:size(b, 2) % for each repeat
                c = b(:, j)';
                [pks, locs] = findpeaks(c);
                if ~isempty(pks)
                    peakv = [peakv; max(pks) t];
%                     peakv = [peakv; pks(:) t*ones(length(pks), 1)];
                end
            end
        end
    end
end

%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
x = -0.2:0.05:1;
figure;
for i = 1:7
    subplot(2, 4, i); hold on
    data = peakv(peakv(:, 2)==i, 1);
    data(isnan(data))=[];
%     estimated_p = fitBinomial(data);
    [f, x, flo] = ecdf(data,'Function','cdf','Alpha',0.05,'Bounds','on');
    
    shadePlot(x, f, abs(flo-f), 'k', 0.5);
    plot(zeros(1, 2), [0 1], '--k');
%     ecdf(peakv(peakv(:, 2)==i, 1),'Function','cdf','Alpha',0.01,'Bounds','on');
    box off
    title(BCTypelabels{i});
    xlim([-0.5 1.01]);
    xticks(-0.5:0.5:1);
    xticklabels({'-0.5', '0', '0.5', '1'});
    xlabel('Peak amplitude');
    ylim([0 1.01]);
    yticks(0:0.25:1);
    yticklabels({'0', '', '0.5', '', '1'});
    ylabel('Cumulative probability');
end
subplot(2, 4, 8); hold on
data = peakv(:, 1);
data(isnan(data))=[];
estimated_p = fitBinomial(data(data>0));
[f, x, flo] = ecdf(data,'Function','cdf','Alpha',0.05,'Bounds','on');
shadePlot(x, f, abs(flo-f), 'k', 0.5);
 plot(zeros(1, 2), [0 1], '--k');
y = cdf('Binomial', x, 1, estimated_p);
plot(x, y, 'b');
y = cdf('Uniform', x, 0, 1);
plot(x, y, 'r');
box off
xlim([-0.5 1.01]);
xticks(-0.5:0.5:1);
xticklabels({'-0.5', '0', '0.5', '1'});
xlabel('Peak amplitude');
ylim([0 1.01]);
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylabel('Cumulative probability');
title('overall');
legend({'Observed data', '95% confidence', '', 'Binomial', 'Uniform'});
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sSupFig3x_NoSpikeinBCs', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);

%% OFF
cts = 2;
Vars = DataROIsDfnrm;
celltype = [50 51 57 58 6 7 89];
peakv = [];
T = (-BasDurFrm:StmDurFrm)/Fz;
Twin = T>0.1 & T<1.5;

% figure;
ii = 0;
for t = 1:7
    ctyp = celltype(t);
    rois = find(Typ == ctyp);
    Vars = DataROIsDfnrm;
    for k = 1:5 % for each spot size
        crois = rois(RepReli(k, 1, rois)>0.5);
        a = squeeze(Vars(Twin, :, k, cts, crois));
        [~, sids] = sort(squeeze(RepReli(k, 1, crois)), 'descend');
        for i = 1:length(crois) % for each experiment
            b = squeeze(a(:, :, sids(i)));
            ii = ii + size(b, 2) ;
            for j = 1:size(b, 2) % for each repeat
                c = b(:, j)';
                [pks, locs] = findpeaks(c);
                if ~isempty(pks)
                    [pks, maxid] =max(pks);
                    peakv = [peakv; pks t locs(maxid) k, crois(sids(i))];
                end
            end
        end
    end
end


%%
figure;
for i = 1:7
    subplot(2, 4, i); hold on
    data = peakv(peakv(:, 2)==i, 1);
    data(isnan(data))=[];
    [f, x, flo] = ecdf(data,'Function','cdf','Alpha',0.05,'Bounds','on');
    shadePlot(x, f, abs(flo-f), 'k', 0.5);
    box off
    title(BCTypelabels{i});
    
    plot(zeros(1, 2), [0 1], '--k');
    xlim([-1 1.01]);
    xticks(-1:0.5:1);
    xticklabels({'-1', '', '0', '', '1'});
    ylim([0 1.01]);
    yticks(0:0.25:1);
    yticklabels({'0', '', '0.5', '', '1'});
    ylabel('Cumulative probability');
    xlabel('Peak amplitude');
end
subplot(2, 4, 8); hold on
data = peakv(:, 1);
data(isnan(data))=[];
[f, x, flo] = ecdf(data,'Function','cdf','Alpha',0.05,'Bounds','on');
shadePlot(x, f, abs(flo-f), 'k', 0.5);
xlim([-0.8 1.01]);
plot(zeros(1, 2), [0 1], '--k');
xlim([-1 1.01]);
xticks(-1:0.5:1);
xticklabels({'-1', '', '0', '', '1'});
xlabel('Peak amplitude');
box off
ylim([0 1.01]);
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylabel('Cumulative probability');
title('overall');
sgtitle('OFF');

%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sSupFig3x_NoSpikeinBCs_OFF', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);

