% goal of this script is to run a boostrapping to estimate how much the
% variation of each block related to the total 
% following by run @Response_BCPlexusComparison & @IdentifiedBCComparison
%%
[~, g] = max(SliceResp_type == BCTypes, [], 2);
%% Test Temporal distribution
nIter =  2000;
y = SimpleTemporalRio;
upbids = 1:4;
lwbids = 5:7;
up_pars = nchoosek(upbids, 2);
lw_pars = nchoosek(lwbids, 2);
nupbids = length(upbids);
nlwbids = length(lwbids);
nup_pars = size(up_pars, 1);
nlw_pars = size(lw_pars, 1);
upb = nan(nIter, nup_pars);
lwb = nan(nIter, nlw_pars);
for i = 1:nIter
    % upper block
    ty = nan(1, nupbids);
    for j = 1:nupbids
        ty(j) = randsample(y(g==upbids(j)), 1);
    end
    for j = 1:nup_pars
        upb(i, j) = sqrt(diff(ty(up_pars(j, :))).^2);
    end
    
     % lower block
    ty = nan(1, nlwbids);
    for j = 1:nlwbids
        ty(j) = randsample(y(g==lwbids(j)), 1);
    end
    for j = 1:nlw_pars
        lwb(i, j) = sqrt(diff(ty(lw_pars(j, :)-4)).^2);
    end
    
    fprintf('progress...%d/%d', i, nIter);
end
%%

[~, p2, k2] = kstest2(upb(:), lwb(:), 'Tail', 'smaller')
%%
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(7);
upColor = mean(Colors(1:4, :), 1);
lwColor = mean(Colors(5:7, :), 1);
figure; hold on
[f,x,flo,fup] = ecdf(upb(:),'Alpha',0.01,'Bounds','on');
shadePlot(x, f, abs(flo-f), upColor, 0.3, 1);
[f,x,flo,fup] = ecdf(lwb(:),'Alpha',0.01,'Bounds','on');
shadePlot(x, f, abs(flo-f), lwColor, 0.3, 1);
xlabel('Transience difference (averaged)');
legend({'', 'Upper', '',  'Lower'})

xlim([min([upb(:); lwb(:)]) max([upb(:); lwb(:)])]);
xticks(0:1:2);
xticklabels({'0', '1', '2'});
ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0', '0.5', '1'});
title('Transience variation between blocks');


%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure3_TransienceBlockVariation_ECDF', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%% Test Surround distribution
nIter =  2000;
y = SimpleCSRio;
upbids = 1:4;
lwbids = 5:7;
up_pars = nchoosek(upbids, 2);
lw_pars = nchoosek(lwbids, 2);
nupbids = length(upbids);
nlwbids = length(lwbids);
nup_pars = size(up_pars, 1);
nlw_pars = size(lw_pars, 1);
upb = nan(nIter, nup_pars);
lwb = nan(nIter, nlw_pars);
for i = 1:nIter
    % upper block
    ty = nan(1, nupbids);
    for j = 1:nupbids
        ty(j) = randsample(y(g==upbids(j)), 1);
    end
    for j = 1:nup_pars
        upb(i, j) = sqrt(diff(ty(up_pars(j, :))).^2);
    end
    
     % lower block
    ty = nan(1, nlwbids);
    for j = 1:nlwbids
        ty(j) = randsample(y(g==lwbids(j)), 1);
    end
    for j = 1:nlw_pars
        lwb(i, j) = sqrt(diff(ty(lw_pars(j, :)-4)).^2);
    end
    
    fprintf('progress...%d/%d', i, nIter);
end
%%
[~, p2, k2] = kstest2(upb(:), lwb(:), 'Tail', 'larger')
%%
figure; hold on
[f,x,flo,fup] = ecdf(upb(:),'Alpha',0.01,'Bounds','on');
shadePlot(x, f, abs(flo-f), upColor, 0.3, 1);
[f,x,flo,fup] = ecdf(lwb(:),'Alpha',0.01,'Bounds','on');
shadePlot(x, f, abs(flo-f), lwColor, 0.3, 1);
legend({'', 'Upper', '',  'Lower'})
xlabel('Surround strength difference (averaged)');
xlim([min([upb(:); lwb(:);]) max([upb(:); lwb(:); ])]);
xticks(0:0.8:1.6);
xticklabels({'0', '0.8', '1.6'});
ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0', '0.5', '1'});
title('Surround strength variation between blocks');
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure3_SurroundStrengthBlockVariation_ECDF', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);