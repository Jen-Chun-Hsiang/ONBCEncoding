clear; close all; clc;
%%
FileIdentifier = 'BC';
PathName = uigetdir('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BC diversity morphology data\', 'Select the IPL profile folder');
FileNames = dir(fullfile(PathName, '*.csv'));    
FileIds = find(contains({FileNames.name},FileIdentifier));
% TypeDistinct = {'BC5i_o', 'BC6', 'BC7', 'XBC', 'BC8_9', 'BC5t'};
% BCTypeNum = [5, 6, 7, 58, 89, 57];
% TypeLabels = {'BC5io', 'BC6', 'BC7', 'XBC', 'BC8-9', 'BC5t'};
%%
TypeDistinct = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC8_9'};
BCTypeNum = [50, 51, 57 58 6, 7, 89];
TypeLabels = TypeDistinct;
%%
nFile = length(FileIds);
%%
DTable = nan(nFile, 4);
clear Profile BCType BCDay
for f = 1:nFile
    close all; clc;
    FileName = FileNames(FileIds(f)).name;
    Words = strsplit(FileName, '-');
    Day = Words{1};
    Cel = Words{2};
    Words = strsplit(Words{3}, '.');
    Type = Words{1};
    Words = strsplit(FileName, '.');
    SheetName = Words{1};
    [DocNum, DocStr] = xlsread([PathName '\' FileName],SheetName, 'A:D');

    cIds = all(~isnan(DocNum), 2);
    Profile{f} = DocNum(cIds, :);
    BCType{f} = Type;
    BCDay{f} = Day;
    ctype = nan;
    if isempty(BCTypeNum(strcmpi(TypeDistinct, Type)))
        ctype = nan;
    else
        ctype = BCTypeNum(strcmpi(TypeDistinct, Type));
    end
    DTable(f, :) = [str2double(Day), str2double(Cel),...
        ctype, size(DocNum, 2)];
end
% a = cellfun(@(x) BCTypeNum(strcmpi(TypeDistinct, x)), BCType);
% DTable(:, 3) = a(:);
%%
PathName = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BC diversity morphology data\';
FileName = 'CollectionBCIPLData.mat';
%%
save([PathName FileName], 'Profile', 'BCType', 'BCDay', 'DTable');

%%
load([PathName FileName], 'Profile', 'BCType', 'BCDay', 'DTable');

%%
xi = linspace(0, 100, 101)/100;
nCel = length(Profile);
sprofile = nan(nCel, length(xi));
for i = 1:nCel
    x = Profile{i}(:, 2)/100;
    v = Profile{i}(:, 3);
    v = v/sum(v(:));
    sprofile(i, :) = interp1(x, v, xi, 'pchip');
end
%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(length(BCTypeNum));
transparency = 0.5;
mV = nan(length(BCTypeNum), size(sprofile, 2));
sem = nan(length(BCTypeNum), size(sprofile, 2));
% DispBC = [5, 57, 58, 6, 7, 89];
DispBC = [50, 51, 57, 58, 6, 7, 89];
for i = 1:length(BCTypeNum)
   subplot(1, length(BCTypeNum), i);
   cIds =  find(DTable(:, 3) == DispBC(i));
   cN = length(cIds);
   mV(i, :) = mean(sprofile(cIds, :), 1);
   sem(i, :) = std(sprofile(cIds, :), [], 1)/sqrt(cN);
end
close all
x = round((1-xi)*100);
figure; 
for i = 1:length(BCTypeNum)
   subplot(1, length(BCTypeNum), i);
%    plot(repmat(round((1-xi)*100), length(BCTypeNum), 1)', mV', 'Color', 0.5*ones(1, 3));
   shadePlot(x, mV(i, :), sem(i, :),Colors(i, :), transparency, 2)
   xlim([40 100]);
   box off
%    if i == 1
%        yticks(0:0.02:0.04);
%        yticklabels({'0', '0.02', '0.04'});
%    end
   yticks(0:0.02:0.04);
   yticklabels({'0', '0.02', '0.04'});
   xticks(40:30:100);
   xticklabels({'40', '70', '100'});
   xlabel('IPL depth(%)');
   ylabel('GFP intensity (norm.)');
   title(TypeLabels(BCTypeNum == DispBC(i)));
   ylim([0 0.05]);
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure1_BCtypeBCIPLDepth', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
IPLdepth = 1-xi;
% IPLprofile = (mV-min(mV, [], 2))./range(mV, 2);
IPLprofile = mV;

%%
PathName = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BC diversity morphology data\';
FileName = 'BCTypeIPLProfile.mat';
save([PathName FileName], 'IPLdepth', 'IPLprofile', 'DispBC', 'TypeLabels');
