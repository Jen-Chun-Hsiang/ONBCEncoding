close all;  clc;
Fz = 52;
SelectedClip = 1:11;
%% Simple testing Copy "a" from excel sheets
if exist('a')
    Frms = a(:, 1);
    MarkX = a(:, 3);
    PhysX = [a(1, 5) a(end, 5)];
    PhysX = MarkX*(PhysX(end)-PhysX(1))/(MarkX(end)-MarkX(1));
    iFrms = Frms(1):Frms(end);
    iPhysX = interp1(Frms, PhysX, iFrms, 'spine');
    
    Speed = diff(iPhysX)/(1/Fz);
    Speed = [Speed(1) Speed];
    figure;
    yyaxis left
    plot(Frms, PhysX, 'k'); hold on
    plot(iFrms, iPhysX, 'b'); hold on
    ylabel('Position (arbi. unit)');
    xlabel('Frame number');
    yyaxis right
    plot(iFrms, Speed);
    ylabel('Speed (feet/second)');
    keyboard;
end

%%
FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\NaturalVideoSpeedMeasurement.xlsx';
[DocNum, DocStr] = xlsread(FilNam,'Sheet1', 'A1:G311');
ColName2Ind = @(String) find(strcmpi(DocStr(1, :), String));
ClmClip =  ColName2Ind('ClipNum');
ClmAct =  ColName2Ind('ActionNum');
ClmMark = ColName2Ind('MarkerId');
ClmFrm = ColName2Ind('FrmNum');
ClmMarkX = ColName2Ind('MarkerPositionX');
ClmPhysi = ColName2Ind('PhysiPosition');

NumClip = length(SelectedClip);
ColName2Ind('PhysiPosition')

SpeedTable =[];
for i = 1:NumClip
    CurClip = DocNum(DocNum(:, ClmClip) == SelectedClip(i), :);
    Actions = unique(CurClip(:, ClmAct));
    NumAct = length(Actions);
    for j = 1:NumAct
        CurAct = CurClip(CurClip(:, ClmAct) == Actions(j), :);
        Events = unique(CurAct(:, ClmMark));
        NumEvent = length(Events);
        for k = 1:NumEvent
            CurEvent = CurAct(CurAct(:, ClmMark) == Events(k), :);
            
            Frms = CurEvent(:, ClmFrm);
            MarkX = CurEvent(:, ClmMarkX);
            PhysX = [CurEvent(1, ClmPhysi) CurEvent(end, ClmPhysi)];
            PhysX = MarkX*(PhysX(end)-PhysX(1))/(MarkX(end)-MarkX(1));
            iFrms = Frms(1):Frms(end);
            iPhysX = interp1(Frms, PhysX, iFrms, 'pchip');
            
            Speed = diff(iPhysX)/(1/Fz);
            Speed = [Speed(1) Speed];
            DataL = length(Speed);
            
            SpeedTable = [SpeedTable; Speed(:), iPhysX(:), SelectedClip(i)*ones(DataL, 1),  Actions(j)*ones(DataL, 1),...
                Events(k)*ones(DataL, 1) [1; zeros(DataL-1, 1)]];
        end
    end
    fprintf('Progress... %d / %d(Clip) \n', i, NumClip);
end

%%
Feet2cm = 30.48;
x = SpeedTable(:, 1)*Feet2cm;
hlims = [min(x) max(x)];
[f, xi] = ksdensity(x);
% figure; 
% plot(xi, f); 
close all
figure; 
h = histogram(x, 30);
x = reshape(repmat(h.BinEdges, 2, 1), [], 1);
y = reshape(repmat( h.Values, 2, 1), [], 1);
plot(x(2:end-1), y, 'k');
box off
% h.Normalization = 'Probability';
% h.EdgeColor = 'b';
% h.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:60:180);
xticklabels({'0', '60', '120', '180'});
yticks(0:400:800);
yticklabels({'0', '400', '800'});
ylabel('Number of frames');
xlabel('Speed (cm /s)');
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sSupFig6_NaturalMovie_movingspeeddistribution', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
