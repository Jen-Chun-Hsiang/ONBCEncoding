close all;clear;clc;
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions\OpenSource');
figure('Name','sankey demo4_1','Units','normalized','Position',[.05,.2,.5,.56])
% (Unequal inflow and outflow data)
links={'Injected at P4-P6(378 mice)','dead before imaging',0.35;
       'Injected at P4-P6(378 mice)','incorrect genetype',1.82;
       'Injected at P4-P6(378 mice)','Imaged (161 mice)',1.61;
       'Imaged (161 mice)','damaged retina',0.06;
       'Imaged (161 mice)','inappropriate expression',0.79;
       'Imaged (161 mice)','Recorded (76 mice - 223 cells)',0.76;
       'Recorded (76 mice - 223 cells)','unidentifiable',0.07*0.76/2.23; %0.03
       'Recorded (76 mice - 223 cells)','insufficient response quality (rod bipolar)',0.12*0.76/2.23; %0.12
       'Recorded (76 mice - 223 cells)','insufficient response quality (others)',1.47*0.76/2.23; % 1.47
       'Recorded (76 mice - 223 cells)','cone bipolar (41 mice - 57 cells)',0.57*0.76/2.23;
       };
% links{16,3}=.5;
% (Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));
SK.ColorList = [zeros(2, 3); 0.9*ones(2, 3); zeros(1, 3); 0.9*ones(2, 3); zeros(1, 3); 0.9*ones(3, 3)];
SK.NodeList={'Injected at P4-P6(378 mice)','Imaged (161 mice)','incorrect genetype','dead before imaging',...
    'Recorded (76 mice - 223 cells)','inappropriate expression','damaged retina',...
    'cone bipolar (41 mice - 57 cells)','insufficient response quality (others)',...
    'insufficient response quality (rod bipolar)','unidentifiable'};
% SK.NodeList={'Injected (400 mice)','> 4.5 months','dead before imaging','wrong genetype','Imaged (161 mice)',...
%     'damaged retina','no isolation','only-ventrally expressed','off-targert expression','weak expression',...
%     'no expression','Recorded (223 cells)',...
%     'unclassifiable','bad responses (RBC)','bad response (Others)','Included (61 cells',...
%     'low quality','no cone contact data','Final images (47 cells - 34 mice)'};
SK.RenderingMethod='interp';
SK.draw();
% 'Imaged (161 mice)','no isolation',0.09;
%        'Imaged (161 mice)','only-ventrally expressed',0.11;
%        'Imaged (161 mice)','off-targert expression',0.15;
%        'Imaged (161 mice)','weak expression',0.17;
%'Included (41 mice - 61 cells)','low quality',0.04*0.74/2.23;
%'Included (41 mice - 61 cells)','no cone contact data',0.1*0.74/2.23;
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure1_SankeyDiagramBCExperiment', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);