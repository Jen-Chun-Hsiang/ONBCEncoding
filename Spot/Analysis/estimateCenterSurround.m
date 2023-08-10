function [params, fitcurve, costv, errs] = estimateCenterSurround(inputs, options)

nROI = size(inputs, 1);

IsPCA = 0;
RmvPeakRange = [1:5 14:30];
IsDebug = 0;
DataOrgnization = [7 30];%[7 30][7 33]
nTime = DataOrgnization(2);
EstimationType = 'DoG'; % 'DoG', 'Amp'
Dmt = [10 25 50 100 200 300 400];
AmpDiameterSets = {[2 3], [6 7]};
ClusterDiameterSets = 1:7;% [4 5 6 7];
TemporalSets = {[7:10], [15:18]};
ContrastTemporalSets = {5:16, 23:29};
if nargin == 2
    if isfield(options, 'IsPCA'), IsPCA = options.IsPCA;end
    if isfield(options, 'RmvPeakRange'), RmvPeakRange = options.RmvPeakRange;end
    if isfield(options, 'IsDebug'), IsDebug = options.IsDebug;end
    if isfield(options, 'DataOrgnization'), DataOrgnization = options.DataOrgnization;end
    if isfield(options, 'Dmt'), Dmt = options.Dmt;end
    if isfield(options, 'EstimationType'), EstimationType = options.EstimationType;end
end


costv = nan(nROI, 1);
fitcurve = nan(nROI, length(Dmt));
errs = zeros(1, 4);
for j = 1:nROI
    % Get [Time, Dim]
    CellTraces = reshape(inputs(j, :), DataOrgnization(2), []);
    if size(CellTraces, 2) ~= DataOrgnization(1)
        error('Data structure is not matched!');
    end
    
    if IsPCA
        [Coef, Score, ~, ~, explained, mu] = pca(CellTraces);
        Sco = Score;
        Sco(:, 2:end) = 0;
        Coe = Coef;
        Coe(:, 2:end) = 0;
        RCellTraces = (Sco*Coe'+mu);
    else
        RCellTraces = CellTraces;
    end
    
    FocusRCellTraces = RCellTraces;
    FocusRCellTraces(RmvPeakRange, :) = 0;
    [~, Peak] = max(FocusRCellTraces(:));
    [PeakId(1), PeakId(2)] = ind2sub(size(RCellTraces), Peak);
    TPoints = PeakId(1)-1:PeakId(1)+2;
    
    switch lower(EstimationType)
        case 'dog'
            if j == 1
                params = nan(nROI, 5);
            end
            %%
            RDmt = Dmt/1000;
            try
                FitCurve = mean(RCellTraces(TPoints, :), 1);
                fitcurve(j, :) = FitCurve;
                Gauss1D = @(x, s) exp(-x.^2/(2*s^2))/(s*sqrt(2*pi));
                CostF = @(w) sqrt(mean((cumsum(w(2)*Gauss1D(RDmt, w(1)) - w(4)*Gauss1D(RDmt, w(3)) + w(5))-FitCurve).^2));% + w(4)/7;
                % Opt.MaxIter = 20;
                % Opt.Method = 'newton0lbfgs';
                % Opt.numDiff = 1;
                % [OptW,f,exitflag,output] = minFunc(CostF, [0.02, 0.2, 0.05, 0.1, min(FitCurve)]', Opt);
                %         [OptW,fval] = fminsearch(CostF,[0.02, 0.2, 0.05, 0.1, min(FitCurve)]');
                [params(j, :),~] = fmincon(CostF, [0.02, 0.2, 0.05, 0.1, min(FitCurve)], [], [], [], [],...
                    [0       0     0    0             -2],...
                    [5       2     5    2              2]);
                
                costv(j) = CostF(params);
                clc
            catch
                errs(j) = 1;
            end
            
            if IsDebug
                close all; figure;
                subplot(1, 3, 1);
                plot(CellTraces);
                subplot(1, 3, 2);
                plot(RCellTraces);
                hold on;
                plot(PeakId(1), RCellTraces(PeakId(1), PeakId(2)), 'rx');
                subplot(1, 3, 3);
                plot(RDmt,FitCurve,'b');
                hold on;plot(RDmt, params(j, 2)*Gauss1D(RDmt, params(j, 1)), 'k');
                hold on;plot(RDmt, -params(j, 4)*Gauss1D(RDmt, params(j, 3)), 'r');
                hold on; plot(RDmt, cumsum(params(j, 2)*Gauss1D(RDmt, params(j, 1)) -params(j, 4)*Gauss1D(RDmt, params(j, 3))+params(j, 5)), 'k.');
                
                legend({'Tuning', 'Center', 'Surround', 'Fitting Model'});
            end
        case 'amp'
            try
                mTraces = mean(RCellTraces(TPoints, :), 1);
                fitcurve(j, :) = mTraces;
                if j == 1, params = nan(nROI, 1);end
                mTraces = mTraces-min(mTraces);
                CenterAmp = mean(mTraces(AmpDiameterSets{1}));
                SurroundAmp = mean(mTraces(AmpDiameterSets{2}));
                params(j) = (CenterAmp-SurroundAmp)/(CenterAmp+SurroundAmp+eps);
%                 params(j) = CenterAmp/SurroundAmp;
            catch
                errs(j) = 1;
            end
        case 'linearity'
            if j == 1
                params = nan(nROI, 1);
            end
            mTraces = mean(RCellTraces(:, 2:4), 2);
            ONAmp = mean(mTraces(ContrastTemporalSets{1}));
            OFFAmp = mean(mTraces(ContrastTemporalSets{2}));
%             params(j) = (OFFAmp)/(abs(ONAmp)+abs(OFFAmp)+eps);
            params(j) = -OFFAmp*2/(abs(OFFAmp)+abs(ONAmp)+eps);
        case 'surround'
            if j == 1
                params = nan(nROI, 1);
            end
            mTraces = mean(RCellTraces(6:17, :), 1);
            CenterAmp = mean(mTraces(AmpDiameterSets{1}));
            SurroundAmp = mean(mTraces(AmpDiameterSets{2}));
%             params(j) = SurroundAmp;
%             params(j) = (SurroundAmp)/(abs(CenterAmp)+abs(SurroundAmp)+eps);
            params(j) = -(SurroundAmp-CenterAmp)/(abs(CenterAmp)+eps);
        case 'transience'
            if j == 1
                params = nan(nROI, 1);
            end
            mTraces = mean(RCellTraces(:, 2:4), 2);
            EarlyAmp = mean(mTraces(5:9));
            LateAmp = mean(mTraces(14:18));
%             params(j) = LateAmp/(abs(EarlyAmp)+abs(LateAmp)+eps);
            params(j) = (EarlyAmp - LateAmp)/(EarlyAmp+eps);
    end
    fprintf('processing... ROI (%d/%d) \n', j, nROI);
end


end