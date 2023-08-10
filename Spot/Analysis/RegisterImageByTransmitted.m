clear; close all; clc;

%% Analysis Parameters
Day = '111522';
RegisterMethod = 2;
CleanMethod = 1;
AnchorFrameStep = 50; 
IsOrderReversed = 0;
IsPlot = 0;
%% Load Data
PathName = uigetdir('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicPreprocessing\Files\Preprocessed', 'Select the MatFile folder');
FileNames = dir(fullfile(PathName, '*.mat'));
FileIds = find(contains({FileNames.name},Day) & contains({FileNames.name},'_3.mat'));
nFile = length(FileIds);

SaveFolder = strsplit(PathName, 'Preprocessed');
SaveFolder = SaveFolder{1};
PathName = [PathName '\'];

if IsOrderReversed
    ReadSequence = fliplr(1:nFile);
else
    ReadSequence = 1:nFile;
end
for f = ReadSequence
    FileName = FileNames(FileIds(f)).name;
    Words = strsplit(FileName, ['_' Day '_']);
    Topic = Words{1};
    Words = strsplit(Words{2}, '_');
    Cel = floor(str2double(Words{1})/1000);
    Exp = mod(str2double(Words{1}),1000);
    
    RegistName = sprintf('%s_Translation_%s_%d%03d.mat',Topic, Day, Cel, Exp );
    if exist([SaveFolder '\Translation\' RegistName], 'file')
        disp('File has been registrated');
    else
        % Load calcium signals
        FilNam = sprintf('%s_%s_%d%03d_1',Topic, Day, Cel, Exp );
        load([PathName FilNam '.mat'], 'I');
        
        % Load Structual / Anatomy Image
        ImgAtm = load([PathName FileName]);
        ImgAtm = ImgAtm.I;
        
        % Get basic parameters
        NumFrm = size(I, 3);
        
        % Denoise / Sharpen images
        switch CleanMethod
            case 1
                ImgAtm = medfilt3(ImgAtm);
            case 2
                for i=1:NumFrm
                    ImgAtm(:, :, i) = imgaussfilt(ImgAtm(:, :, i), 1);
                    ImgAtm(:, :, i) = imsharpen(ImgAtm(:, :, i), 'Radius',2,'Amount',1);
                end
        end
        
        % Define the fixed frame
        MidFrm = round(NumFrm/2);
        
        % Start Registrator
        switch RegisterMethod
            case 1 % Cross correlation
                usfac = 100;
                addpath('./Functions');
                
                %
                Pixel_X = size(I, 1);
                Pixel_Y = size(I, 2);
                
                % Fourier Transform
                I_fft = nan(Pixel_X, Pixel_Y, NumFrm);
                ImgAtm_fft = nan(Pixel_X, Pixel_Y, NumFrm);
                
                for i = 1:NumFrm
                    I_fft(:, :, i) = fft2(I(:,:, i));
                    ImgAtm_fft(:, :, i) = fft2(ImgAtm(:,:, i));
                end
                
                FFTY = nan(NumFrm, Pixel_Y);
                FFTX = nan(NumFrm, Pixel_X);
                for i = 1:NumFrm
                    FFTX(i, :) = mean(abs(ImgAtm_fft(:, :, i)), 2);
                    FFTY(i, :) = mean(abs(ImgAtm_fft(:, :, i)), 1);
                end
                
                % Exclusion by X, Y
                [~, Score] = pca((FFTX-min(FFTX, [], 2))./range(FFTX, 2));
                figure; scatter(Score(:, 1), Score(:, 2), [], Colors);
                [x,y] = ginput(2);
                hold on; plot([min(x) min(x) max(x) max(x) min(x)], [min(y) max(y) max(y) min(y) min(y)]);
                ExcludedFrmsX = Score(:, 1) > max(x) | Score(:, 1) < min(x) |...
                    Score(:, 2) > max(y) | Score(:, 2) < min(y);
                
                [~, Score] = pca((FFTY-min(FFTY, [], 2))./range(FFTY, 2));
                figure; scatter(Score(:, 1), Score(:, 2), [], Colors);
                [x,y] = ginput(2);
                hold on; plot([min(x) min(x) max(x) max(x) min(x)], [min(y) max(y) max(y) min(y) min(y)]);
                ExcludedFrmsY = Score(:, 1) > max(x) | Score(:, 1) < min(x) |...
                    Score(:, 2) > max(y) | Score(:, 2) < min(y);
                ExcludedFrms = ExcludedFrmsX | ExcludedFrmsY;
                figure; scatter(Score(~ExcludedFrms, 1), Score(~ExcludedFrms, 2), [], 'b');hold on
                scatter(Score(ExcludedFrms, 1), Score(ExcludedFrms, 2), [], 'k');hold on
                legend({'Included', 'Excluded'});
                
                if ExcludedFrms(round(NumFrm/2))
                    error('Middle frame is excluded!');
                end
                
                RegistChange = zeros(NumFrm, 4);
                for i = 1:NumFrm
                    if ~ExcludedFrms(i)
                        [output, Greg, Params] = DFT_Registration(ImgAtm_fft(:,:, MidFrm),ImgAtm_fft(:, :, i),usfac);
                        RegistChange(i, :) = output(:)';
                        I(:,:,i) = abs(ifft2(I_fft(:, :, i).*Params.EntrywiseProduct*Params.MatrixProduct));
                        ImgAtm(:,:,i) = abs(ifft2(ImgAtm_fft(:, :, i).*Params.EntrywiseProduct*Params.MatrixProduct));
                    end
                end
                figure; scatter(RegistChange(:, 3), RegistChange(:, 4), 2*ones(NumFrm,1), Colors);
                xlim([-2, 2]);
                ylim([-2, 2]);
                save([SaveFolder RegistName], 'ImgAtm', 'I', 'ExcludedFrms','RegistChange');
                clear ImgAtm I ExcludedFrms RegistChange
            case 2 % Intensity-Based
                % Ignore the artifact
                ImgAtm = ImgAtm(1:end-2, 6:end-5, :);
                
                % Convert Image to grayscale (0~1)
                RangeImgAtm = range(ImgAtm(:));
                MinImgAtm = min(ImgAtm(:));
                ImgAtm = (ImgAtm-MinImgAtm)/RangeImgAtm;
                Fxd = ImgAtm(:,:,MidFrm);
                
                % Optimizer
                [Opt, Met]             = imregconfig('Multimodal');
                Opt.Epsilon            = 1.5e-6;
                Opt.MaximumIterations  = 300;
                Opt.InitialRadius      = 6.25e-4;
                
                % Estimate the center point, and adjust the translation
                % coordinates
                ProgressTransform = nan(3,3,NumFrm);
                for i = 1:NumFrm
                    Mov = ImgAtm(:,:,i);
                    RefRgt = imregtform(Mov, Fxd, 'translation', Opt, Met);
                    ProgressTransform(:, :, i) = RefRgt.T;
                    if mod(i, 20)==0
                        clc
                        fprintf('Files: %d/%d \n', f, nFile);
                        fprintf('\t Step: %d/%d (estimation)', 1, 2);
                        fprintf('\t\t Frames: %d/%d \n', i, NumFrm);
                    end
                end
                xT = medfilt1(squeeze(ProgressTransform(3, 1, :)), 3);
                yT = medfilt1(squeeze(ProgressTransform(3, 2, :)), 3);
                CenterT = [mean(xT) mean(yT)];
                xT = xT - CenterT(1);
                yT = yT - CenterT(2);
                %% draw the estimate registation
                if IsPlot
                    close all
                    figure('visible','off');
                    subplot(1, 2, 1);
                    xstep = cumsum(diff(xT))+xT(1);
                    ystep = cumsum(diff(yT))+yT(2);
                    Colors = jet(length(ystep)-1);
                    for i = 1:length(xstep)-1
                        plot(xstep(i:i+1), ystep(i:i+1), 'color', Colors(i, :)); hold on
                    end
                    subplot(1, 2, 2);
                    xstep = cumsum(diff(xT))+xT(1);
                    ystep = cumsum(diff(yT))+yT(1);
                    for i = 1:length(xstep)-1
                        plot(xstep(i:i+1), ystep(i:i+1), 'color', Colors(i, :)); hold on
                    end
                    FigFilNam = sprintf('%s/%s_TranslationEstimation_%d%03d', Day, Topic, Cel, Exp);
                    saveas(gcf,[SaveFolder 'Figures\' FigFilNam '.png']);
                end
                %%                
                for i = 1:NumFrm
                    RefRgt = affine2d([ ...
                        ProgressTransform(1, :, i);...
                        ProgressTransform(2, :, i); ...
                        xT(i) yT(i) 1]);
                    Mov = ImgAtm(:,:,i);
                    ImgAtm(:,:,i) = imwarp(Mov,RefRgt,'OutputView',imref2d(size(Mov)));
                    Mov = I(:,:,i);
                    I(:,:,i) = imwarp(Mov,RefRgt,'OutputView',imref2d(size(Mov)));
                    if mod(i, 20)==0
                        clc
                        fprintf('Files: %d/%d \n', f, nFile);
                        fprintf('\t Step: %d/%d (registration)', 2, 2);
                        fprintf('\t\t Frames: %d/%d \n', i, NumFrm);
                    end
                end
                
                save([SaveFolder 'Translation\' RegistName],'ImgAtm', 'I', 'CenterT', 'ProgressTransform');
                clear ImgAtm I CenterT ProgressTransform
            case 3 % Intensity-Based- sequential
                % Ignore the artifact 
%                 ImgAtm(:, [1:15 end-14:end], :) = nan;
%                 ImgAtm(end-1:end, [1:25 end-24:end], :) = nan;
                ImgAtm = ImgAtm(1:end-2, 16:end-15, :);
                
                % Convert Image to grayscale (0~1)
                RangeImgAtm = range(ImgAtm(:));
                MinImgAtm = min(ImgAtm(:));
                ImgAtm = (ImgAtm-MinImgAtm)/RangeImgAtm;
                Fxd = ImgAtm(:,:,MidFrm);
                
                % Optimizer
                [Opt, Met]             = imregconfig('Multimodal');
                Opt.Epsilon            = 1.5e-6;
                Opt.MaximumIterations  = 300;
%                 optimizer.GrowthFactor = 1.0001;
                Opt.InitialRadius      = 6.25e-4;
                
                TimeIters = nan(NumFrm, 2);
                ProgressTransform = nan(3,3,NumFrm);
                AnchorFrm = find(mod(1:NumFrm, AnchorFrameStep)==0);
                AnchorPair = nan(NumFrm, 1);
                AnchorTransform = nan(3, 3, length(AnchorFrm));
                tic
                for i = 1:NumFrm
                    [~, AnchorId] = min(abs(AnchorFrm-i));
                    AnchorPair(i) = AnchorId;
%                     if ~ismember(i, AnchorFrm)
%                         Fxd = ImgAtm(:,:, AnchorFrm(AnchorId));                        
%                     else
%                         if AnchorId ~= length(AnchorFrm)
%                             Fxd = ImgAtm(:,:, AnchorFrm(AnchorId+1));
%                         else
%                             Fxd = ImgAtm(:,:, end);
%                         end
%                     end
                    Mov = ImgAtm(:,:,i);
                    RefRgt = imregtform(Mov, Fxd, 'translation', Opt, Met);
                    ProgressTransform(:, :, i) = RefRgt.T;
                    if ismember(i, AnchorFrm)
                        MovT = imwarp(Mov,RefRgt,'OutputView',imref2d(size(Mov)));
                        RefRgt = imregtform(Fxd, MovT, 'rigid', Opt, Met);
                        AnchorTransform(:, :, AnchorFrm==i) = RefRgt.T;
                    end
%                     RefRgt = imregtform(Mov, Fxd, 'affine', Opt, Met);
%                     MovImg = I(:,:,i);
%                     I(:,:,i) = imwarp(MovImg,RefRgt,'OutputView',imref2d(size(Fxd)));
%                     ImgAtm(:,:,i) = imwarp(Mov,RefRgt,'OutputView',imref2d(size(Fxd)));
                    if mod(i, 50)==0
                        clc
                        fprintf('Files: %d/%d \n', f, nFile);
                        fprintf('\t Frames: %d/%d \n', i, NumFrm);
                    end
                    TimeIters(i, :) = [i, toc];                    
                end
                keyboard;
                save([SaveFolder RegistName],'ImgAtm', 'I', 'TimeIters');
                clear ImgAtm I TimeIters
            otherwise
                error('Not a defined method for registration');
        end
    end
end