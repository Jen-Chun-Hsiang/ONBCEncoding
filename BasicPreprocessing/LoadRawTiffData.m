clear; close all; clc;

%% Analysis Parameters
Day = '111522';
Deinterleave = 4;
ChGFP = 1;
Rig = 2; % 1: old scientifica, 2: hyperscope
%% Load Data
PathName = uigetdir(pwd, 'Select a folder');
FileNames = dir(fullfile(PathName, '*.tif'));
FileIds = find(contains({FileNames.name},Day));
nFile = length(FileIds);

SaveFolder = strsplit(PathName, 'RawData');
SaveFolder = [SaveFolder{1} '\MatFile\'];
PathName = [PathName '\'];
for f = 1:nFile
    FileName = FileNames(FileIds(f)).name;
    Words = strsplit(FileName, ['_' Day '_']);
    Topic = Words{1};
    switch Rig
        case 1
            Words = strsplit(Words{2}, '.');
            Cel = floor(str2double(Words{1})/1000);
            Exp = mod(str2double(Words{1}), 1000);
            HandleDuplicate = 1;
        case 2
            Words = strsplit(Words{2}, '_');
            Cel = str2double(Words{1});
            Words = strsplit(Words{2}, '.');
            Exp = str2double(Words{1});
            HandleDuplicate = 0;
    end
    
    SaveFileName = sprintf('%s_%s_%d%03d_%d',Topic, Day, Cel, Exp, Deinterleave);
    if exist([SaveFolder SaveFileName '.mat'], 'file')
        disp('Tiff data has been loaded');
    else
        Info = imfinfo([PathName FileName]);
        NumFrm = numel(Info);
        Pixel_x = Info(1).Height;
        Pixel_y = Info(1).Width;
        if NumFrm == 25600
            HandleDuplicate =0;
        end
        if HandleDuplicate
            switch Pixel_y
                case 256
                    switch Pixel_x
                        case 40
                            DuplicaPos = 33:40;
                            CopyPos = 9:16;
                            FinalPixel_x = 1:32;
                        case 80
                            DuplicaPos = 65:80;
                            CopyPos = 49:64;
                            FinalPixel_x = 1:64;
                    end
                case 128
                    switch Pixel_x
                        case 80
                            DuplicaPos = 65:80;
                            CopyPos = 17:32;
                            FinalPixel_x = 1:64;
                    end
                case 64
                    HandleDuplicate = 0;
            end
        end
        if mod(NumFrm, Deinterleave) ~= 0
            error('The number of deinterleave is not matched to the overall frame!');
        else
            NumFrm = NumFrm / Deinterleave;
        end
        clear Info
        disp('Loading raw data from tiff file...');
        
        I_GCh = nan(Pixel_x, Pixel_y, NumFrm);
        I_RCh = nan(Pixel_x, Pixel_y, NumFrm);
        I_Struct = nan(Pixel_x, Pixel_y, NumFrm);
        if Deinterleave == 4
            I_Stim = nan(Pixel_x, Pixel_y, NumFrm);
        end
        TifLink = Tiff([PathName FileName], 'r');
        for i=1:NumFrm
            TifLink.setDirectory((i-1)*Deinterleave+1);
            I_GCh(:,:,i) = TifLink.read();
            TifLink.setDirectory((i-1)*Deinterleave+2);
            I_RCh(:,:,i) = TifLink.read();
            TifLink.setDirectory((i-1)*Deinterleave+3);
            I_Struct(:,:,i) = TifLink.read();
            
            if Deinterleave == 4
                TifLink.setDirectory((i-1)*Deinterleave+4);
                I_Stim(:,:,i) = TifLink.read();
            end
            if mod(i, 20) == 0
                clc
                fprintf('File: %d/%d \n', f, nFile);
                fprintf('\t Frame: %d/%d \n', i, NumFrm);
            end
        end
        TifLink.close();
        
        switch ChGFP
            case 1
                I = I_GCh;
            case 2
                I = I_RCh;
        end
        if HandleDuplicate
            I(CopyPos, :, :) = (I(CopyPos, :, :)+I(DuplicaPos, :, :))/2;
            I = I(FinalPixel_x, :, :);
        end
        SaveFileName = sprintf('%s_%s_%d%03d_%d',Topic, Day, Cel, Exp, 1);
        save([SaveFolder SaveFileName '.mat'], 'I');
        
        %         I = I_RCh;
        %         SaveFileName = sprintf('%s_%s_%d%03d_%d',Topic, Day, Cel, Exp, 2);
        %         save([SaveFolder SaveFileName '.mat'], 'I');
        
        I = I_Struct;
        if HandleDuplicate
            I(CopyPos, :, :) = (I(CopyPos, :, :)+I(DuplicaPos, :, :))/2;
            I = I(FinalPixel_x, :, :);
        end
        SaveFileName = sprintf('%s_%s_%d%03d_%d',Topic, Day, Cel, Exp, 3);
        save([SaveFolder SaveFileName '.mat'], 'I');
        
        if Deinterleave == 4
            I = I_Stim;
            if HandleDuplicate
                I(CopyPos, :, :) = (I(CopyPos, :, :)+I(DuplicaPos, :, :))/2;
                I = I(FinalPixel_x, :, :);
            end
            SaveFileName = sprintf('%s_%s_%d%03d_%d',Topic, Day, Cel, Exp, 4);
            save([SaveFolder SaveFileName '.mat'], 'I');
            clear I_Stim
        end
        
        disp('Finished!');
    end
end