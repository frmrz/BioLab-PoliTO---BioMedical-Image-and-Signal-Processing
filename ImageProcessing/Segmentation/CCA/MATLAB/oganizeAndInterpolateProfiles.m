close all
clear
clc

%% script for the organization and interpolation of profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% part 1: organize and interpolate manual profiles from multi-op
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% folders

dataDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA';
centersDir = fullfile(dataDir,'DATASET_Multicenter/') ;

outProfilesDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT1-ORIGINAL';
outInterpProfilesDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT1-INTERPOLATED';

outProfilesDirBis = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT1B-ORIGINAL';
outInterpProfilesDirBis = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT1B-INTERPOLATED';

outProfilesDir2 = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT2-ORIGINAL';
outInterpProfilesDir2 = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT2-INTERPOLATED';

outProfilesDir3 = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT3-ORIGINAL';
outInterpProfilesDir3 = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT3-INTERPOLATED';

outputMerge = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT-CONSENSUS';
%% files

centers = dir(centersDir);

%% loop across centers
originalProfiles = struct('name',[],'center',[],'GT1_LI',[],'GT1_MA',[], ...
    'GT1B_LI',[],'GT1B_MA',[], ...
    'GT2_LI',[],'GT2_MA',[], ...
    'GT3_LI',[],'GT3_MA',[]);
interpProfiles = struct('name',[],'center',[],'GT1_LI',[],'GT1_MA',[], ...
    'GT1B_LI',[],'GT1B_MA',[], ...
    'GT2_LI',[],'GT2_MA',[], ...
    'GT3_LI',[],'GT3_MA',[], ...
    'cons_LI',[],'cons_MA',[]);
q=1;
% figure,

for ii = 3 : length(centers)
    fprintf('%d\n',ii);
    tempCenter = fullfile(centersDir,centers(ii).name);

    tempCFdir = fullfile(tempCenter,"CF/");

    try
        tempGT1dir = fullfile(tempCenter,"GT1-LIMA-Profiles/");
        tempGT1Bdir = fullfile(tempCenter,"GT1bis-LIMA-Profiles/");
    
        CFs = dir(tempCFdir);
        GT1 = dir(tempGT1dir);
        GT1B = dir(tempGT1Bdir);
        imgs = dir(fullfile(tempCenter,"IMAGES"));

    
        tempGT2dir = fullfile(tempCenter,"GT2-LIMA-Profiles/");
        tempGT3dir = fullfile(tempCenter,"GT3-LIMA-Profiles/");
    
        GT2 = dir(tempGT2dir);
        GT3 = dir(tempGT3dir);
    catch
        fprintf('No GT2/GT3 directories')
    end

    for jj = 3 : length(imgs)

        filename = imgs(jj).name(1:end-4)

        if ~contains(filename,'Thumb') && ~contains(filename,'fig0')

            img = imread(fullfile(imgs(jj).folder,imgs(jj).name));

            try
                tempCF = load(fullfile(tempCFdir,CFs(jj).name));
            catch
                tempCF = 0.06;
            end
    
            % load profiles
            GT1_LI = []; GT1_MA = [];
            GT1B_LI = []; GT1B_MA = [];
            GT2_LI = []; GT2_MA = [];
            GT3_LI = []; GT3_MA = [];

            GT1_LI_interp = []; GT1_MA_interp = [];
            GT1B_LI_interp = []; GT1B_MA_interp = [];
            GT2_LI_interp = []; GT2_MA_interp = [];
            GT3_LI_interp = []; GT3_MA_interp = [];

            %% GT1
            try
                GT1_LI = load(fullfile(tempGT1dir,[filename '-GTLI.txt']));
                GT1_MA = load(fullfile(tempGT1dir,[filename '-GTMA.txt']));
                [GT1_LI_interp, GT1_MA_interp] = LI_MA_interp(GT1_LI,GT1_MA);
    
                if mean(GT1_MA_interp(2,:)) < mean(GT1_LI_interp(2,:))
                    tempSwap = GT1_LI_interp;
                    GT1_LI_interp = GT1_MA_interp;
                    GT1_MA_interp = tempSwap;
                end
            
            catch
                fprintf('No GT1 profile in file %s\n',filename)
            end
    
            %% GT1B
            try
                GT1B_LI = load(fullfile(tempGT1Bdir,[filename '-GTLI.txt']));
                GT1B_MA = load(fullfile(tempGT1Bdir,[filename '-GTMA.txt']));
                [GT1B_LI_interp, GT1B_MA_interp] = LI_MA_interp(GT1B_LI,GT1B_MA);
    
                if mean(GT1B_MA_interp(2,:)) < mean(GT1B_LI_interp(2,:))
                    tempSwap = GT1B_LI_interp;
                    GT1B_LI_interp = GT1B_MA_interp;
                    GT1B_MA_interp = tempSwap;
                end
            
            catch
                fprintf('No GT1B profile in file %s\n',filename)
            end
            
                    %% GT2
            try
                GT2_LI = load(fullfile(tempGT2dir,[filename '-GTLI.txt']));
                GT2_MA = load(fullfile(tempGT2dir,[filename '-GTMA.txt']));
                [GT2_LI_interp, GT2_MA_interp] = LI_MA_interp(GT2_LI,GT2_MA);

                if mean(GT2_MA_interp(2,:)) < mean(GT2_LI_interp(2,:))
                    tempSwap = GT2_LI_interp;
                    GT2_LI_interp = GT2_MA_interp;
                    GT2_MA_interp = tempSwap;
                end
                    
            catch
                fprintf('No GT2 profile in file %s\n',filename)
            end
            
                    %% GT3
            try
                GT3_LI = load(fullfile(tempGT3dir,[filename '-GTLI.txt']));
                GT3_MA = load(fullfile(tempGT3dir,[filename '-GTMA.txt']));
                [GT3_LI_interp, GT3_MA_interp] = LI_MA_interp(GT3_LI,GT3_MA);
    
                if mean(GT3_MA_interp(2,:)) < mean(GT3_LI_interp(2,:))
                    tempSwap = GT3_LI_interp;
                    GT3_LI_interp = GT3_MA_interp;
                    GT3_MA_interp = tempSwap;
                end
            catch
                fprintf('No GT3 profile in file %s\n',filename)
            end
           
                
            %% make consensus
            [cons_LI,cons_MA,cvLI,cvMA,flag_high_bias] = consensusGT(img,GT1_LI_interp,GT1_MA_interp,...
                GT1B_LI_interp,GT1B_MA_interp,...
                GT2_LI_interp,GT2_MA_interp,...
                GT3_LI_interp,GT3_MA_interp);
            
            %% save profiles
            % write struct
    
            temp.name = filename;
            temp.center = centers(ii).name;
    
            temp.GT1_LI = GT1_LI;        temp.GT1_MA = GT1_MA;
            temp.GT1B_LI = GT1B_LI;      temp.GT1B_MA = GT1B_MA;
            temp.GT2_LI = GT2_LI;      temp.GT2_MA = GT2_MA;
            temp.GT3_LI = GT3_LI;      temp.GT3_MA = GT3_MA;
    
            originalProfiles(q) = temp;
    
            clear temp
    
            temp.name = filename;
            temp.center = centers(ii).name;
    
            temp.GT1_LI = GT1_LI_interp;        temp.GT1_MA = GT1_MA_interp;

            if ~isempty(GT1B_LI_interp) && ~isempty(GT1_LI_interp) && sum(GT1B_LI_interp(2,:))-sum(GT1_LI_interp(2,:))~=0
                temp.GT1B_LI = GT1B_LI_interp;      temp.GT1B_MA = GT1B_MA_interp;
            else 
                temp.GT1B_LI = []; temp.GT1B_MA = [];
            end
            
            temp.GT2_LI = GT2_LI_interp;      temp.GT2_MA = GT2_MA_interp;
            temp.GT3_LI = GT3_LI_interp;      temp.GT3_MA = GT3_MA_interp;
            temp.cons_LI = cons_LI;      temp.cons_MA = cons_MA;

            interpProfiles(q) = temp;
            
            q=q+1;

            clear temp 
            clear GT1 GT1_MA GT1B_LI GT1B_MA GT2_LI GT2_MA GT3_LI GT3_MA
            clear GT1_LI_interp GT1_MA_interp GT1B_LI_interp GT1B_MA_interp GT2_LI_interp GT2_MA_interp GT3_LI_interp GT3_MA_interp

            GT1_LI = []; GT1_MA = [];
            GT1B_LI = []; GT1B_MA = [];
            GT2_LI = []; GT2_MA = [];
            GT3_LI = []; GT3_MA = [];

            GT1_LI_interp = []; GT1_MA_interp = [];
            GT1B_LI_interp = []; GT1B_MA_interp = [];
            GT2_LI_interp = []; GT2_MA_interp = [];
            GT3_LI_interp = []; GT3_MA_interp = [];

            %% check visually
            
%             Msk = uint8(roipoly(img,[cons_LI(1,:),flip(cons_MA(1,:))], ...
%                 [cons_LI(2,:),flip(cons_MA(2,:))])).*255;
%             ov = labeloverlay(img,Msk,"Transparency",0.8);
%             imshow(ov),hold on
%             imshow(img),hold on
% 
%             try
%             plot(GT1_LI_interp(1,:),GT1_LI_interp(2,:),'r','LineWidth',2.5),drawnow,hold on,
%             plot(GT1_MA_interp(1,:),GT1_MA_interp(2,:),'--r','LineWidth',2.5),drawnow,hold on,
%             catch
%             end
%             try
%             plot(GT1B_LI_interp(1,:),GT1B_LI_interp(2,:),'g','LineWidth',2),drawnow,hold on,
%             plot(GT1B_MA_interp(1,:),GT1B_MA_interp(2,:),'--g','LineWidth',2),drawnow,hold on,
%             catch
%             end
%             try
%             plot(GT2_LI_interp(1,:),GT2_LI_interp(2,:),'c','LineWidth',1.5),drawnow,hold on,
%             plot(GT2_MA_interp(1,:),GT2_MA_interp(2,:),'--c','LineWidth',1.5),drawnow,hold on,
%             catch
%             end
%             try
%             plot(GT3_LI_interp(1,:),GT3_LI_interp(2,:),'m','LineWidth',1),drawnow,hold on,
%             plot(GT3_MA_interp(1,:),GT3_MA_interp(2,:),'--m','LineWidth',1),drawnow,hold on,
%             catch
%             end
%             try
%             plot(cons_LI(1,:),cons_LI(2,:),'y','LineWidth',0.5),drawnow,hold on,
%             plot(cons_MA(1,:),cons_MA(2,:),'--y','LineWidth',0.5),drawnow,hold on,
%             catch
%             end
% 
% LI_A1 = GT1_LI_interp
% MA_A1 = GT1_MA_interp
% LI_A1B = GT1B_LI_interp
% MA_A1B = GT1B_MA_interp
% LI_A2 = GT2_LI_interp
% MA_A2 = GT2_MA_interp
% LI_A3 = GT3_LI_interp
% MA_A3 = GT3_MA_interp
            %% write files
        
%             write_txt_file(GT1_LI,[filename '-GTLI.txt'],outProfilesDir)
%             write_txt_file(GT1_MA,[filename '-GTMA.txt'],outProfilesDir)
    
%             write_txt_file(GT1_LI_interp,[filename '-GTLI.txt'],outInterpProfilesDir)
%             write_txt_file(GT1_MA_interp,[filename '-GTMA.txt'],outInterpProfilesDir)

%             write_txt_file(GT1B_LI,[filename '-LI.txt'],outProfilesDirBis)
%             write_txt_file(GT1B_MA,[filename '-MA.txt'],outProfilesDirBis)
%     
%             write_txt_file(GT1B_LI_interp,[filename '-LI.txt'],outInterpProfilesDirBis)
%             write_txt_file(GT1B_MA_interp,[filename '-MA.txt'],outInterpProfilesDirBis)

%             write_txt_file(GT2_LI,[filename '-LI.txt'],outProfilesDir2)
%             write_txt_file(GT2_MA,[filename '-MA.txt'],outProfilesDir2)
    
%             write_txt_file(GT2_LI_interp,[filename '-LI.txt'],outInterpProfilesDir2)
%             write_txt_file(GT2_MA_interp,[filename '-MA.txt'],outInterpProfilesDir2)

%             write_txt_file(GT3_LI,[filename '-LI.txt'],outProfilesDir3)
%             write_txt_file(GT3_MA,[filename '-MA.txt'],outProfilesDir3)
    
%             write_txt_file(GT3_LI_interp,[filename '-LI.txt'],outInterpProfilesDir3)
%             write_txt_file(GT3_MA_interp,[filename '-MA.txt'],outInterpProfilesDir3)
% 
%             if ~(mean(cvLI(isfinite(cvLI))) < 0.25 && mean(cvMA(isfinite(cvMA))) < 0.25 || flag_high_bias == 1)
%                 write_txt_file(cons_LI,[filename '-LI.txt'],outputMerge);
%                 write_txt_file(cons_MA,[filename '-MA.txt'],outputMerge);
%             end
        end
    end
end

% count images with all profiles
tableProfiles = struct2table(Struct_Empty_To_Nan(interpProfiles));

consTable = table();
gt1bbTable = table();
gt2Table = table();

q = 1;
q1 = 1;
q2 = 1;
q3 = 1;

for ii = 1 : height(tableProfiles)
    % count consensus
    temp = tableProfiles{ii,11};
    if ~isnan(sum(temp{1,1},"all"))
        consTable(q,:) = tableProfiles(ii,:);
        q=q+1;
    end
    % count GT1b
    temp = tableProfiles{ii,5};
    if ~isnan(sum(temp{1,1},"all"))
        gt1bbTable(q1,:) = tableProfiles(ii,:);
        q1=q1+1;
    end
    % count GT2 & GT3
    temp1 = tableProfiles{ii,7};
    temp2 = tableProfiles{ii,9};
    if ~isnan(sum(temp1{1,1},"all")) && ~isnan(sum(temp2{1,1},"all"))
        gt2Table(q3,:) = tableProfiles(ii,:);
        q3=q3+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% part 2: extract and organize profiles from unet masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% resultsDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/RESULTS/PYTHON';
% originalImagesDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/IMAGES';
% CFDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/CF';
% profilesDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES';
% profilesOvDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES-OVERLAY';
% profilesDEbugDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES-DEBUG';
% experiments = dir(resultsDir);
% cfs = dir(CFDir);
% 
% set(0,'DefaultFigureVisible','on');
% 
% % loop across experiments
% for ii = 17 : length(experiments)
%     tic
%     experiments(ii).name
%     resizedOutputDir = fullfile(resultsDir,experiments(ii).name,"Test_output");
%     outResMasks = dir(resizedOutputDir);
% 
% %     if isempty(outResMasks)
% %         resizedOutputDir = fullfile(resultsDir,experiments(ii).name,"Test_output");
% %         outResMasks = dir(resizedOutputDir);
% %     end
% 
%     destination = fullfile(profilesDir,experiments(ii).name);
% 
% %     if ~isfolder(destination)
% %         mkdir(destination);
% 
%         % loop across images
%         for jj = 1698 : length(outResMasks)
%             filename = outResMasks(jj).name(1:end-4),
%     
%             msk = imread(fullfile(resizedOutputDir,outResMasks(jj).name));
%             imgOriginal = imread(fullfile(originalImagesDir,outResMasks(jj).name));
%     
%             if size(imgOriginal,3) == 3
%                 imgOriginal = rgb2gray(imgOriginal);
%             end
%     
%             try
%                 CF = load(fullfile(CFDir,[filename '_CF.txt']));    
%             catch
%                 CF = 0.06;
%             end
%     
%             %% find image area and resize parameters
%             [row,col] = size(imgOriginal);
%     
%             try
%                 [Bounds,Img_cropped,~,~]=find_US_Image_area(imgOriginal);    % first crop based on info
%                 [row_c,col_c] = size(Img_cropped);
%             catch
%                 fprintf(' Error in resizing image %s, cropping with alternative function\n',filename);
%                 [Img_cropped,~,minX,minY] = fm_autocrop(imgOriginal,[]);
%                 Bounds(1) = minY;
%                 Bounds(2) = minY + size(Img_cropped,1)-1;
%                 Bounds(3) = minX;
%                 Bounds(4) = minX + size(Img_cropped,2)-1;
%                 [row_c,col_c] = size(Img_cropped);
%             end
%     
%             row_final = Bounds(2)-Bounds(1);
%             col_final = Bounds(4)-Bounds(3);
%     
%             rr = CF/0.0747;
%             First_crop = Img_cropped;
%             Img_cropped = imresize(Img_cropped,[row_c*rr,col_c*rr]); % resize to standard resolution
%             [row_r,col_r] = size(Img_cropped);
%     
%             % find ratio of i-th image
%             ratio=col_r/row_r;
%             to_mean = 1 - ratio;
%     
%             % padding to square
%             if to_mean > 0
%                 col_pad = round(row_r*1);
%                 col_margin = round((col_pad-col_r)/2);
%                 
%                 Img_padded = padarray(Img_cropped,[0 col_margin],0,'both');
%             elseif to_mean < 0
%                 row_pad = round(col_r/1);
%                 row_margin = round((row_pad-row_r)/2);
%                 
%                 Img_padded = padarray(Img_cropped,[row_margin 0],0,'both');
%             else 
%                 Img_padded = Img_cropped;
%             end
%             
%             [row_final1,col_final1]=size(Img_padded); % size of this image would be resized to 
%             ratio_final=col_final1/row_final1;      % 480x480
%             
%             %% resize OUTPUT mask accordingly to parameters
%             % resize to ratio final
%             out_last = imresize(msk,[row_final1,col_final1]);
%     
%             fig_debug = figure, hold on,
%             ov_ratioFinal_check = labeloverlay(Img_padded,uint8(out_last>127.5));
%             subplot(2,2,1),imshow(ov_ratioFinal_check),title('Last resize');
%     
%             % crop to "padding to square"
%             [rf,cf]=size(out_last);
%     
%             if exist('col_margin')
%                 out_crop = out_last(:,col_margin:cf-col_margin-1);
%             elseif exist('row_margin')
%                 out_crop = out_last(row_margin:rf-row_margin-1,:);
%             else
%                 out_crop = out_last;
%             end
%     
%             ov_ratioFinal_check = labeloverlay(Img_cropped,uint8(out_crop>127.5));
%             subplot(2,2,2),imshow(ov_ratioFinal_check),title('First pad');
%     
%             % resize to original resolution
%             out_roi = imresize(out_crop,[row_c,col_c]);
%     
%             ov_ratioFinal_check = labeloverlay(First_crop,uint8(out_roi>127.5));
%             subplot(2,2,3),imshow(ov_ratioFinal_check),title('First resize');
%     
%             % pad to original size
%             OUT_final  = padarray(out_roi,Bounds(1)-1,'pre');     OUT_final  = padarray(OUT_final,row-Bounds(2),'post');
%             OUT_final  = padarray(OUT_final',Bounds(3)-1,'pre')';   OUT_final  = padarray(OUT_final',col-Bounds(4),'post')';
%             
%             if size(OUT_final) ~= size(imgOriginal)
%                 fprintf(' Error in resizing image %\n',filename);
%                 OUT_final = imresize(OUT_final,[row,col]);
%             end
%                     
%             ov_ratioFinal_check = labeloverlay(imgOriginal,uint8(OUT_final>127.5));
%             subplot(2,2,4),imshow(ov_ratioFinal_check),title('Original');
%             saveas(fig_debug,fullfile(profilesDEbugDir,experiments(ii).name,[filename '.jpg']));
%     
%             %% visualize resizing result
%     %         ov = labeloverlay(imgOriginal,logical(OUT_final),"Transparency",0.7);
%     %         imshow(ov);
%     
%             clear Img_cropped Img_padded out_roi out_crop out_last row_margin col_margin row_pad col_pad
%             %% extract profiles from mask
%             try
%                 OUT_final = OUT_final>255*0.5; st = regionprops(OUT_final,'Area'); if ~isempty(st) ; OUT_final = imfill(bwareaopen(OUT_final,max([st(:).Area])-1),'holes');end
%                 [LI_AUTO,MA_AUTO]=getLIMAfromMask(OUT_final);
%         
% %                 destination = fullfile(profilesDir,experiments(ii).name);
%         
% %                 if ~isfolder(destination)
% %                     mkdir(destination);
% %                 end
%         
%     %             h=figure("Visible","off");
%     %             imshow(imgOriginal),hold on,
%     %             plot(LI_AUTO(:,1),LI_AUTO(:,2),'r','LineWidth',1),drawnow,hold on,
%     %             plot(MA_AUTO(:,1),MA_AUTO(:,2),'g','LineWidth',1),drawnow,hold on,
%     %             saveas(h,fullfile(profilesOvDir,experiments(ii).name,[filename '.jpg']));
%     %             close all
%     
% %                 write_txt_file(LI_AUTO,[filename,'-LI.txt'],destination);
% %                 write_txt_file(MA_AUTO,[filename,'-MA.txt'],destination);
%     
%                 clear LI_AUTO MA_AUTO filename OUT_final
%             catch
%                 if sum(OUT_final,'all')==0
%                     fprintf('Problem in image %s of experiment %s : no object\n',filename,experiments(ii).name);
%                     clear filenname OUT_final
%                 else
%                     fprintf('Problem in image %s of experiment %s : other problem\n',filename,experiments(ii).name);
%                     clear filenname OUT_final
%                 end
%             end
%     
%         end
%         toc
% %     end
% end
