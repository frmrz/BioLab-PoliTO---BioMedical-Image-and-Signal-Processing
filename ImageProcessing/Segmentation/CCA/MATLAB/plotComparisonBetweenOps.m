%% script for plotting profiles on images

close all
clear
clc

mainDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE';

imgDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","IMAGES");
profilesDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","PROFILES");
gtFolder = fullfile(profilesDir,"GT-CONSENSUS");

folders = dir(profilesDir);
imgs = dir(imgDir);

toPlot = [4,8,10,26];

subPos = [[1,2,5,6];[3,4,7,8];[9,10,13,14];[11,12,15,16]];

for jj = 14 : length(imgs)

    filename = extractBefore(imgs(jj).name,'.png')
    img = imread(fullfile(imgDir,imgs(jj).name));

    h=figure;
    h.Position = [ 0 0 1800 900];
    h.Visible = 'on';

    for ii = 1 : length(toPlot)
        idz = toPlot(ii);

        tempFolder = fullfile(profilesDir,folders(idz).name);

        try
            LI = TurnColumn(load(fullfile(tempFolder,strrep(imgs(jj).name,'.png','-LI.txt'))));
            MA = TurnColumn(load(fullfile(tempFolder,strrep(imgs(jj).name,'.png','-MA.txt'))));
    
            LI_GT = TurnColumn(load(fullfile(gtFolder,strrep(imgs(jj).name,'.png','-LI.txt'))));
            MA_GT = TurnColumn(load(fullfile(gtFolder,strrep(imgs(jj).name,'.png','-MA.txt'))));
    
%             subplot(4,4,subPos(ii,:))
            imshow(img),hold on,
%             title(sprintf('%s',folders(idz).name));
    
            plot(LI(:,1),LI(:,2),'c','LineWidth',3),drawnow,hold on,
            plot(MA(:,1),MA(:,2),'y','LineWidth',3),drawnow,hold on,

%             saveas(h, '/media/francesco/DEV001/PROJECT-CUBS-SPIE/PAPER/FIGURES/FIG1-poster/GT3.png');
            export_fig(h,'/media/francesco/DEV001/PROJECT-CUBS-SPIE/PAPER/FIGURES/FIG1-poster/GT3.png','-png','-r300');
%             plot(LI_GT(:,1),LI_GT(:,2),'c','LineWidth',0.5),drawnow,hold on,
%             plot(MA_GT(:,1),MA_GT(:,2),'y','LineWidth',0.5),drawnow,hold on
hold off
ii=3
        catch
            fprintf('Error for file %s', filename);
        end

    end

%     if exist('h','var')
%         saveas(h,fullfile('/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/RESULTS/MATLAB/PROFILE-OVERLAY-v1',imgs(jj).name));
%     end
       
    close all
    clear h
end