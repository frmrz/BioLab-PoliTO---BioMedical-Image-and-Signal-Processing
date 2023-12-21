%% script for plotting profiles on images

close all
clear
clc

mainDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/';

imgDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","IMAGES");
profilesDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","PROFILES");
gtFolder = fullfile(profilesDir,"GT-CONSENSUS");

folders = dir(profilesDir);
imgs = dir(imgDir);

toPlot = [22:24,26];

subPos = [[1,2,5,6];[3,4,7,8];[9,10,13,14];[11,12,15,16]];

for jj = 3 : length(imgs)

    filename = extractBefore(imgs(jj).name,'.png');
    img = imread(fullfile(imgDir,imgs(jj).name));

    h=figure;
    h.Position = [ 0 0 2000 2000];
    h.Visible = 'off';

    for ii = 1 : length(toPlot)
        idz = toPlot(ii);

        tempFolder = fullfile(profilesDir,folders(idz).name);

        try
            LI = TurnColumn(load(fullfile(tempFolder,strrep(imgs(jj).name,'.png','-LI.txt'))));
            MA = TurnColumn(load(fullfile(tempFolder,strrep(imgs(jj).name,'.png','-MA.txt'))));
    
            LI_GT = TurnColumn(load(fullfile(gtFolder,strrep(imgs(jj).name,'.png','-LI.txt'))));
            MA_GT = TurnColumn(load(fullfile(gtFolder,strrep(imgs(jj).name,'.png','-MA.txt'))));
    
%             subplot(4,4,subPos(ii,:)),
            imshow(img),hold on,
%             title(sprintf('%s',folders(idz).name));
    
%             plot(LI(:,1),LI(:,2),'r','LineWidth',1),drawnow,hold on,
%             plot(MA(:,1),MA(:,2),'g','LineWidth',1),drawnow,hold on,
    
            plot(LI_GT(:,1),LI_GT(:,2),'c','LineWidth',1),drawnow,hold on,
            plot(MA_GT(:,1),MA_GT(:,2),'y','LineWidth',1),drawnow,hold on
        catch
        end

    end

    if exist('h','var')
        saveas(h,fullfile('/media/francesco/DEV001/PROJECT-CUBS-SPIE//RESULTS/MATLAB/PROFILE_IUS23',imgs(jj).name));
    end
       
    close all
    clear h
end