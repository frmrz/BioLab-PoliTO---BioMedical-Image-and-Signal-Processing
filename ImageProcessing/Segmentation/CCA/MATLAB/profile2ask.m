clear 
clc
close all

main = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST';

Profiles_dir = fullfile('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT3-INTERPOLATED');
Images_dir = fullfile('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/IMAGES');

savePath = fullfile('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/MASKS-GT3');

imgs = dir(Images_dir);

 
for i = 1:size(imgs,1)
     i
    if ~imgs(i).isdir && ~strcmp(imgs(i).name,'.DS_Store') 
        try 
            filename = strtok(imgs(i).name,'.'); 
                    
            %  profiles
            LI = load(fullfile(Profiles_dir,[filename '-LI.txt']));
            MA = load(fullfile(Profiles_dir,[filename '-MA.txt']));
            
            img = imread(fullfile(Images_dir,[filename '.png']));
                    
    %         Mask = uint8(roipoly(img,[LI(:,1),flip(MA(:,1))],[LI(:,2),flip(MA(:,2))])).*255;
            Mask = uint8(roipoly(img,[LI(1,:),flip(MA(1,:))],[LI(2,:),flip(MA(2,:))])).*255;
            
            imwrite(Mask,fullfile(savePath,[filename '.png']));
        catch
            fprintf('\n Error in image %s', filename)
        end
    end
end