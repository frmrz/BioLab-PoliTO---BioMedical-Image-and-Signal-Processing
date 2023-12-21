% FILEPATH: /media/francesco/DEV001/CENTRAL_PROJECTS_REPOSITORY/BioLab-PoliTO---BioMedical-Image-and-Signal-Processing/ImageProcessing/Segmentation/CCA/MATLAB/profile2Mask.m
% This script converts profiles into masks for a set of images.
% It reads profiles from a directory, loads corresponding images, and generates masks based on the profiles.
% The masks are then saved in a specified directory.
% 
% Inputs:
% - Profiles_dir: Directory path where the profiles are stored.
% - Images_dir: Directory path where the images are stored.
% - savePath: Directory path where the generated masks will be saved.
% 
% Output:
% - None
% 
% Example usage:
% Profiles_dir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT3-INTERPOLATED';
% Images_dir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/IMAGES';
% savePath = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/MASKS-GT3';
% profile2Mask(Profiles_dir, Images_dir, savePath);
%
% Note: This script assumes that the profiles and images have specific naming conventions.
% It expects the profiles to be in the format '<filename>-LI.txt' and '<filename>-MA.txt',
% and the images to be in the format '<filename>.png'.
% Any exceptions or errors encountered during the process will be displayed in the command window.

clear 
clc
close all

main = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST';

Profiles_dir = fullfile('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT3-INTERPOLATED');
Images_dir = fullfile('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/IMAGES');

savePath = fullfile('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/MASKS-GT3');

imgs = dir(Images_dir);

 
for i = 1:size(imgs,1)
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