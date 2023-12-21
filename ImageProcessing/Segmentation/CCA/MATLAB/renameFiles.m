close all
clear

% mainDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES';
mainDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DATASET_Multicenter/DATABASE-HongKong/IMAGES';
folders = dir(mainDir);
% for jj = 3 : length(folders)
    % Get all text files in the current folder
%     files = dir(fullfile(mainDir,folders(jj).name));
    files = dir(mainDir);
    % Loop through each file 
    for id = 3:length(files)
        % Get the file name 
        [a, f,ext] = fileparts(files(id).name);
        rename = [strrep(f,' ','') '.bmp']; 
        movefile(fullfile(mainDir,files(id).name), fullfile(mainDir,rename)); 
%         if contains(f,'fig0')
%     %         rename = [strrep(f,' ','') '.png']; 
% %             rename = [strrep(f,'GT','') ext]; 
% %             movefile(files(id).name, rename); 
%             delete(fullfile(files(id).folder,files(id).name))
%         end
    end
% end