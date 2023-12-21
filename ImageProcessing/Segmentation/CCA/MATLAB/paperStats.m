% this script is a utility script for retrieving the stats for the paper

clear  
close all

%% define images and GT directories


Data_fd = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/';

Train_fd = fullfile(Data_fd,'FULL');
imgs_fd = fullfile(Train_fd,"IMAGES/");
CF_fd  = fullfile(Train_fd,'CF');

imgs = dir(imgs_fd);
cfs = dir(CF_fd);

for i = 3 : length(cfs)
    i
    CFs(i-2,1)=load(fullfile(cfs(i).folder,cfs(i).name));
    IMGs = imread(fullfile(imgs(i).folder,imgs(i).name));

    if size(IMGs,3)>1
        [row,col] = size(rgb2gray(IMGs));
    else 
        [row,col] = size(IMGs);
    end

    rr(i-2,1)= col/row;
end

mean(CFs)
std(CFs)

mean(rr)
std(rr)

%% test
load('/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/interOP.mat')
Test_fd = fullfile(Data_fd,'TEST');
imgs_fd = fullfile(Test_fd,"IMAGES/");
CF_fd  = fullfile(Test_fd,'CF');

files = dir(imgs_fd);
q = 1;

fxs = table2struct(gt2Table(:,1));
for i = 1 :length(fxs)
    fxs(i).name = [fxs(i).name '.png'];
end
    
flag =0 ;
for i = 3 : length(fxs)
    i
    try
        CFs(i-2,1)=load(fullfile(CF_fd,strrep(fxs(i).name,'.png','_CF.txt')));
    catch
        CFs(i-2,1) = 0.06;
        flag = flag+1;
    end
    IMGs = imread(fullfile(imgs_fd,fxs(i).name));

    if size(IMGs,3)>1
        [row,col] = size(rgb2gray(IMGs));
    else 
        [row,col] = size(IMGs);
    end

    rr(i-2,1)= col/row;
end

mean(CFs)
std(CFs)

mean(rr)
std(rr)