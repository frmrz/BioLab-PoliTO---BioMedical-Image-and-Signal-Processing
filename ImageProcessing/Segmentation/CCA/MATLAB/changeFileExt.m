clear

imgs = dir(cd)

for ii = 3 : length(imgs)
    filename= imgs(ii).name;
    filename = filename(1:end-4);
    img = imread(fullfile(imgs(ii).folder,imgs(ii).name));
    imwrite(img,fullfile(imgs(ii).folder,[filename '.png']))
end