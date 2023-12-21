close all
clear

main = '/media/francesco/FMZ_archive/POLITO/PROJECT-CUBS-SPIE';

imgs_fd = fullfile(main,'DATA','ORIGINAL','images');
imgs_files = del_junk(dir(imgs_fd));

cf_fd = fullfile(main,'DATA','CF');

Profiles_dir_A1 = fullfile(main,'DATA','ORIGINAL','PROFILES','Manual-A1');
Profiles_dir_A1s = fullfile(main,'DATA','ORIGINAL','PROFILES','Manual-A1s');
Profiles_dir_TUM = fullfile(main,'DATA','ORIGINAL','PROFILES','Computerized-TUM_DE');

for jj = 1:length(imgs_files)

    filename = strtok(imgs_files(jj).name,'.'); 

    img = imread(fullfile(imgs_files(jj).folder,imgs_files(jj).name));
    CF = load(fullfile(cf_fd,[filename '_CF.txt']));

    %%  Get SNR from GT
    li = load(fullfile(Profiles_dir_A1,[filename '-LI.txt']));
    ma = load(fullfile(Profiles_dir_A1,[filename '-MA.txt']));

    li(2,:) = round(li(2,:));
    ma(2,:) = round(ma(2,:));

    thick = mean(ma(:,2)-li(:,2))*1.5;

    Mask_signal = uint8(roipoly(img,[li(1,:),flip(ma(1,:))],[li(2,:),flip(ma(2,:))]));
    Mask_noise = uint8(roipoly(img,[li(1,:),flip(ma(1,:))],[li(2,:)-thick,flip(ma(2,:))-thick]));

    ms = single(img(logical(Mask_signal)==1));
    mn = single(img(logical(Mask_noise)==1));

    s.ps = sum(ms.^2); s.pn = sum(mn.^2);

    SNR = 20*log10(s.ps/s.pn);

    %% A1 profiles and mask

    li_A1 = li; ma_A1 = ma; mask_A1 = Mask_signal; 
    [Stat_A1]= LI_MA_stats(li,ma,li_A1,ma_A1,CF,filename,'mm');
    
    %% A1s profiles and mask
    li = load(fullfile(Profiles_dir_A1s,[filename '-LI.txt']));
    ma = load(fullfile(Profiles_dir_A1s,[filename '-MA.txt']));

    Mask_signal = uint8(roipoly(img,[li(1,:),flip(ma(1,:))],[li(2,:),flip(ma(2,:))]));

    li_A1s = li; ma_A1s = ma; mask_A1s = Mask_signal; 
    [Stat_A1s]= LI_MA_stats(li_A1,ma_A1,li_A1s,ma_A1s,CF,filename,'mm');
    diceA1s = dice(logical(mask_A1),logical(mask_A1s));

    %% TUM profiles and mask
    li = load(fullfile(Profiles_dir_TUM,[filename '-LI.txt']))';
    ma = load(fullfile(Profiles_dir_TUM,[filename '-MA.txt']))';

    Mask_signal = uint8(roipoly(img,[li(1,:),flip(ma(1,:))],[li(2,:),flip(ma(2,:))]));

    li_TUM = li; ma_TUM = ma; mask_TUM = Mask_signal; 
    [Stat_TUM]= LI_MA_stats(li_A1,ma_A1,li_TUM,ma_TUM,CF,filename,'mm');
    dice_TUM = dice(logical(mask_A1),logical(mask_TUM));


    %% summary

    tempStat.filename = filename;
    tempStat.CF = CF;
    tempStat.SNR = SNR;

    tempStat.DiceA1s = diceA1s;
    tempStat.DiceTUM = dice_TUM;

    tempStat.BiasA1s = Stat_A1s.BiasIMT;
    tempStat.BiasTUM = Stat_TUM.BiasIMT;

    tempStat.HLI_A1s = Stat_A1s.HM_LI;
    tempStat.HLI_TUM = Stat_TUM.HM_LI;

    tempStat.HMA_A1s = Stat_A1s.HM_MA;
    tempStat.HMA_TUM = Stat_TUM.HM_MA;

    Stats(jj) = tempStat;
end


%% evaluate which manual segmentation to get

diceA1s = [Stats(:).DiceA1s];
diceTUM = [Stats(:).DiceTUM];
figure,
subplot(121),histogram(diceA1s);
subplot(122),histogram(diceTUM);

corStat = struct2table(Stats);
corStat = corStat(:,3:end);

names = fieldnames(Stats);
corStat.Properties.VariableNames = names(3:end);

R = corrplot(corStat);

p50_diceA1 = prctile(diceA1s,75);



