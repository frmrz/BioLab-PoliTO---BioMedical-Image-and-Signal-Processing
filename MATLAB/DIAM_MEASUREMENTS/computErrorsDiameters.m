close all
clear

%% load excel with computed diameters and do stats

mainFolder = '/media/francesco/DEV001/PROJECT-GITHUB/Image_Processing/Segmentation/Ultrasound-OpticNerveSheaths/DATA/EXCEL';
DiamTable = readtable(fullfile(mainFolder,'DiameterMeasurements.xlsx'));

% tabWrite = table();
% tabWrite.Properties.VariableNames = {'Subset','OND Error','OND Abs. Error','Pval','ICC',...
%                                      'ONSD Error','ONSD Abs. Error','Pval','ICC',};

%% all machines
maskVSmanual = DiamTable(:,[1,2,3,8,9,10]);
unetVSmanual = DiamTable(:,[1,4,5,8,9,10]);
fmzVSmanual = DiamTable(:,[1,6,7,8,9,10]);

stat(1) = statTable(maskVSmanual);
stat(2) = statTable(unetVSmanual);
stat(3) = statTable(fmzVSmanual);
stat(1).name = 'msk';stat(2).name = 'net';stat(3).name = 'fmz';

tabWrite = formatStat(stat);
% writetable(tabWrite,fullfile(main_directory,'diametersErrors.xlsx'),'Sheet','full')
clear stat tabWrite

%% mach 1
DiamTable1 = DiamTable(1:207,:);

maskVSmanual = DiamTable1(:,[1,2,3,8,9,10]);
unetVSmanual = DiamTable1(:,[1,4,5,8,9,10]);
fmzVSmanual = DiamTable1(:,[1,6,7,8,9,10]);

stat(1) = statTable(maskVSmanual);
stat(2) = statTable(unetVSmanual);
stat(3) = statTable(fmzVSmanual);
stat(1).name = 'msk';stat(2).name = 'net';stat(3).name = 'fmz';

tabWrite = formatStat(stat);
% writetable(tabWrite,fullfile(main_directory,'diametersErrors.xlsx'),'Sheet','1')
clear stat tabWrite

%% mach 2
DiamTable1 = DiamTable(208:321,:);

maskVSmanual = DiamTable1(:,[1,2,3,8,9,10]);
unetVSmanual = DiamTable1(:,[1,4,5,8,9,10]);
fmzVSmanual = DiamTable1(:,[1,6,7,8,9,10]);

stat(1) = statTable(maskVSmanual);
stat(2) = statTable(unetVSmanual);
[stat(3)] = statTable(fmzVSmanual);
stat(1).name = 'msk';stat(2).name = 'net';stat(3).name = 'fmz';

tabWrite = formatStat(stat);
% writetable(tabWrite,fullfile(main_directory,'diametersErrors.xlsx'),'Sheet','2')
clear stat tabWrite

%% mach 1
DiamTable1 = DiamTable(321:414,:);

maskVSmanual = DiamTable1(:,[1,2,3,8,9,10]);
unetVSmanual = DiamTable1(:,[1,4,5,8,9,10]);
fmzVSmanual = DiamTable1(:,[1,6,7,8,9,10]);

stat(1) = statTable(maskVSmanual);
stat(2) = statTable(unetVSmanual);
stat(3) = statTable(fmzVSmanual);
stat(1).name = 'msk';stat(2).name = 'net';stat(3).name = 'fmz';

tabWrite = formatStat(stat);
% writetable(tabWrite,fullfile(main_directory,'diametersErrors.xlsx'),'Sheet','3')
clear stat tabWrite

%% mach 1
DiamTable1 = DiamTable(414:end,:);

maskVSmanual = DiamTable1(:,[1,2,3,8,9,10]);
unetVSmanual = DiamTable1(:,[1,4,5,8,9,10]);
fmzVSmanual = DiamTable1(:,[1,6,7,8,9,10]);

stat(1) = statTable(maskVSmanual);
stat(2) = statTable(unetVSmanual);
[stat(3)] = statTable(fmzVSmanual);
stat(1).name = 'msk';stat(2).name = 'net';stat(3).name = 'fmz';

tabWrite = formatStat(stat);
% writetable(tabWrite,fullfile(main_directory,'diametersErrors.xlsx'),'Sheet','4')
clear stat tabWrite


function [stat,ErrorOND,ErrorONSD] = statTable(tempTable)

    tempTable.Properties.VariableNames = {'ID', 'OND', 'ONSD', 'ONDgt', 'ONSDgt', 'Machine'};
    
    for i = 1 : height(tempTable)
        ErrorOND(i,1) = tempTable(i,2).OND-tempTable(i,4).ONDgt;
        ErrorONSD(i,1) = tempTable(i,3).ONSD-tempTable(i,5).ONSDgt;
        OND(i,1) = tempTable(i,2).OND;
        ONDgt(i,1) = tempTable(i,4).ONDgt;
        ONSD(i,1) = tempTable(i,3).ONSD;
        ONSDgt(i,1) = tempTable(i,5).ONSDgt;
    end
    tempOND = OND; tempONSD = ONSD;
    tempONDgt = ONDgt; tempONSDgt = ONSDgt;

    OND = OND(~isnan(tempOND) & ~isnan(tempONDgt));
    ONDgt = ONDgt(~isnan(tempOND) & ~isnan(tempONDgt));
    ONSD = ONSD(~isnan(tempONSD) & ~isnan(tempONSDgt));
    ONSDgt = ONSDgt(~isnan(tempONSD) & ~isnan(tempONSDgt));

    stat.meanErrorOND = mean(ErrorOND,'omitnan');    
    stat.stdErrorOND = std(ErrorOND,'omitnan');
    stat.meanErrorONSD = mean(ErrorONSD,'omitnan');
    stat.stdErrorONSD = std(ErrorONSD,'omitnan');

    stat.meanAbsErrorOND = mean(abs(ErrorOND),'omitnan');    
    stat.stdAbsErrorOND = std(abs(ErrorOND),'omitnan');
    stat.meanAbsErrorONSD = mean(abs(ErrorONSD),'omitnan');
    stat.stdAbsErrorONSD = std(abs(ErrorONSD),'omitnan');

    stat.meanSqrdErrorOND = mean((ErrorOND).^2,'omitnan');    
    stat.stdSqrdErrorOND = std((ErrorOND).^2,'omitnan');
    stat.meanSqrdErrorONSD = mean((ErrorONSD).^2,'omitnan');
    stat.stdSqrdErrorONSD = std((ErrorONSD).^2,'omitnan');

    stat.PvalOND = signrank(OND,ONDgt);
    stat.PvalONSD = signrank(ONSD,ONSDgt);
    
    stat.iccOND = ICC([OND,ONDgt],'1-1');
    stat.iccONSD = ICC([ONSD,ONSDgt],'1-1');

    stat.notMeasured = height(tempTable)-length(OND);
end

function [writeTab] = formatStat(stat)
    tempTab = struct2table(stat);
    writeTab = struct('name',[],'OND_Error',[],'ONSD_Error',[],'OND_absError',[],'ONSD_absError',[],...
        'PvalOND',[],'PvalONSD',[],'iccOND',[],'iccONSD',[],'notMeasured',[]);

    for i = 1 : height(tempTab)

        temp.name = stat(i).name;
        temp.OND_Error = [sprintf('%.3f',stat(i).meanErrorOND) char(177) sprintf('%.3f',stat(i).stdErrorOND)];
        temp.ONSD_Error = [sprintf('%.3f',stat(i).meanErrorONSD) char(177) sprintf('%.3f',stat(i).stdErrorONSD)];
        temp.OND_absError = [sprintf('%.3f',stat(i).meanAbsErrorOND) char(177) sprintf('%.3f',stat(i).stdAbsErrorOND)];
        temp.ONSD_absError = [sprintf('%.3f',stat(i).meanAbsErrorONSD) char(177) sprintf('%.3f',stat(i).stdAbsErrorONSD)];        temp.PvalOND = stat(i).PvalOND;
        temp.PvalONSD = stat(i).PvalONSD;
        temp.iccOND = stat(i).iccOND;
        temp.iccONSD = stat(i).iccONSD;
        temp.notMeasured = stat(i).notMeasured;

        writeTab = [writeTab;temp];
    end
        writeTab = struct2table(writeTab(2:end));

end
