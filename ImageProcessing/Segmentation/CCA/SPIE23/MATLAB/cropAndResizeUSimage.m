% this script takes the manual segmentation and creates mask for the cnn training

clear  
close all

%% define images and GT directories
Data_fd = '/media/francesco/DEV001/PROJECT-CUBS-DSP/DATA480/TEST/';

dir_databases = dir(Data_fd);
cont = 0;

Img_fd = fullfile(Data_fd,'IMAGES-common');
CF_fd  = fullfile(Data_fd,'CF');
       
Mask_fd  = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/MASKS-GT1';
Save_resize_mask = '/media/francesco/DEV001/PROJECT-CUBS-DSP/DATA480/TEST/MASKS/GT1';
Profile_dir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT1-INTERPOLATED';

Save_debug = fullfile(Data_fd,'DEBUG');
 
Img_dir = dir(Img_fd);
Mask_dir = dir(Mask_fd);
CF_dir = dir(CF_fd);

myMap = [ 0.75 0 0; 1 1 0; 0 0.75 0];

%% loop to extract and store masks

for i = 3 : length(Img_dir)
    if ~Img_dir(i).isdir && ~strcmp(Img_dir(i).name,'.DS_Store') ...
            && ~contains(Img_dir(i).name,'Thumbs')

        cont = cont + 1;

        if isfile(fullfile(CF_fd,[Img_dir(i).name(1:end-5) '_CF.txt']))
            CF = load(fullfile(CF_fd,[Img_dir(i).name(1:end-5) '_CF.txt']));                                 
        else
            CF = 0.06;
        end

        %% Image
        filename = Img_dir(i).name(1:end-4); % PNG if TEST
%         filename = Img_dir(i).name(1:end-5); % TIFF if TRAIN

        Img = imread(fullfile(Img_dir(i).folder,Img_dir(i).name));

        %% Mask
        try
            LI = load(fullfile(Profile_dir, [filename '-LI.txt']));
            MA = load(fullfile(Profile_dir, [filename '-MA.txt']));
            Msk = uint8(roipoly(Img,[LI(1,:),flip(MA(1,:))],...
                [LI(2,:),flip(MA(2,:))])).*255;
    
            fprintf('%s\n',filename);
    
            if size(Img,3) == 3
                Img = rgb2gray(Img);
            end
    
            %% find image area and resize parameters
            [row,col] = size(Img);
    
            try
                [Bounds,Img_cropped,~,~]=find_US_Image_area(Img);    % first crop based on info
                [row_c,col_c] = size(Img_cropped);
            catch
                fprintf(' Error in resizing image %s, cropping with alternative function\n',filename);
                [Img_cropped,~,minX,minY] = fm_autocrop(Img,[]);
                Bounds(1) = minY;
                Bounds(2) = minY + size(Img_cropped,1)-1;
                Bounds(3) = minX;
                Bounds(4) = minX + size(Img_cropped,2)-1;
                [row_c,col_c] = size(Img_cropped);
            end     
    
            Msk_cropped = Msk(Bounds(1):Bounds(2),Bounds(3):Bounds(4));
            
            [row(cont),col(cont)] = size(Img_cropped);
    
            rr = CF/0.0747;
            
            Img_cropped = imresize(Img_cropped,[row(cont)*rr,col(cont)*rr]);
            Msk_cropped = imresize(Msk_cropped,[row(cont)*rr,col(cont)*rr],'nearest');
    
            [row_r,col_r] = size(Img_cropped);
    
            % find ratio of i-th image
            ratio=col_r/row_r;
            to_mean = 1 - ratio;
            
            % padding
            if to_mean >= 0
                col_pad = round(row_r*1);
                col_margin = round((col_pad-col_r)/2);
                
                Img_padded = padarray(Img_cropped,[0 col_margin],0,'both');
                Msk_padded = padarray(Msk_cropped,[0 col_margin],0,'both');
            else
                row_pad = round(col_r/1);
                row_margin = round((row_pad-row_r)/2);
                
                Img_padded = padarray(Img_cropped,[row_margin 0],0,'both');
                Msk_padded = padarray(Msk_cropped,[row_margin 0],0,'both');
            end
            
            [Img_padded1,~,~] = padding_rectangular(Img_cropped);
            [Msk_padded1, ~, ~] = padding_rectangular(Msk_cropped);
    
            if numel(Img_padded1) ~= numel(Img_padded1)
                aaaaa=0;
            end
    
            % resizing
            Img_final = imresize(Img_padded,[480 480]);
            Msk_final = imresize(Msk_padded,[480 480],'nearest');
            
            imwrite(Msk_final,fullfile(Save_resize_mask,[filename,'.png']))
            imwrite(Msk,fullfile(Mask_fd,[filename,'.png']))

            clear Bounds Img Img_padded Msk_padded col_pad col_margin row_pad row_margin row_r col_r ratio LI MA
            close all
        catch
            fprintf('No profile in image %s',filename);
        end
        
    end
end

% save('TestSet','images','row','col','row_r','col_r','n','dims','rr_vec','Bounds_all');