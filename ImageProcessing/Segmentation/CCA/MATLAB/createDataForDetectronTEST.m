_% this script takes the manual segmentation and creates mask for the cnn training

clear  
close all

%% define images and GT directories


Data_fd = cd;

dir_databases = dir(Data_fd);
cont = 0;

for j = 3:size(dir_databases,1)

    if contains(dir_databases(j).name,'DATABASE')

        Img_fd = fullfile(Data_fd,dir_databases(j).name,'IMAGES');
        CF_fd  = fullfile(Data_fd,dir_databases(j).name,'CF');
        Profiles_dir = fullfile(Data_fd,dir_databases(j).name,'GT1-LIMA-Profiles');
        Profiles_dir2 = fullfile(Data_fd,dir_databases(j).name,'GT1bis-LIMA-Profiles');
        Profiles_dir3 = fullfile(Data_fd,dir_databases(j).name,'GT2-LIMA-Profiles');
        
%         Mask_fd  = create_dir(fullfile(Data_fd,dir_databases(j).name,'MASKS'));
%         Save_resize_mask = create_dir(fullfile(Data_fd,dir_databases(j).name,'MASKS-RESIZED'));
%         Save_resize_image = create_dir(fullfile(Data_fd,dir_databases(j).name,'IMAGES-RESIZED'));
%         debug_dir = create_dir(fullfile(Data_fd,dir_databases(j).name,'GT-debug'));

        Img_dir = dir(Img_fd);
%         Mask_dir = dir(Mask_fd);
        CF_dir = dir(CF_fd);
        
        myMap = [ 0.75 0 0; 1 1 0; 0 0.75 0];
        
        %% loop to extract and store masks
        
        ws=[];hs=[];
        
        for i = 3 : length(Img_dir)
            if ~Img_dir(i).isdir && ~strcmp(Img_dir(i).name,'.DS_Store') ...
                    && ~contains(Img_dir(i).name,'Thumbs') % 
            %if contains(Img_dir(i).name,'AAAPE.')%
                
                cont = cont + 1;

                if isfile(fullfile(CF_fd,[Img_dir(i).name(1:end-4) '_CF.txt']))
                    CF = load(fullfile(CF_fd,[Img_dir(i).name(1:end-4) '_CF.txt']));                                 
                else
                    CF = 0.06;
                end
        
                %% Image
                if ~contains(dir_databases(j).name,'Japan')
                    filename = strtok(Img_dir(i).name,'.');
                else
                    filename = Img_dir(i).name(1:end-4);
                end

                images{cont} = filename;

                fprintf('%s\n',filename);

                if contains(Img_dir(i).name(end-4:end),'dcm')
                    Img = dicomread(fullfile(Img_fd,Img_dir(i).name));
                else
                    Img = imread(fullfile(Img_fd,Img_dir(i).name));
                end

                if size(Img,3) == 3
                    Img = rgb2gray(Img);
                end

                [name,number] = strtok(filename,' ');
                if ~isempty(number)
                    filename = [name number(2:end)];
                end

                % Create mask
                if isfile(fullfile(Profiles_dir,[filename '-GTLI.txt']))
                    LI = load(fullfile(Profiles_dir,[filename '-GTLI.txt']));
                    MA = load(fullfile(Profiles_dir,[filename '-GTMA.txt']));
                elseif isfile(fullfile(Profiles_dir2,[filename '-GTLI.txt']))
                    LI = load(fullfile(Profiles_dir2,[filename '-GTLI.txt']));
                    MA = load(fullfile(Profiles_dir2,[filename '-GTMA.txt']));
                elseif isfile(fullfile(Profiles_dir3,[filename '-GTLI.txt']))
                    LI = load(fullfile(Profiles_dir3,[filename '-GTLI.txt']));
                    MA = load(fullfile(Profiles_dir3,[filename '-GTMA.txt']));
                end

                %Sort LI and MA based on x coordinates
                [~,indLI] = sort(LI(:,1));
                [~,indMA] = sort(MA(:,1));
                LI = LI(indLI,:);
                MA = MA(indMA,:);

%                 fig = figure('units','normalized','outerposition',[0 0 1 1],'visible','off');
%                 imshow(Img),hold on, plot(LI(:,1),LI(:,2),'r','Linewidth',2);
%                 hold on, plot(MA(:,1),MA(:,2),'g','Linewidth',2);
%                 export_fig ('-jpg',fullfile(debug_dir,[filename '-GTLIMA.jpg']));
%                 close(fig);
        
                Msk = uint8(roipoly(Img,[LI(:,1); flip(MA(:,1))],[LI(:,2); flip(MA(:,2))])).*255;
                
                imwrite(Msk,fullfile(Mask_fd,[filename '.png']));
                
                %Msk = imread(fullfile(Mask_fd,[filename '.png']));
                        
                %% find image area and resize
                try
                    [Bounds,Img_cropped,~,~]=find_US_Image_area(Img);      
                catch
                    [Img_cropped,~,minX,minY] = fm_autocrop(Img,[]);
                    Bounds(1) = minY;
                    Bounds(2) = minY + size(Img_cropped,1)-1;
                    Bounds(3) = minX;
                    Bounds(4) = minX + size(Img_cropped,2)-1;
                end

                Msk_cropped = Msk(Bounds(1):Bounds(2),Bounds(3):Bounds(4));
                if sum(sum(Msk_cropped)) == 0
                    [Img_cropped,~,minX,minY] = fm_autocrop(Img,[]);
                    Bounds(1) = minY;
                    Bounds(2) = minY + size(Img_cropped,1)-1;
                    Bounds(3) = minX;
                    Bounds(4) = minX + size(Img_cropped,2)-1;
                    Msk_cropped = Msk(Bounds(1):Bounds(2),Bounds(3):Bounds(4));
                end

                Bounds_all(cont,:) = Bounds;
        
        %         figure,subplot(221),imshow(Msk_cropped),subplot(222),imshow(Msk),
        %                subplot(223),imshow(Img_cropped),subplot(224),imshow(Img)
                
                [row(cont),col(cont)] = size(Img_cropped);
                rr = CF/0.0747;
                rr_vec(cont) = rr;
                
                Img_cropped = imresize(Img_cropped,[row(cont)*rr,col(cont)*rr]);
                Msk_cropped = imresize(Msk_cropped,[row(cont)*rr,col(cont)*rr],'nearest');
        
                [row_r(cont),col_r(cont)] = size(Img_cropped);
                
                [Img_padded, n(cont), dims(cont,:)] = padding_rectangular(Img_cropped);
                [Msk_padded, ~, ~] = padding_rectangular(Msk_cropped);
                
                % resizing
                Img_final = imresize(Img_padded,[480 480]);
                Msk_final = imresize(Msk_padded,[480 480],'nearest');
                
                imwrite(Msk_final,fullfile(Save_resize_mask,[filename,'.png']))
                imwrite(Img_final,fullfile(Save_resize_image,[filename,'.png']))
                
                clear Bounds Img Img_padded Msk_padded 
                close all
                
            end
        end
    end
end

save('TestSet','images','row','col','row_r','col_r','n','dims','rr_vec','Bounds_all');