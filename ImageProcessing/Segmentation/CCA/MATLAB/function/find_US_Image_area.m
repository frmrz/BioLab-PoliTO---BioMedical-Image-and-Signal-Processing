
function [Bounds,Img_cropped,row,col]=find_US_Image_area(Img)

        [row,col] = size(Img);
        Ior = Img;

%         figure,imshow(Img)
        Img = imopen(Img,strel('disk',3));
        
        %% 1) find image area by COL entropy
        
        [first_col,last_col] = find_max_connected_entro_cols(Img,0);      
        
        %% zero detected column
        Img=Ior;
        Img(:,1:first_col)=0;
        Img(:,last_col:end)=0;
%         figure,imshow(Img);
        
        %% 2) find image area by ROW entropy

        Img = imopen(Img,strel('disk',2));
        
        %% 2.1) suppress boundaries
        entro_row = zeros(1,row);

        for r = 1 : size(Img,1)
            entro_row(r) = entropy(Img(r,:));
        end
        
        % suppress rows before all_zero_rows in first 0.15 of the image   

        row_relative_magnitude = sum( Img(1:round(row*0.15),:)<=25,2);
        row_margin_entro = entro_row(1,1:round(row*0.15));
        
        all_zero_row = row_relative_magnitude>col*0.90 & (row_margin_entro<=4)';
        
        last_all_zero = find(all_zero_row==1,1,'last');
        Img(1:last_all_zero,:)=0;

        clear entro_row
        
        % measure row entropy

        entro_row = zeros(1,row);
                
        for r = 1 : size(Img,1)
            entro_row(r) = entropy(Img(r,:));
        end
        
        % suppress rows after all_zero_rows in last 0.2 of the image

        row_relative_magnitude = sum( Img(round(row*0.8):end,:)<=25,2);
        row_margin_entro = entro_row(1,round(row*0.8):end);
        
        all_zero_row = row_relative_magnitude>col*0.90 & (row_margin_entro<=4)';
        not_row = ~all_zero_row;
        longest = max(accumarray(nonzeros((cumsum(~not_row)+1).*not_row),1));
        
        if longest > 35 
            idx = find(all_zero_row==0,1,'last');
            all_zero_row(1:idx)=0;
        end
        
        last_all_zero = find(all_zero_row==1,1,'first');
        Img(round(row*0.8)+last_all_zero:end,:)=0;

        %% 2.2) find usable image boundaries with entropy
        % find entropy of rows
        clear entro_row entro_col
        
        entro_row = zeros(1,row);
        
        for r = 1 : size(Img,1)
            entro_row(r) = entropy(Img(r,:));
        end 
        
        % find differnce between consecutive entropies
        diff_row = abs(diff(entro_row));
        
        % set the central 50% to zero diff
        diff_row(1,round(row*0.3):round(row*0.7)) = 0;
        
        % find max in first half ROW
        diff_row_first = diff_row;
        diff_row_first(1,round(row*0.5):end) = 0;            
        [jump,first_row] = max(diff_row_first);
        if jump<0.5 || entro_row(first_row)>4
            first_row = find(entro_row>2,1,'first');
        end

        % find max in second half ROW
        diff_row_second = diff_row;
        diff_row_second(1:round(row*0.5)) = 0;
        [jump,last_row] = max(diff_row_second);
        if (jump<0.5 || entro_row(last_row)>4) && jump<2
            last_row = find(entro_row>2,1,'last');
        end
        
        %% 3) Crop image based on entro
        Bounds = [first_row,last_row,first_col,last_col];
        Img_cropped = Ior(first_row:last_row,first_col:last_col);
        
%         figure,
%         subplot(321),imshow(Ior);
%         subplot(322),imshow(Ior);
%         hold on,rectangle('Position',[Bounds(3),Bounds(1),Bounds(4)-Bounds(3),Bounds(2)-Bounds(1)],...
%                           'EdgeColor','r');
%                       
%         subplot(323),plot(flip(entro_row),1:row);
%         subplot(324),plot(1:col,entro_col);      
%         subplot(325),plot(flip(diff_row),1:row-1);
%         subplot(326),plot(1:col-1,diff_col);    
%         
%         pause(0.001)

end

function [first_col,last_col] = find_max_connected_entro_cols(Img,show)
        
        [~,col] = size(Img);
        
        entro_col = zeros(1,col);
        
        for r = 1 : col
            entro_col(r) = entropy(Img(:,r));
        end
        
        % define threshold for entro col
        T = otsuthresh(entro_col); T = T + 0.05*T;
        
        % find largest connected area for entro > T
        entro_up = entro_col > T*max(entro_col); % sections with entro > T
        diff_entro_up = diff(entro_up); 
        
        rising_edge = find(diff_entro_up==1);
        found = 0;
        % define first and last col
        if ~isempty(rising_edge)
            longest = max(accumarray(nonzeros((cumsum(~entro_up)+1).*entro_up),1));
            for q = 1 : length(rising_edge)
                if sum(diff_entro_up(1,rising_edge(q)+1:rising_edge(q)+longest-1)~=0)==0
                    first_col = rising_edge(q);
                    last_col = first_col+longest;
                    found=1;
                    break
                end
            end
        else if isempty(rising_edge) || max(accumarray(nonzeros((cumsum(~entro_up)+1).*entro_up),1)) < 100
                first_col = 1;
                last_col = col;
            end
        end
        
        % first adjust
        if found == 1
            stop = 0;
            while stop == 0 && first_col > 0.05*col
                new_first = find(diff_entro_up(1,first_col- round(col*0.05):first_col-1)==1,1,'first');
                if ~isempty(new_first) && (first_col - new_first)>1
                    first_col = first_col - new_first;
                else
                    stop =1;
                end
            end
        end
        
        % last adjust
        all_up=0;
        if found == 1
            stop = 0;
            while stop == 0 && last_col < 0.95*col
                new_last = find(diff_entro_up(1,last_col+1:last_col + round(col*0.05))==-1,1,'first');
                if ~isempty(new_last)
                    last_col = last_col + new_last;
                else if (diff_entro_up(1,last_col+1)==1) && sum(diff_entro_up(1,last_col+2:last_col + round(col*0.05)))==0
                        last_col = last_col + round(col*0.05);
                        all_up=1;
                    else if all_up==1 && sum(entro_up(1,last_col+1:last_col + round(col*0.05)))==round(col*0.05) % && sum(diff_entro_up(1,last_col+1:last_col + round(col*0.05)))==0 && 
                            last_col = last_col + round(col*0.05);
                else
                    stop =1; all_up=0;
                        end
                    end
                end
            end
        end
        
        % show results
        if show == 1
            figure,
            subplot(211),imshow(Img),hold on,xline(last_col,'r'),hold on,xline(first_col,'r')
            subplot(212),plot(entro_col),hold on,plot(1:col,ones(1,col)*T*max(entro_col)),hold on
                         xline(first_col,'r'),hold on,xline(last_col,'r')
        end
end