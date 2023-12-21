function st = del_junk(st)
    i = 1;
    
    while i < height(st)
        if sum(strcmp(st(i).name,{'.','..','.DS_Store'}))>0
            if i == 1
                st = st(i+1:end);
            else
                st = [st(1:i-1);st(i+1:end)];
            end
        else 
            i = i + 1;
        end
    end
    
end