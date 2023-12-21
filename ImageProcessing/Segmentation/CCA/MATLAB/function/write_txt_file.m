function write_txt_file(vector,filename,directory)

%Function that writes a vector with NUMERICAL DATA to a file with a certain
%filename in a certain directory


[r c]=size(vector);

if ~isdir(directory)
    mkdir(directory)
end

cd(directory)
if strcmp(filename(end-3:end),'.txt')
    fp=fopen(filename,'w'); 
else
    fp=fopen([filename '.txt'],'w');
end

if ~isempty(vector)
    for i=1:r
       for j=1:c
           fprintf(fp,'%f ',vector(i,j));
       end
       fprintf(fp,'\n');
    end
end

fclose(fp);