% write_txt_file - Function that writes a vector with numerical data to a file with a certain filename in a certain directory.
%
% Syntax: write_txt_file(vector, filename, directory)
%
% Inputs:
%   - vector: The input vector containing numerical data.
%   - filename: The name of the file to be created.
%   - directory: The directory where the file will be created.
%
% Outputs:
%   None
%
% Example:
%   vector = [1 2 3; 4 5 6; 7 8 9];
%   filename = 'output';
%   directory = '/path/to/directory';
%   write_txt_file(vector, filename, directory);
%
% Note: This function assumes that the input vector contains numerical data only.

function write_txt_file(vector,filename,directory)

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