
% Struct_Empty_To_Nan is a function that replaces empty values in a structure array with NaN.
%
% Input:
%   - A: The input structure array
%
% Output:
%   - structArray: The modified structure array with empty values replaced by NaN
%
% Example:
%   A = struct('field1', [], 'field2', [1 2 3], 'field3', []);
%   structArray = Struct_Empty_To_Nan(A);
%
%   The resulting structArray will be:
%   structArray = struct('field1', NaN, 'field2', [1 2 3], 'field3', NaN);
%
% Author: Francesco Marzola
%
function [structArray] = Struct_Empty_To_Nan(A)
    % Check the direction of the structure array
    if size(A,1) > size(A,2)
        dir = 1; % Vertical direction
    else
        dir = 2; % Horizontal direction
    end
    
    % Convert the structure array to cell array
    if dir == 1
        B = struct2cell(A);
    else
        B = struct2table(A);
        B = table2cell(B);
    end
    
    % Find empty cells and replace them with NaN
    empties = cellfun('isempty', B);
    B(empties) = {NaN};
    
    % Get the field names of the input structure array
    fields = fieldnames(A);
    
    % Convert the modified cell array back to structure array
    structArray = cell2struct(B, fields, dir);
end