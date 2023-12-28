
% TurnRow - Converts the input vector x into a row vector X.
%
% Syntax: X = TurnRow(x)
%
% Inputs:
%   x - Input vector (column or row format)
%
% Outputs:
%   X - Row vector
%
% Example:
%   x = [1; 2; 3]; % column format
%   X = TurnRow(x); % X = [1 2 3]
%
%   x = [1 2 3]; % row format
%   X = TurnRow(x); % X = [1 2 3]
%

function X = TurnRow(x)


S = size(x);
if S(1) < S(2) && S(1) > 1
    %x is in row format
    X = x;
    flag = 0;
elseif S(1) == 1 || S(1) >= S(2)
    X = x';
    flag = 1;
end