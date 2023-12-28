
% TurnColumn - Converts a vector into a column vector
%
%   X = TurnColumn(x) converts the input vector x into a column vector X.
%   If x is already in column format, nothing is done. If x is in row
%   format, it is transposed and assigned to X.
%
%   Inputs:
%       x - Input vector
%
%   Outputs:
%       X - Column vector
%
%   Example:
%       x = [1; 2; 3];
%       X = TurnColumn(x);

function X = TurnColumn(x)

S = size(x);
if S(1) > S(2) && S(1) > 1
    %x is in row format
    X = x;
    flag = 0;
elseif S(1) == 1 || S(1) <= S(2)
    X = x';
    flag = 1;
end