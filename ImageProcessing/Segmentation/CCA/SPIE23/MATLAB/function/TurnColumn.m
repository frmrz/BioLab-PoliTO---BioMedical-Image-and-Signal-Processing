function X = TurnColumn(x)
%
% This routine simply turns the input vector x into a row vector X.
% If x is in column format, it is transposed and assigned to X.
% If x is in row format, nothing is done.

S = size(x);
if S(1) > S(2) && S(1) > 1
    %x is in row format
    X = x;
    flag = 0;
elseif S(1) == 1 || S(1) <= S(2)
    X = x';
    flag = 1;
end