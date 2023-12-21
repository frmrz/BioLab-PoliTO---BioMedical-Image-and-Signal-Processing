function [structArray] = Struct_Empty_To_Nan(A)

    if size(A,1)>size(A,2); dir=1; else dir=2; end
    
    if dir == 1
        B = struct2cell(A);
    else
        B = struct2table(A);
        B = table2cell(B);
    end
    empties = cellfun('isempty',B);
    B(empties) = {NaN};
    fields = fieldnames(A);
    structArray = cell2struct(B, fields, dir);
    
end