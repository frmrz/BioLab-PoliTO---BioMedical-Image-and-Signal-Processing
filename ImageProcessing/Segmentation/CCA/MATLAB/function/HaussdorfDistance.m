function [dH] = HaussdorfDistance(A,B)

A = TurnColumn(A);
B = TurnColumn(B);

% A with respect to B
dist    = zeros(length(B),1);
mindist = zeros(length(A),1);

for i=1:length(A)    
    pt = A(i,:);
    for j=1:length(B)
        dist(j) = sqrt((pt(1,1)-B(j,1))^2 + (pt(1,2)-B(j,2))^2);
    end
    mindist(i)=min(dist);
end

maxdistA = max(mindist);

% B with respect to A
for i=1:length(B)
    pt = B(i,:);
    for j=1:length(A)
          dist(j) = sqrt((pt(1,1)-A(j,1))^2 + (pt(1,2)-A(j,2))^2);
    end
    mindist(i)=min(dist);
end

maxdistB = max(mindist);

dH = max(maxdistA,maxdistB);