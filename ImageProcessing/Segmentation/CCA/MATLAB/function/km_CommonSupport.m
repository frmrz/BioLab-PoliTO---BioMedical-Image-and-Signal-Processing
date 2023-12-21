function [A,B,flag] = km_CommonSupport(C,D)

C=TurnRow(C);
D=TurnRow(D);

[Cx, indC, ~] = unique(C(1,:));
Cy = C(2,indC);

[Dx, indD, ~] = unique(D(1,:));
Dy = D(2,indD);

clear C D
C = [Cx; Cy];
D = [Dx; Dy];

CS = [max([C(1,1) D(1,1)]) min([C(1,end) D(1,end)])];

if CS(1) <= CS(2)
   
  if C(1,1)<CS(1)
      A1y=spline(C(1,:),C(2,:),CS(1));
      A1=[CS(1);A1y];
  else
        A1=[];
  end
   
  if D(1,1)<CS(1)
      B1y=spline(D(1,:),D(2,:),CS(1));
      B1=[CS(1);B1y];
  else
        B1=[];
  end
   
  if C(1,end)>CS(2)
      Aendy=spline(C(1,:),C(2,:),CS(2));
      Aend=[CS(2);Aendy];
  else
        Aend=[];
  end
   
  if D(1,end)>CS(2)
      Bendy=spline(D(1,:),D(2,:),CS(2));
      Bend=[CS(2);Bendy];
  else
        Bend=[];
  end
       
   ind = find(C(1,:) <= CS(2) & C(1,:) >= CS(1));
   A = [A1 C(:,ind) Aend];
   ind = find(D(1,:) <= CS(2) & D(1,:) >= CS(1));
   B = [B1 D(:,ind) Bend];
   flag = 1;   %there is common support
   
   [Ax indAx]=unique(A(1,:));
   Ay=A(2,indAx);
   
   [Bx indBx]=unique(B(1,:));
   By=B(2,indBx);
   
   clear A B
   A=[Ax;Ay];
   B=[Bx;By];
   
   %% added for SPIE
   A = TurnColumn(A);
   B = TurnColumn(B);

else
    disp('sorry ... no common support!');
    A = C;
    B = D;
    flag = 0;   %there is no common support
end