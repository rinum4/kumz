function [ LB ] = JeongCheckLB2(S, U, p, d, su, C1_S, C2_S)
% F2, (1-2-1-2), reentr, sdst, T
% Jeong (2014)
% LB - не универсальная!
%важно! в данной функции сменилась модель хранения su

[nrs,ncs] = size(S);
LB=0;

if isempty(U)
     return;
end    

U2 = JeongFindU2(S, U); % set of second-pass sub-jobs included in U 
V = [S(nrs,ncs) U]; %set of sub-jobs that includes all sub-jobs in U and the last sub-job in partial sequence S

Un=unique(U); %тест

p1=p; %p(:,Un);
p_min = inf;
ind_min = [0 0];
pUm1=[];
pUm2=[];

for i=1:numel(Un)
    ind = find(S==Un(i));
    
    if isempty(ind)
        % в S работы нет
        if p1(2,Un(i))<p_min
            p_min = p1(2,Un(i));
            ind_min = [Un(i) 1]; 
        end    
        pUm1(i) = p1(1,Un(i));
        pUm2(i) = p1(2,Un(i));
    else    
        % в S работа уже есть 
        if p1(4,Un(i))<p_min
            p_min = p1(4,Un(i));
            ind_min = [Un(i) 2]; 
        end    
        pUm1(i) = p1(3,Un(i));
        pUm2(i) = p1(4,Un(i));
    end
end
W = U; %set of sub-jobs that includes all sub-jobs in U except a sub-job with minimum processing time on machine 2 in U 
%if ind_min(2)==2
[w,I] = find(W==ind_min(1),1);
W(I) = [];
%end    
A = JeongFindA(S, U);  %set of second-pass sub-jobs of jobs whose first-pass subjobs are in S

pr1 = sort(pUm1);

pr2 = sort(pUm2);

%Un = unique(U);
o = [];
o2=[];
for i=1:numel(Un); 
    i_s = Un(i); % i*
    ind = find(A==i_s);
    %
    if ~isempty(ind)
        % для i* 1-й заход уже есть в S
        o(1,i_s) = p(3,i_s); %1-2-(1)
        sVi= JeongFindSUi2(su, V, i_s);
        o(2,i_s) = sVi + p(4,i_s);
        
        o2(1,i_s) = p(3,i_s); %1-2-(1)
        sWi= JeongFindSUi2(su, W, i_s);
        o2(2,i_s) = sWi + p(4,i_s);
    else
        % для i* заходов в S нет
        o(1,i_s) = 2*p(1,i_s); %(1)-2-...
        sVi= JeongFindSUi2(su, V, i_s);
        o(2,i_s) = 2*sVi + 2*p(2,i_s);
        
        o2(1,i_s) = 2*p(1,i_s); %(1)-2-...
        sWi= JeongFindSUi2(su, W, i_s);
        o2(2,i_s) = 2*sWi + 2*p(2,i_s);
    end
end

ou21 = o(1,U2); %M1
ou22 = o(2,U2); %M2
or1=sort(ou21); %M1
or2=sort(ou22); %M2

o2u21 = o2(1,U2); %M1
o2u22 = o2(2,U2); %M2
o2r1=sort(o2u21); %M1
o2r2=sort(o2u22); %M2

LB=0;
dU2 = d(U2);
dr_sort = sort(dU2);

for r = 1:numel(U2)
    Or1 = sum(or1(1:r)); %M1
    Or2 = sum(or2(1:r)); %M2
    if r>1
        O2r2 = sum(o2r2(1:r-1));
    else
        O2r2 = 0;
    end    
    s1 = C1_S + Or1 + pr2(1);
    %мое дополнение
    %if U2(r)==S(nrs,ncs)
    %    s2 = C2_S +Or2 + pr1(1);
    %else
        s2 = C2_S +Or2;
    %end    
    s3 = C1_S +  pr1(1) + pr2(1) + O2r2;
    %C(r) = max([C1_S + Or1 + pr2(1); C2_S +Or2; C1_S +  pr1(1) + pr2(1) + O2r2]);
    C(r) = max([s1 s2 s3]);
    
    LB = LB + max (0, C(r) - dr_sort(r));
end    

end
