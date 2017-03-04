function [ A ] = JeongFindA(S_in, U_in)
% найдем все вторые вхождения работы в U
%  
k=1;
S = unique(S_in);
%в этот раз оптимизируем распределение памяти
A=zeros(1,numel(S));
for i=1:numel(S)
    ind = find(U_in==S(i));
    if ~isempty(ind)
        A(k) =S(i);
        k=k+1;
    end
end    
[i] = find(A);
A=A(1:max(i));
end

