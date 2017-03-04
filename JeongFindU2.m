function [ U2 ] = JeongFindU2(S_in, U_in)
% найдем все вторые вхождения работы в U из S
%  
k=1;
U = unique(U_in);
%в этот раз оптимизируем распределение памяти
U2=zeros(1,numel(U));
for i=1:numel(U)
    ind = find(S_in==U(i));
    if ~isempty(ind)
        U2(k) =U(i);
        k=k+1;
    end
end    
[i] = find(U2);
U2=U2(1:max(i));
end

