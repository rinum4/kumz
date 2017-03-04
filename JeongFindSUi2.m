function [ sUi_min ] = JeongFindSUi2(su, U_in , i_in)
% найдем все вторые вхождения работы в U
%важно! в данной функции сменилась модель хранения su

%  
%k=1;
U = unique(U_in);
sUi_min = inf;
for i=1:numel(U)
    %U1 = U;
    %U1[i]=[];
    
%     if U(i)==i_in %вопрос пропускать ли самого себя? пробуем не пропускать
%         continue;
%     end    
    %sUi = min(min(su{U(i)}(:,i_in))); %U1 %индексировать по ячейкам можно, но макс не определен
    sUi = min(min(su{1}(U(i),i_in)),min(su{2}(U(i),i_in)));%min(min(su{i1}(:,U))); %индексировать по ячейкам можно, но мин не определен
    sUi = min(sUi,min(su{3}(U(i),i_in)));
    sUi = min(sUi,min(su{4}(U(i),i_in)));
    
    if sUi<sUi_min
        sUi_min = sUi;
    end    
    
end    

end

