%T_min = inf;
%S_min = [];
%C_min = inf;
%S_new_par = zeros(numel(S));
clear S_new_par;
T_S_new=ones(1,numel(S))*inf;
%C2_S_new=ones(1,numel(S))*inf;

parfor l=1:numel(S) +1 %parfor %по тестам выходит, что пар цикл мешает оставим эту реализацию... пока
    if l==1
        S_new_par(l,:) = [sr S];
    elseif l== numel(S) +1
        S_new_par(l,:) = [S sr];
    else
        S_new_par(l,:) = [S(1:l-1) sr S(l:numel(S))];
    end
    %здесь хорошо бы проверить условия доминирования, но алгоритм и так быстрый
    S_new = [S_new_par(l,:) S0(ri+1:numel(S0))];
    S_new_full = JeongFullS(S_new);
    %S_new_full = JeongFullSp(S_new);
    %важно! вместо d - d1
    %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
    % важно! добавил r
    %[T_S_new(l),~,C2_S_new(l), ~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
    [T_S_new(l),~,~,~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
end
[T_min, ind_m] = min(T_S_new);
T_min_old = T_min;
%S=S_min;
S = S_new_par(ind_m,:);

%--------------------------------------------------------------------------

T_min = inf;
S_min = [];
C_min = inf;
%T_S_new=ones(1,numel(S))*inf;
%C2_S_new=ones(1,numel(S))*inf;
for l=1:numel(S) +1
    if l==1
        S_new_par = [sr S];
    elseif l== numel(S) +1
        S_new_par = [S sr];
    else
        S_new_par = [S(1:l-1) sr S(l:numel(S))];
    end
    %здесь хорошо бы проверить условия доминирования, но алгоритм и так быстрый
    S_new = [S_new_par S0(ri+1:numel(S0))];
    S_new_full = JeongFullS(S_new);
    %S_new_full = JeongFullSp(S_new);
    %важно! вместо d - d1
    %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
    % важно! добавил r
    [T_S_new,~,C2_S_new, ~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
    if T_S_new<T_min  || (T_S_new==T_min && C2_S_new<C_min)
        T_min = T_S_new;
        S_min = S_new_par;
        C_min = C2_S_new;
    end
end
S=S_min;
T_min_old = T_min;

%--------------------------------------------------------------------------

%T_min = inf;
%S_min = [];
%C_min = inf;
%S_new_par = zeros(numel(S));
%clear S_new_par;
T_S_new=ones(1,numel(S))*inf;
%C2_S_new=ones(1,numel(S))*inf;
for l=1:numel(S) +1 %parfor %по тестам выходит, что пар цикл мешает оставим эту реализацию... пока
    if l==1
        S_new_par = [sr S];
    elseif l== numel(S) + 1
        S_new_par = [S sr];
    else
        S_new_par = [S(1:l-1) sr S(l:numel(S))];
    end
    %здесь хорошо бы проверить условия доминирования, но алгоритм и так быстрый
    S_new = [S_new_par S0(ri+1:numel(S0))];
    S_new_full = JeongFullS(S_new);
    %S_new_full = JeongFullSp(S_new);
    %важно! вместо d - d1
    %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
    % важно! добавил r
    %[T_S_new(l),~,C2_S_new(l), ~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
    [T_S_new(l),~,~,~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
end
[T_min, ind_m] = min(T_S_new);
T_min_old = T_min;
%S=S_min;
if ind_m == 1
    S = [sr S];
elseif ind_m == numel(S) + 1
    S = [S sr];
else
    S = [S(1:l-1) sr S(l:numel(S))];
end
%S = S_new_par(ind_m,:);