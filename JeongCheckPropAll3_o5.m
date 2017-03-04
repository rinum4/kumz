function [fl] = JeongCheckPropAll3_o5(S, U, i1, j1, l1, p, n, m, ma, d, su, end_time_ex, C1_0, r)
% F2, (1-2-1-2), reentr, sdst, T
% Jeong (2014)
% prop 5
%важно! в данной функции сменилась модель хранения su

%важно: S - старое расписание в укороченном виде
%             i1, j1,l1 претенденты на включение в расписание
%             U - работы вне частичного расписания

fl=0;
if i1==j1 || j1==l1 || i1==l1 % l* не должно равняться j2
    return;
end

%проверим на всякий случай 3-ьи вхождения
ind = find(S==i1);
if numel(ind)==4
    return;
end
ind2 = find(S==j1);
if numel(ind2)==4
    return;
end
ind3 = find(S==l1);
if numel(ind3)==4
    return;
end

S1=[i1 i1 j1 j1];
%[~,C1_ij,C2_ij]= XieTargFun10(p, S1, n, m, ma, d, su, S, end_time_ex); %d(G)
[~,C1_ij,C2_ij]= XieTargFun10r(p, S1, n, m, ma, d, su, S, end_time_ex, r); %d(G)
S2=[j1 j1 i1 i1];
%[~,C1_ji,~]= XieTargFun10(p, S2, n, m, ma, d, su, S, end_time_ex); %d(G)
[~,C1_ji,~]= XieTargFun10r(p, S2, n, m, ma, d, su, S, end_time_ex, r); %d(G)

sjU_max = max(max(su{1}(j1,U)),max(su{2}(j1,U)));%max(max(su{j1}(:,U))); %индексировать по ячейкам можно, но макс не определен
sjU_max = max(sjU_max,max(su{3}(j1,U)));
sjU_max = max(sjU_max,max(su{4}(j1,U)));

S3=[i1 i1 j1 j1 l1 l1];
%[~,C1_ijl,C2_ijl]= XieTargFun10(p, S3, n, m, ma, d, su, S, end_time_ex); %d(G1)
[~,C1_ijl,C2_ijl]= XieTargFun10r(p, S3, n, m, ma, d, su, S, end_time_ex, r); %d(G1)
S4=[j1 j1 i1 i1 l1 l1];
%[~,~,C2_jil]= XieTargFun10(p, S4, n, m, ma, d, su, S, end_time_ex); %d(G1)
[~,~,C2_jil]= XieTargFun10r(p, S4, n, m, ma, d, su, S, end_time_ex, r); %d(G1)

%для i это должен быть 1-й проход
if isempty(ind)
    %для j это должен быть 1-й проход
    if isempty(ind2)
        if C2_ij + sjU_max <= C1_ijl
            fl=1;
            return;
        end
    elseif ~isempty(ind2) && numel(ind2)==2
        if C2_ij + sjU_max <= C1_ijl && C2_ij <= d(j1)
            fl=1;
            return;
        end
    end
elseif ~isempty(ind) && numel(ind)==2
    if isempty(ind2)
        if end_time_ex(2,i1) <= C1_0 && ... %не универсально - нужно время окончания 1-го прохода на 2-1 машине
                C2_ij + sjU_max <= C1_ijl
            fl=1;
            return;
        end
    elseif ~isempty(ind2) && numel(ind2)==2
        if C1_ij <= C1_ji && ...
                C2_ijl <= C2_jil && ...
                C2_ij<=d(j1)
            fl=1;
            return;
        end
    end
end
if C1_ij <= C1_ji && ...
        C2_ijl <= C2_jil && ...
        C2_ij<=C1_ijl
    fl=1;
    return;
end
if C1_ij <= C1_ji && ...
        C2_ij + sjU_max <= C1_ijl
    fl=1;
    return;
end
end
