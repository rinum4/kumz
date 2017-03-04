function [ T_S_new ] = S_xchange19(ii,S,ri,sr2,p, n, m, ma, d, su, r, mq, mq_u)
T_S_new=ones(1,ri)*inf; %numel(S)
parfor jj=1:ri %parfor бесмысленнен
    S_new_par = S;
    tmp = S_new_par(ii);
    S_new_par(ii) = S_new_par(jj);
    S_new_par(jj) = tmp;
    S_new = [S_new_par sr2];
    %S_new_full = JeongFullS(S_new); %
    %S_new_full = JeongFullSp(S_new); %
    %важно! вместо d - d1
    %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
    % важно! добавил r
    [T_S_new(jj),~,~,~]= KUMZTargFun19f(p, S_new, n, m, ma, d, su, mq, r, mq_u); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
end
end

