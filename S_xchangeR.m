function [ T_S_new ] = S_xchangeR(ii,S,ri,sr2,p1, n, m1, ma1, d1, su1, r)
T_S_new=ones(1,ri)*inf;
for jj=1:ri %parfor ������������
    S_new_par = S;
    tmp = S_new_par(ii);
    S_new_par(ii) = S_new_par(jj);
    S_new_par(jj) = tmp;
    S_new = [S_new_par sr2];
    S_new_full = JeongFullS(S_new); %
    %S_new_full = JeongFullSp(S_new); %
    %�����! ������ d - d1
    %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % ��� 1-2-1-2 ���������, ��� ������ ��������� %(p1, S, n, m1, ma1, d, su1);
    % �����! ������� r
    [T_S_new(jj),~,~,~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % ��� 1-2-1-2 ���������, ��� ������ ��������� %(p1, S, n, m1, ma1, d, su1);
end
end

