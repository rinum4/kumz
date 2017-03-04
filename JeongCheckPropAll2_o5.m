function [fl] = JeongCheckPropAll2_o5(S, U, i1, j1, p, n, m, ma, d, su, end_time_ex, r)
% F2, (1-2-1-2), reentr, sdst, T
% Jeong (2014)
% prop 1
%�����! � ������ ������� ��������� ������ �������� su

%�����: S - ������ ���������� � ������ ����
%             i1, j1 ����������� �� ��������� � ����������
%             U - ������ ��� ���������� ����������

fl=0;
if i1==j1
    return;
end

%�������� �� ������ ������ 3-�� ���������
ind = find(S==i1);
if numel(ind)==4
    return;
end
ind2 = find(S==j1);
if numel(ind2)==4
    return;
end

S1=[i1 i1 j1 j1];
%[T1_i,C1_i,C2_i]= XieTargFun10(p, S1, n, m, ma, d, su, S, end_time_ex); %d(G)
[T1_i,C1_i,C2_i]= XieTargFun10r(p, S1, n, m, ma, d, su, S, end_time_ex, r); %d(G)

sjU_max = max(max(su{1}(j1,U)),max(su{2}(j1,U)));%max(max(su{j1}(:,U))); %������������� �� ������� �����, �� ���� �� ���������
sjU_max = max(sjU_max,max(su{3}(j1,U)));
sjU_max = max(sjU_max,max(su{4}(j1,U)));

S2=[j1 j1 i1 i1];

%[T2_j,C1_j,C2_j]= XieTargFun10(p, S2, n, m, ma, d, su, S, end_time_ex); %d(G)
[T2_j,C1_j,C2_j]= XieTargFun10r(p, S2, n, m, ma, d, su, S, end_time_ex, r); %d(G)
siU_min = min(min(su{1}(i1,U)),min(su{2}(i1,U)));%min(min(su{i1}(:,U))); %������������� �� ������� �����, �� ��� �� ���������
siU_min = min(siU_min,min(su{3}(i1,U)));
siU_min = min(siU_min,min(su{4}(i1,U)));

S3=[j1 j1];
%[~,~,C2]= XieTargFun10(p, S3, n, m, ma, d, su, S, end_time_ex);
[~,~,C2]= XieTargFun10r(p, S3, n, m, ma, d, su, S, end_time_ex, r);

%��� i ��� ������ ���� 1-� ������
if isempty(ind)
    %��� j ��� ������ ���� 2-� ������, �� ������� �� ������)))
    if ~isempty(ind2) && numel(ind2)==2
        if C2_i + sjU_max <= C2_j + siU_min && C2_i <= d(j1)
            fl=1;
            return;
        end
    end
elseif ~isempty(ind) && numel(ind)==2
    if ~isempty(ind2) && numel(ind2)==2
        if C1_i <= C1_j && C2_i + sjU_max <= C2_j + siU_min && C2_i<=d(j1)
            fl=1;
            return;
        end
    end
    nn=numel(unique(U)) - 2; %���� ��� � U1 ������ ���� 2-� �����, ���� ��������� i2 � j2 ��� ����� 2
    if T1_i <= T2_j && ...
            (T2_j - T1_i)>nn*max(C1_i-C1_j,C2_i-C2_j+(sjU_max-siU_min))
        fl=1;
        return;
    end
end
if C1_i<=C1_j && ...
        C2_i + sjU_max <= C2_j + siU_min && ...
        C2_i <= C2
    fl=1;
    return;
end

end
