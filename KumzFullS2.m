function [S] = KumzFullS2(rS,ma)
% ��� ������������� ������� ���������� ������������� ����������� ����������
% ������ ��������������� ����������
%�����! ���������� ��������������� �� ��

[~,nc] = size(rS);
%S=[];
nk = numel(ma);

for i=1:nk
        cur_job = rS;
        
        %S=[S cur_job];
        S(1+(i-1)*nc:nc*i)=cur_job;
end
end
