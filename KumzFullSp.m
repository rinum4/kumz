function [S] = KumzFullSp(rS,ma)
% ��� ������������� ������� ���������� ������������� ����������� ����������
% ������ ��������������� ����������

[~,nc] = size(rS);
S=[];

%for i=1:nr %nr - ������ 1
    parfor j=1:nc %parfor
        cur_job = ones(1,numel(ma))*rS(1,j);
        
        S=[S cur_job];
    end
%end
end
