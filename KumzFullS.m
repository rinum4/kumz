function [S] = KumzFullS(rS,ma)
% ��� ������������� ������� ���������� ������������� ����������� ����������
% ������ ��������������� ����������

[~,nc] = size(rS);
%S=[];
nk = numel(ma);
S=zeros(1,nk*nc);

%for i=1:nr %������ ����������� 1
    for j=1:nc %parfor
        cur_job = ones(1,nk)*rS(1,j);
        
        %S=[S cur_job];
        S(1+(j-1)*nk:nk*j)=cur_job;
    end
%end
end
