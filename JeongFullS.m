function [S] = JeongFullS(rS)
% для универсальной функции необходимо преобразовать укороченное расписание
%
%n=numel(ma);
[~,nc] = size(rS);
%S=[];
S=zeros(1,nc*2);
%for i=1:nr %nr =1 всегда
    for j=1:nc
        %cur_job = rS(1,j);
        %S=[S cur_job cur_job];
        S(2*j-1)=rS(1,j);
        S(2*j)=rS(1,j);
    end
%end
end

