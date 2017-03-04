function [S] = JeongFullSp(rS)
% для универсальной функции необходимо преобразовать укороченное расписание
%
%n=numel(ma);
[~,nc] = size(rS);
S=[];

%for i=1:nr %nr =1 всегда
    parfor j=1:nc
        cur_job = rS(1,j);
        tmp = [cur_job cur_job];
        S=[S tmp];
    end
%end
end

