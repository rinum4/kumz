function [S] = KumzFullSp(rS,ma)
% для универсальной функции необходимо преобразовать укороченное расписание
% полное перестановочное расписание

[~,nc] = size(rS);
S=[];

%for i=1:nr %nr - всегда 1
    parfor j=1:nc %parfor
        cur_job = ones(1,numel(ma))*rS(1,j);
        
        S=[S cur_job];
    end
%end
end
