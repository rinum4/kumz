function [S] = KumzFullS(rS,ma)
% для универсальной функции необходимо преобразовать укороченное расписание
% полное перестановочное расписание

[~,nc] = size(rS);
%S=[];
nk = numel(ma);
S=zeros(1,nk*nc);

%for i=1:nr %всегда размерность 1
    for j=1:nc %parfor
        cur_job = ones(1,nk)*rS(1,j);
        
        %S=[S cur_job];
        S(1+(j-1)*nk:nk*j)=cur_job;
    end
%end
end
