function [S] = KumzFullS12(S_1212,ma)
% для универсальной функции необходимо преобразовать укороченное расписание
% полное перестановочное расписание

%расписание вида 1-2-1-2

[~,nc] = size(S_1212);
nk=numel(ma);
%S=[];
S=zeros(1,nc/2*nk);
%S2_U=[];
S2_U = zeros(1,nc);
%for i=1:nr
ne=1;
    for j=1:nc %parfor
        ind = find(S2_U==S_1212(1,j), 1);
        if isempty(ind)
            cur_job = ones(1,6)*S_1212(1,j);
            S(ne:ne+5)=cur_job;
            ne = ne+6;
        else
            cur_job = ones(1,13)*S_1212(1,j);
            S(ne:ne+12)=cur_job;
            ne = ne+13;
        end
        %S=[S cur_job];
        %S2_U = [S2_U S_1212(1,j)];
        S2_U(j) = S_1212(1,j);
    end
    
%end
end
