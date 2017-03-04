function [S] = KumzFullS6(S1,S2,ma)
% для универсальной функции необходимо преобразовать укороченное расписание
% полное перестановочное расписание
% расписание вида S1(1-2-1-2)
% ma=[(1 3 1)-S1 (2 3 4)-S1 (1 3 1)-S2 (2 3 4)-S2 (1 3 1)-S2 (2 3 4 5)-S2];
%S121_1=[4 5 4 5 3 1 3 1 2 2];%n=5
%S121_2 =[5 4 5 4 3 2 3 2 1 1];%n=5

%S1 - достаточно 1го вхождения, т.к. второе вхождение не наобходимо
S1_un=unique(S1,'stable');
[~,nc] = size(S1_un);
%S=[];
nk=numel(ma);
S=zeros(1,nc*nk);
%for i=1:nr
ne=1;
    for j=1:nc
        cur_job = ones(1,6)*S1_un(1,j);
        
        %S=[S cur_job];
        S(ne:ne+5)=cur_job;
        ne = ne+6;
    end
%end

[~,nc] = size(S2);
%S=[];
%S2_U=[];
S2_U = zeros(1,nc);

%for i=1:nr
    for j=1:nc
        ind = find(S2_U==S2(1,j), 1);
        if isempty(ind)
            cur_job = ones(1,6)*S2(1,j);
            S(ne:ne+5)=cur_job;
            ne = ne+6;
        else
            cur_job = ones(1,7)*S2(1,j);
            S(ne:ne+6)=cur_job;
            ne = ne+7;
        end
        
        %S=[S cur_job];
        %S2_U = [S2_U S2(i,j)];
        S2_U(j) = S2(1,j);
    end
%end

end
