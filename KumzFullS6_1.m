function [S] = KumzFullS6_1(S1,S2)
% для универсальной функции необходимо преобразовать укороченное расписание
% полное перестановочное расписание
% расписание вида S1(1-2-1-2)
% ma=[(1 3 1)-S1 (2 3 4)-S1 (1 3 1)-S2 (2 3 4)-S2 (1 3 1)-S2 (2 3 4 5)-S2];
%S121_1=[4 5 4 5 3 1 3 1 2 2];%n=5
%S121_2 =[5 4 5 4 3 2 3 2 1 1];%n=5

%S1 - достаточно 1го вхождения, т.к. второе вхождение не наобходимо
[~,nc] = size(S1);
S=[];
%S1_U=[];

%for i=1:nr
ne=1;
    for j=1:nc
        %ind = find(S2_U==S2(i,j), 1);
        %if isempty(ind)
            cur_job = ones(1,6)*S1(1,j);
        %else
        %    cur_job = ones(1,7)*S2(i,j);
        %end
        %S2_U = [S2_U S2(i,j)];
        %S=[S cur_job];
        S(ne:ne+5)=cur_job;
        ne = ne+6;
    end
%end

[~,nc] = size(S2);
%S=[];
%S2_U=[];

%for i=1:nr
    for j=1:nc
        %ind = find(S2_U==S2(i,j), 1);
        %if isempty(ind)
        %    cur_job = ones(1,6)*S2(i,j);
        %else
            cur_job = ones(1,7)*S2(1,j);
        %end
        %S2_U = [S2_U S2(i,j)];
        %S=[S cur_job];
        S(ne:ne+6)=cur_job;
        ne = ne+7;
    end
%end

end
