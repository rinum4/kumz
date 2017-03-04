function [S] = KumzFullS11(S1_s,x,S2_e,S2_U)
% для универсальной функции необходимо преобразовать укороченное расписание
% полное перестановочное расписание
% расписание вида S1(1-2-1-2)
%ma=[(1 3 1)-S1 (2 3 4)-S1 (1 3 1)-x (2 3 4)-x (1 3 1)-S2 (2 3 4 5)-S2];
%S121_1=[4 5 4 5 3 1 3 1 2 2];%n=5
%S121_2 =[5 4 5 4 3 2 3 2 1 1];%n=5

[~,nc] = size(x);
S=[];
%упрощаем - x всегда вектор
%for i=1:nr
    for j=nc:-1:1 %parfor
        ind = find(S2_U==x(1,j), 1);
        if isempty(ind)
            ind2 = find(S==x(1,j), 1);
            if isempty(ind2)
                cur_job = ones(1,7)*x(1,j);
            else
                cur_job = ones(1,6)*x(1,j);
            end
        else    
            cur_job = ones(1,6)*x(1,j);
        end    
        S=[cur_job S];
    end
%end

S = [S1_s S S2_e];

end
