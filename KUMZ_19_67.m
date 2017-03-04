function [S] = KUMZ_19_67(S_in,n)
% функция преобразразования 19 центрового расписания в 31 центровой -
% добавление тележек

c_j = zeros(1,n);
S=[];

for i=1:numel(S_in)
   cur_job = S_in(1,i); 
   c_j(cur_job) = c_j(cur_job) + 1;
   switch c_j(cur_job)
       case {1,3,7,9,13,15}
           S = [S cur_job cur_job cur_job cur_job cur_job cur_job cur_job cur_job cur_job];
       otherwise    
           S = [S cur_job];
   end        
end    
end

