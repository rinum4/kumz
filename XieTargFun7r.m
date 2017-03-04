function [sumT, Cmax1, Cmax2, end_time_ex] = XieTargFun7r(p, S, n, m, ma, d, su, r)
% универсальная функция для Fm | reentr | Sum{Tmax}
%важно! в данной функции сменилась модель хранения su
%важно! в данной функции начинаем использовать r- ограничение начало не
%ранее

%важно: здесь S - полное расписание M1-M2-M3-...M1-M2-M3-...

cur_job_state=zeros(n,1);
last_job_ws = zeros(m,1); %3 РЦ
ma_ex = 1:numel(ma);

%для оптимизации памяти сразу сгенерим массив
end_time = zeros(max(ma),n);
end_time_ex = zeros(numel(ma),n);

[nr,nc] = size(S);
prev_job = zeros(max(ma),1);
for i=1:nr
    for j=1:nc
        cur_job = S(i,j);
        cur_job_state(cur_job) =  cur_job_state(cur_job) + 1;
        
        cur_ws = ma(cur_job_state(cur_job));
        cur_ws_ex = ma_ex(cur_job_state(cur_job));
        
        if cur_ws==0
            cur_ws = m;
        end
                
        if  cur_job_state(cur_job) ==1
            %r важно только для 1 работы заказа
            if prev_job(cur_ws) == 0
                %1-й заказ
                %st_t = last_job_ws(cur_ws);  
                st_t = r(cur_job); %max(st_t, r(cur_job));
                end_time(cur_ws, cur_job) = st_t + p(cur_job_state(cur_job), cur_job); %last_job_ws(cur_ws)
            else
                st_t = last_job_ws(cur_ws);  
                st_t = max(st_t, r(cur_job));
                end_time(cur_ws, cur_job) = st_t + p(cur_job_state(cur_job), cur_job) + ...
                                              su{cur_ws_ex}(cur_job, prev_job(cur_ws)); % last_job_ws(cur_ws) su{cur_ws_ex}(prev_job(cur_ws), cur_job) su{prev_job(cur_ws)}(cur_ws_ex, cur_job)
            end    
            %last_job_ws(cur_ws)= end_time(cur_ws, cur_job);
        else
            prev_ws = ma(cur_job_state(cur_job)-1);
            if prev_job(cur_ws) == 0
                end_time(cur_ws, cur_job) = max(end_time(prev_ws, cur_job),last_job_ws(cur_ws)) + p(cur_job_state(cur_job), cur_job);
            else    
                end_time(cur_ws, cur_job) = max(end_time(prev_ws, cur_job)- su{cur_ws_ex}(cur_job, prev_job(cur_ws)),last_job_ws(cur_ws))+...
                                               p(cur_job_state(cur_job), cur_job) + ...
                                               su{cur_ws_ex}(cur_job, prev_job(cur_ws)); % su{cur_ws_ex}(prev_job(cur_ws), cur_job) su{prev_job(cur_ws)}(cur_ws_ex, cur_job)
            end                               
            %last_job_ws(cur_ws)= end_time(cur_ws, cur_job);
        end
        last_job_ws(cur_ws)= end_time(cur_ws, cur_job);
        end_time_ex(cur_ws_ex, cur_job) = end_time(cur_ws, cur_job);
        prev_job(cur_ws) = cur_job;
    end
end
%снизим размерность
[i,j] = find(end_time_ex);
end_time_ex = end_time_ex(:,1:max(j));
end_time_ex = end_time_ex(1:max(i),:);

[i,j] = find(end_time);
end_time = end_time(:,1:max(j));
end_time = end_time(1:max(i),:);

d1=d(1:max(j));

T=max(0,end_time(ma(numel(ma)),:)-d1);
sumT = sum(T);
Cmax1 = end_time(prev_ws, cur_job);
Cmax2 = end_time(cur_ws, cur_job);
end

