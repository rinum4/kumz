function [sumT, Cmax1, Cmax2, end_time_ex] = XieTargFun10r(p, S, n, m, ma, d, su, S_in,end_time_ex_in, r)
% универсальная функция для Fm | reentr | Sum{Tmax}
%важно! в данной функции сменилась модель хранения su

%важно: здесь S - полное расписание M1-M2-M3-...M1-M2-M3-...
%важно: здесь S_in - полное расписание M1-M2-M3-...M1-M2-M3-...

ma_ex = 1:numel(ma);
cur_job_state=zeros(n,1);
end_time = zeros(max(ma),n);
end_time_ex = zeros(numel(ma),n);
last_job_ws = zeros(m,1); %3 РЦ
prev_job = zeros(max(ma),1);

%для оптимизации памяти сразу сгенерим массив
%не универсально - только для М1-М2-М1-М2 и для perm расписания
[nrs,ncs]=size(S_in);
if nrs>0 && ncs>0 
    lust_s = S_in(nrs,ncs);
    
    for i=1:n
        if end_time_ex_in(3,i)>0
            end_time(:,i) = end_time_ex_in(3:4,i);
            if i== lust_s
                last_job_ws=end_time_ex_in(3:4,i);
            end    
        else
            end_time(:,i) = end_time_ex_in(1:2,i);
            if i== lust_s
                last_job_ws=end_time_ex_in(1:2,i);
            end    
        end
    end
    
    prev_job = ones(max(ma),1)*S_in(nrs,ncs);
    
    for i=1:nrs
        for j=1:ncs
            cur_job = S_in(i,j);
            cur_job_state(cur_job) =  cur_job_state(cur_job) + 1;
        end
    end
else
    
end

[nr,nc] = size(S);

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
            last_job_ws(cur_ws)= end_time(cur_ws, cur_job);
        else
            prev_ws = ma(cur_job_state(cur_job)-1);
            if prev_job(cur_ws) == 0
                end_time(cur_ws, cur_job) = max(end_time(prev_ws, cur_job),last_job_ws(cur_ws)) + p(cur_job_state(cur_job), cur_job);
            else    
                end_time(cur_ws, cur_job) = max(end_time(prev_ws, cur_job)- su{cur_ws_ex}(prev_job(cur_ws), cur_job),last_job_ws(cur_ws))+...
                                               p(cur_job_state(cur_job), cur_job) + ...
                                               su{cur_ws_ex}(prev_job(cur_ws), cur_job); %su{prev_job(cur_ws)}(cur_ws_ex, cur_job)
            end                               
            
            last_job_ws(cur_ws)= end_time(cur_ws, cur_job);
        end
        end_time_ex(cur_ws_ex, cur_job) = end_time(cur_ws, cur_job);
        prev_job(cur_ws) = cur_job;
    end
end
%снизим размерность
%[i,j] = find(end_time_ex);
%вариант когда возвращает все
if numel(end_time_ex_in)>0
    end_time_ex = max(end_time_ex, end_time_ex_in);
    end_time = max(end_time, end_time_ex_in(3:4,:));
end
%end_time_ex = end_time_ex(:,1:max(j)); %для ex не снижаем размерность
%end_time_ex = end_time_ex(1:max(i),:);

%[i,j] = find(end_time);
%end_time = end_time(:,1:max(j)); %для timw не снижаем размерность
%end_time = end_time(1:max(i),:);

%d1=d(1:max(j)); %вариант со сниж размерностью
d1=d; %вариант без

T=max(0,end_time(ma(numel(ma)),:)-d1);

sumT = sum(T);

Cmax1 = end_time(prev_ws, cur_job);
Cmax2 = end_time(cur_ws, cur_job);
end

