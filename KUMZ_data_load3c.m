function [start_time_full, end_time_full] = KUMZ_data_load3c(p, S, n, m, ma, d, su, mq, r, mq_u, db_data, id_t, id_in)
% универсальная функция для Fm | reentr | Sum{Tmax}
%важно! в данной функции сменилась модель хранения su
%важно! в данной функции начинаем использовать множественность РЦ
%важно! в данной функции начинаем использовать r - начало не ранее
%важно! без начала работ не обойтись
%важно! сохраняем начала работ на рабочем центре
%важно! смена технологии хранения промежутков времени

%важно: здесь S - полное расписание M1-M2-M3-...M1-M2-M3-...

cur_job_state=zeros(n,1); %состояние текущей работы
last_j_full = ones(n,1);
ma_ex = 1:numel(ma);

s_e_ws = cell(m,1);
m_j_s_e_ws = cell(m,1);
for i=1:m %parfor
    s_e_ws{i}=zeros(3,n*18); %18 возможных заходов рц 1
    m_j_s_e_ws{i} = 1;
    if mq_u(i)>1 %для множественных РЦ создадим расшир массив окончаний
        s_e_ws{i}=cell(mq_u(i),1);
        m_j_s_e_ws{i} = ones(mq_u(i),1);
        for j=1:mq_u(i)
            s_e_ws{i}{j}=zeros(3,n*36); %36 возможных заходов рц 6
        end
    end
end

%minq_job_ws = ones(m,1); %если множественный рц, то хранится индекс влож РЦ с мин врем
%minq_job_ws_ex = ones(numel(ma),1); %если множественный рц, то хранится индекс влож РЦ с мин врем

last_job_ws_q = cell(m,1);
for i=1:m %parfor
    if mq_u(i)>1 %для множественных РЦ создадим расшир массив окончаний
        last_job_ws_q{i}=zeros(mq_u(i),n); %окончание последней работы на множ рц (расширенный, с работами)
    end
end

%start_job_ws_q = cell(m,1);
% last_job_ws_q_ex = cell(numel(ma),1);
% for i=1:numel(ma) %parfor
%     if mq(i)>1 %для множественных РЦ создадим расшир массив окончаний
% %       %start_job_ws_q{i}=zeros(mq(i),n);
%         last_job_ws_q_ex{i}=zeros(mq(i),n); %окончание последней работы на множ рц (расширенный, с работами)
%     end
% end

%для оптимизации памяти сразу сгенерим массив
start_time_ex = zeros(numel(ma),n);
start_time_full = zeros(300,n); % время с учетом полной тех карты - %300 с запасом!

end_time = zeros(max(ma),n);
end_time_ex = zeros(numel(ma),n);
end_time_full = zeros(300,n); % время с учетом полной тех карты - %300 с запасом!

[~,nc] = size(S);
prev_job = zeros(max(ma),1);

id_in_u = unique(id_in,'stable');
id_in_u = nonzeros(id_in_u);%find(id_in_u~=0);

%for i=1:nr %лишний цикл nr=1 всегда
%for j=1:nc
j=1;
%last_j = 0;
while j<nc+1 %for меняем на while для управляемости
    cur_job = S(1,j);
    f_srch = id_in_u(cur_job);
    %     if cur_job==24
    %         stop=1;
    %         if cur_job_state(cur_job) == 2
    %             stop = 1;
    %         end
    %     end
    if f_srch==47
        stop=1;
    end
    %cur_job_state(cur_job) =  cur_job_state(cur_job) + 1;
    %для полного расписания нужно понять стадию на которой нах-ся
    %изделие (в этот раз важны и транспортные операции
    if cur_job_state(cur_job) == 0
        id_n = find(id_t==f_srch); %cur_job
        n_id = numel(id_n);
        %stage=1;
        for jj=1:n_id
            if strcmp(db_data{id_n(jj),3},'S01::ST')
                cur_job_state(cur_job) = 1;
                break;
            end % 
            if strcmp(db_data{id_n(jj),3},'CRM::TX')
                cur_job_state(cur_job) = 1;
                break;
            end % 
            if strcmp(db_data{id_n(jj),3},'CRM::RS') 
                cur_job_state(cur_job) = 2;
                break;
            end % 
            if strcmp(db_data{id_n(jj),3},'CTI::ST') || strcmp(db_data{id_n(jj),3},'CTI::TX')
                cur_job_state(cur_job) = 10;
                break;
            end % 
            if strcmp(db_data{id_n(jj),3},'EBC::FX') || strcmp(db_data{id_n(jj),3},'EBC::TX')
                cur_job_state(cur_job) = 20;
                break;
            end % 
            if strcmp(db_data{id_n(jj),3},'G03::CU') || strcmp(db_data{id_n(jj),3},'G03::IN') || strcmp(db_data{id_n(jj),3},'G03::PA') || strcmp(db_data{id_n(jj),3},'G03::PK')
                cur_job_state(cur_job) = 22;
                break;
            end % 
            if strcmp(db_data{id_n(jj),3},'G01::CU') || strcmp(db_data{id_n(jj),3},'G01::IN') || strcmp(db_data{id_n(jj),3},'G01::PA') || strcmp(db_data{id_n(jj),3},'G01::PK') || ...
                    strcmp(db_data{id_n(jj),3},'G02::CU') || strcmp(db_data{id_n(jj),3},'G02::IN') || strcmp(db_data{id_n(jj),3},'G02::PA') || strcmp(db_data{id_n(jj),3},'G02::PK') ||...
                    strcmp(db_data{id_n(jj),3},'FIN01::ST') || strcmp(db_data{id_n(jj),3},'FIN02::ST') || strcmp(db_data{id_n(jj),3},'FIN03::ST')
                cur_job_state(cur_job) = 67;
                break;
            end % 
        end
    else
        id_n = find(id_t==f_srch); %cur_job
        n_id = numel(id_n);
        %if ~strcmp(db_data{id_n(jj),4}{1},'K') && ~strcmp(db_data{id_n(jj),4}{1},'T')
        %    cur_job_state(cur_job) =  cur_job_state(cur_job) + 1;
        %end
        jj=last_j_full(cur_job);
        if jj>n_id
            j=j+1;
            continue;
        end    
        jjj = id_n(jj);
        stage = cur_job_state(cur_job)+1;
        if strcmp(db_data{jjj,3},'CRM::TX') 
            
             if stage ~= 1 && stage ~= 3 && stage ~= 4 && stage ~= 6 && stage ~= 7 && stage ~= 9 && stage ~= 11 && ...
                     stage ~= 13 && stage ~= 14 && stage ~= 16 && stage ~= 17 && stage ~= 19 && stage ~= 23 && stage ~= 25 && ...
                        stage ~= 26 && stage ~= 28 && stage ~= 29 && stage ~= 31 && stage ~= 33 && stage ~= 35 && stage ~= 36 && ...
                      stage ~= 38 && stage ~= 39 && stage ~= 41 && stage ~= 45 && stage ~= 47 && stage ~= 48 && stage ~= 50 && ...  
                    stage ~= 51 && stage ~= 53 && stage ~= 55 && stage ~= 57 && stage ~= 58 && stage ~= 60 && stage ~= 61 && ...    
                      stage ~= 63
                            stage = stage + 1;
             end
            
        elseif strcmp(db_data{jjj,3},'CRM::RS')
            
             if stage ~= 2 && stage ~= 5 && stage ~= 8 && stage ~= 12 && stage ~= 15 && stage ~= 18 && stage ~= 24  && ...
                     stage ~= 27 && stage ~= 30 && stage ~= 34 && stage ~= 37 && stage ~= 40 && stage ~= 46 && stage ~= 49  && ...
                   stage ~= 52 && stage ~= 56 && stage ~= 59 && stage ~= 62 
                            stage = stage + 1;
             end
            
        elseif strcmp(db_data{jjj,3},'CTI::ST')  ||  strcmp(db_data{jjj,3},'CTI::TX')
            
            %забыл зачем это нужно пока закомментирую
            %             if stage == 20 || stage ==42 || stage == 64
            %                 stage = stage - 1;
            %             end
            
            if stage ~= 10 && stage ~= 21 && stage ~= 32 && stage ~= 43 && stage ~= 54 && stage ~= 65
                stage = stage + 1;
            end
            
        elseif strcmp(db_data{jjj,3},'EBC::FX')  ||  strcmp(db_data{jjj,3},'EBC::TX')
            
            if stage ~= 20 && stage ~= 42 && stage ~= 64
                stage = stage + 1;
            end    
            
        elseif  strcmp(db_data{jjj,3},'G03::CU')  ||  strcmp(db_data{jjj,3},'G03::IN')  ||  strcmp(db_data{jjj,3},'G03::PA')  ||  strcmp(db_data{jjj,3},'G03::PK')  %1рц
            
            if stage ~= 22 && stage ~= 44 && stage ~= 66
                stage = stage + 1;
            end
            
        elseif strcmp(db_data{jjj,3},'G01::CU')  ||  strcmp(db_data{jjj,3},'G01::IN')  ||  strcmp(db_data{jjj,3},'G01::PA')  ||  strcmp(db_data{jjj,3},'G01::PK') || ...
                strcmp(db_data{jjj,3},'G02::CU')  ||  strcmp(db_data{jjj,3},'G02::IN')  ||  strcmp(db_data{jjj,3},'G02::PA')  ||  strcmp(db_data{jjj,3},'G02::PK')
            if stage ~= 67
                stage = stage + 1;
            end    
        end
        
        cur_job_state(cur_job) = cur_job_state(cur_job) + 1;
        if cur_job_state(cur_job)>67
            j=j+1;
            continue;
        end
        if stage>cur_job_state(cur_job)
            %если стадия пропускаяется важно записать ее начало и окончание
            cur_ws_ex = ma_ex(cur_job_state(cur_job));
            cur_ws = ma(cur_job_state(cur_job));
            
            if cur_job_state(cur_job)>1
                prev_ws = ma(cur_job_state(cur_job)-1);
                prev_ws_ex = cur_ws_ex - 1;
            else
                prev_ws = 0;
                prev_ws_ex = 0;
            end
            
            if prev_ws_ex > 0
                start_time_ex(cur_ws_ex, cur_job) =  start_time_ex(prev_ws_ex, cur_job);
                end_time_ex(cur_ws_ex, cur_job) = end_time_ex(prev_ws_ex, cur_job);
            end
            if prev_ws > 0
                end_time(cur_ws, cur_job) = end_time(prev_ws, cur_job);
            end    
            j=j+1;
            continue;
        end
    end
    
    cur_ws_ex = ma_ex(cur_job_state(cur_job));
    
    cur_ws = ma(cur_job_state(cur_job));
    %     if cur_ws==2 && p(cur_ws_ex, cur_job)>0 %&& cur_job==8
    %         stop=1;
    %     end
    
    if cur_ws==0
        cur_ws = m;
    end
    
    q_ws = mq(cur_ws_ex); %количество машин в РЦ
    
    if cur_job_state(cur_job)>1
        prev_ws = ma(cur_job_state(cur_job)-1);
        prev_ws_ex = cur_ws_ex - 1;
    else
        prev_ws = 0;
        prev_ws_ex = 0;
    end
    
    %теоретически нужно p(cur_ws_ex, cur_job) переопределить, т.к. p
    %формируется с уcловсностями
    search_sym1 = '';   
    search_sym2 = '';   
    search_sym3 = '';
    search_sym4 = '';
    search_sym5 = '';
    fl_crm_tx = 0;
    switch cur_job_state(cur_job)
        case {10,21,32,43,54,65} %3
            search_sym1 = 'CTI::ST';
            search_sym2 = 'CTI::TX';
        case {1} %6    
            search_sym1 = 'CRM::TX';
            search_sym2 = 'S01::ST';    
        case {3,4,6,7,9,11,13,14,16,17,23,25,26,28,29,31,33,35,36,38,39,45,47,48,50,51,53,55,57,58,60,61} %6
            search_sym1 = 'CRM::TX';
            %search_sym2 = 'S01::ST';    
            fl_crm_tx = 1;
        case {2,5,8,12,15,18,24,27,30,34,37,40,46,49,52,56,59,62} %1
            search_sym1 = 'CRM::RS'; %
        case {19,41,63}% заплатка на случай 1-3-1-3-1
            search_sym1 = 'CRM::RS'; %
            search_sym2 = 'CRM::TX'; %
            search_sym3 = 'CTI::ST';
            search_sym4 = 'CTI::TX';
        case {20,42,64}
            search_sym1 = 'EBC::FX';
            search_sym2 = 'EBC::TX';
        case {22,44,66}
            search_sym1 = 'G03::CU';
            search_sym2 = 'G03::IN';
            search_sym3 = 'G03::PA';
            search_sym4 = 'G03::PK';
            search_sym5 = 'FIN03::ST';
        case 67
            %все оставшиеся
    end
    
    p_curs_ws_cur_job = 0;
    fl_en = 0;
    for jj=last_j_full(cur_job):n_id
        if strcmp(db_data{id_n(jj),3},search_sym1) || strcmp(db_data{id_n(jj),3},search_sym2) || strcmp(db_data{id_n(jj),3},search_sym3) || ...
                strcmp(db_data{id_n(jj),3},search_sym4) || strcmp(db_data{id_n(jj),3},search_sym5) || cur_job_state(cur_job) == 67
            
            p_curs_ws_cur_job = p_curs_ws_cur_job + db_data{id_n(jj),16}; %*db_data{id_n(jj),20} умножать на количество как оказалось не надо!
            fl_en = 1;
            if fl_crm_tx ==1 % тележки CRM:TX сами по себе (не надо склеивать)
                break;
            end    
        else
            if fl_en == 1
                break;
            end    
        end
    end
    
    %вычислим допустимое время последней работы на данном РЦ
    if p(cur_ws_ex, cur_job)>0 || p_curs_ws_cur_job>0
        if prev_ws_ex>0
            s_t_b = end_time_ex(prev_ws_ex, cur_job); %end_time(prev_ws, cur_job);
        else
            s_t_b = r(cur_job);
        end
        
        %s_t_e = s_t_b + p(cur_ws_ex, cur_job)*q_ws; %тут переналадки не катят! + su{cur_ws_ex}(cur_job, prev_job(cur_ws))
        %tmp = p(cur_ws_ex, cur_job)*q_ws*60; %не заывайй про секунды!
        s_t_e = s_t_b + p_curs_ws_cur_job; %тут переналадки не катят! + su{cur_ws_ex}(cur_job, prev_job(cur_ws))
        
        last_job_t = 0;
        %--------------------------------------------------------------------------------------------------------------------------------------
        %важно! данный цикл нужен для множественных РЦ!
        flag=1;
        s_e_ind = 1; % индекс
        t_last_job_t = 0;
        while flag == 1 %
            if q_ws==1
                %[sort_e_t, I_sort] = sort(end_time_ex(cur_ws_ex, :));
                %sort_e_t_prev = start_time_ex(cur_ws_ex, I_sort);
                [sort_e_t, I_sort] = sort(s_e_ws{cur_ws}(2,:));
                sort_e_t_prev = s_e_ws{cur_ws}(1,I_sort);
                sort_e_t_ind = s_e_ws{cur_ws}(3,I_sort);
            else
                %[sort_e_t, I_sort] = sort(last_job_ws_q{cur_ws_ex}(minq_job_ws_ex(cur_ws_ex), :));
                %sort_e_t_prev = start_job_ws_q{cur_ws_ex}(minq_job_ws_ex(cur_ws_ex),  I_sort);
                [sort_e_t, I_sort] = sort(s_e_ws{cur_ws}{s_e_ind}(2,:));
                sort_e_t_prev = s_e_ws{cur_ws}{s_e_ind}(1,I_sort);
                sort_e_t_ind = s_e_ws{cur_ws}{s_e_ind}(3,I_sort);
            end
            
            ind = find(sort_e_t>s_t_e,1);
            nn = numel(sort_e_t);
            if ~isempty(ind)
                fl=0;
                
                for jj = ind-1:nn-1 %ind - 1:n-1
                    %переналадки проверяются тут! т.к. только в этот момент
                    %мы знаем какая работа перед нами
                    %важно! необходимо искать место с 2мя переналадками! с
                    %предыд и послед
                    pr_j = sort_e_t_ind(jj);
                    nx_j = sort_e_t_ind(jj+1);
                    if pr_j > 0
                        t_per = su{cur_ws_ex}(cur_job, pr_j);
                    else
                        t_per = 0;
                    end
                    if nx_j  > 0
                        t_per2 = su{cur_ws_ex}(cur_job, nx_j);
                    else
                        t_per2 = 0;
                    end
                    if (s_t_e - s_t_b + t_per + t_per2) < (sort_e_t_prev(jj+1) - sort_e_t(jj)) && ...
                            s_t_e + t_per < sort_e_t_prev(jj+1) %I_sort(jj)
                        if s_t_b > sort_e_t(jj)
                            last_job_t = s_t_b;
                            %if prev_same_ws_ex > 0
                            %   if sort_e_t(jj)> end_time_ex(prev_same_ws_ex, prev_job(cur_ws)) %cur_job
                            %prev_job(cur_ws) = I_sort(jj);
                            prev_job(cur_ws) =pr_j;
                            %    else
                            %
                            %   end
                            %end
                        else
                            last_job_t = sort_e_t(jj);
                            %prev_job(cur_ws) = I_sort(jj);
                            prev_job(cur_ws) = pr_j;
                        end
                        
                        fl=1;
                        break;
                    end
                end
                if fl==0
                    last_job_t = sort_e_t(nn); %n
                    %prev_job(cur_ws) = I_sort(n);
                    prev_job(cur_ws) = sort_e_t_ind(nn);
                end
            else
                if sort_e_t(nn)>s_t_b
                    last_job_t = sort_e_t(nn);
                    prev_job(cur_ws) = sort_e_t_ind(nn);
                else
                    last_job_t = s_t_b;
                    %if prev_same_ws_ex>0
                    %   if sort_e_t(nn)> end_time_ex(prev_same_ws_ex,  prev_job(cur_ws)) % cur_job
                    %prev_job(cur_ws) = I_sort(n);
                    prev_job(cur_ws) = sort_e_t_ind(nn);
                    %   end
                    %end
                end
                %end
            end
            
            s_e_ind=s_e_ind+1;
            if s_e_ind > q_ws %q_ws  mq_u(cur_ws)
                break;
            end
            if last_job_t == s_t_b
                break;
            end
            if last_job_t == t_last_job_t
                break;
            end
            t_last_job_t = last_job_t;
        end %flag
        %-----------------------------------------------------------------------------------------
        s_e_ind=s_e_ind-1;
    else
        last_job_t = 0;
    end
    
    
    if  cur_job_state(cur_job) ==1 %проверим можно ли начать работу
        for jj=1:n_id
            if strcmp(db_data{id_n(jj),3},'S01::ST') || strcmp(db_data{id_n(jj),3},'CRM::TX') % || strcmp(db_data{id_n(jj),3},'CRM::RS')
                if last_j_full(cur_job) == 1
                    start_time_full(last_j_full(cur_job), cur_job) = last_job_t;
                    if prev_job(cur_ws) == 0 % 1-я работа вообще
                        end_time_full(last_j_full(cur_job) , cur_job) = last_job_t + db_data{id_n(jj),16}; %q_ws %*db_data{id_n(jj),20} умножать на количество как оказалось не надо!
                    else %работа не 1-я необх учесть переналадку (в секундах)
                        end_time_full(last_j_full(cur_job) , cur_job) = last_job_t + db_data{id_n(jj),16} + 60*su{cur_ws_ex}(cur_job, prev_job(cur_ws)); %*db_data{id_n(jj),20} умножать на количество как оказалось не надо!
                    end    %q_ws
                     start_time_ex(cur_ws_ex, cur_job) = start_time_full(last_j_full(cur_job), cur_job);  
                else
                    start_time_full(last_j_full(cur_job), cur_job) = end_time_full(last_j_full(cur_job) -1, cur_job);
                    end_time_full(last_j_full(cur_job) , cur_job) = end_time_full(last_j_full(cur_job) -1, cur_job) + db_data{id_n(jj),16}; %q_ws %*db_data{id_n(jj),20} умножать на количество как оказалось не надо!
                end
                
                last_j_full(cur_job) = last_j_full(cur_job) + 1;
            else
                break;
            end
        end
        %для проверки (тест)
        %tmp1 = last_job_t + 60*p(cur_ws_ex, cur_job)*q_ws; %
        tmp2 = end_time_full(last_j_full(cur_job) -1, cur_job);
        end_time(cur_ws, cur_job) = tmp2; %last_job_ws
    else
        if p(cur_ws_ex, cur_job)>0  || p_curs_ws_cur_job>0
            if end_time(prev_ws, cur_job) == last_job_t %заплатка
                if prev_job(cur_ws)>0
                    if last_job_t - su{cur_ws_ex}(cur_job, prev_job(cur_ws)) < end_time_ex(cur_ws_ex, prev_job(cur_ws)) && last_job_t > end_time_ex(cur_ws_ex, prev_job(cur_ws))
                        last_job_t = end_time_ex(cur_ws_ex, prev_job(cur_ws));
                    else
                        last_job_t = last_job_t - su{cur_ws_ex}(cur_job, prev_job(cur_ws));
                    end
                end
            end
            %для полного расписания важно на каком этапе находится работа
            search_sym3 = '';
            search_sym4 = '';
            search_sym5 = '';
            fl_crm_tx = 0;
            switch cur_job_state(cur_job)
                case {10,21,32,43,54,65} %3
                    search_sym1 = 'CTI::ST';
                    search_sym2 = 'CTI::TX';
                case {1,3,4,6,7,9,11,13,14,16,17,23,25,26,28,29,31,33,35,36,38,39,45,47,48,50,51,53,55,57,58,60,61} %6
                    search_sym1 = 'CRM::TX';
                    %search_sym2 = 'S01::ST';
                    fl_crm_tx = 1;
                case {2,5,8,12,15,18,24,27,30,34,37,40,46,49,52,56,59,62} %1
                    search_sym1 = 'CRM::RS'; %
                case {19,41,63}% заплатка на случай 1-3-1-3-1
                    search_sym1 = 'CRM::RS'; %
                    search_sym2 = 'CRM::TX'; %
                    search_sym3 = 'CTI::ST';
                    search_sym4 = 'CTI::TX';
                case {20,42,64}
                    search_sym1 = 'EBC::FX';
                    search_sym2 = 'EBC::TX';
                case {22,44,66}
                    search_sym1 = 'G03::CU';
                    search_sym2 = 'G03::IN';
                    search_sym3 = 'G03::PA';
                    search_sym4 = 'G03::PK';
                    search_sym5 = 'FIN03::ST';
                case 67
                    %все оставшиеся
            end
            ccount=1;
            fl_en = 0;
            fl_en_STI = 0;
            fl_en_CRM_1 = 0;
            fl_en_CRM_2 = 0;
            fl_en_CRM_3 = 0;
            CTI_en = 0;
            for jj=last_j_full(cur_job):n_id
                if strcmp(db_data{id_n(jj),3},search_sym1) || strcmp(db_data{id_n(jj),3},search_sym2) || strcmp(db_data{id_n(jj),3},search_sym3) || ...
                        strcmp(db_data{id_n(jj),3},search_sym4) || strcmp(db_data{id_n(jj),3},search_sym5) || cur_job_state(cur_job) == 67
                    if ccount == 1
                        if prev_job(cur_ws) > 0 %первая работа на РЦ?   
                            start_time_full(last_j_full(cur_job), cur_job) = max(end_time(prev_ws, cur_job) - 60*su{cur_ws_ex}(cur_job, prev_job(cur_ws)), last_job_t);
                            end_time_full(last_j_full(cur_job) , cur_job) = max(end_time(prev_ws, cur_job) - 60*su{cur_ws_ex}(cur_job, prev_job(cur_ws)), last_job_t) +...
                                db_data{id_n(jj),16} + 60*su{cur_ws_ex}(cur_job, prev_job(cur_ws)); %q_ws %*db_data{id_n(jj),20} умножать на количество как оказалось не надо!
                        else
                            start_time_full(last_j_full(cur_job), cur_job) =max(end_time(prev_ws, cur_job) , last_job_t);
                            end_time_full(last_j_full(cur_job) , cur_job) = max(end_time(prev_ws, cur_job) , last_job_t) +...
                                db_data{id_n(jj),16}; %q_ws %*db_data{id_n(jj),20} умножать на количество как оказалось не надо!
                        end
                        start_time_ex(cur_ws_ex, cur_job) = start_time_full(last_j_full(cur_job), cur_job); 
                    else
                        if (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym3) 
                            CTI_en = CTI_en + 1;
                        end    
                        if (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym3) && fl_en_STI == 0 && CTI_en > 2 
                            %заплатка на случай 1-3-1-3-1 - STI - вход
                            st_STI = end_time_full(last_j_full(cur_job)-1, cur_job);
                            fl_en_STI = 1;
                        elseif (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym4) && fl_en_STI == 1 
                            %заплатка на случай 1-3-1-3-1 - STI - выход
                            end_STI = end_time_full(last_j_full(cur_job)-1, cur_job);
                        elseif  (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym1)  && fl_en_STI == 1 && fl_en_CRM_1 == 0 
                            %заплатка на случай 1-3-1-3-1 - CRM - вход №1
                            st_CRM_1 = end_time_full(last_j_full(cur_job)-1, cur_job);
                            fl_en_CRM_1 = 1;
                        elseif  (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym2)  && fl_en_CRM_1 == 1 
                            %заплатка на случай 1-3-1-3-1 - CRM - выход - №1
                            end_CRM_1 = end_time_full(last_j_full(cur_job)-1, cur_job);
                        elseif  (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym1)  && fl_en_CRM_1 == 1 && fl_en_CRM_2 == 0 
                            %заплатка на случай 1-3-1-3-1 - CRM - вход №2
                            st_CRM_2 = end_time_full(last_j_full(cur_job)-1, cur_job);
                            fl_en_CRM_2 = 1;
                        elseif  (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym2)  && fl_en_CRM_2 == 1 
                            %заплатка на случай 1-3-1-3-1 - CRM - выход - №2
                            end_CRM_2 = end_time_full(last_j_full(cur_job)-1, cur_job);    
                        elseif  (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym1)  && fl_en_CRM_2 == 1 && fl_en_CRM_3 == 0 
                            %заплатка на случай 1-3-1-3-1 - CRM - вход №3
                            st_CRM_3 = end_time_full(last_j_full(cur_job)-1, cur_job);
                            fl_en_CRM_3 = 1;
                        elseif  (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && strcmp(db_data{id_n(jj),3},search_sym2)  && fl_en_CRM_3 == 1 
                            %заплатка на случай 1-3-1-3-1 - CRM - выход - №3
                            end_CRM_3 = end_time_full(last_j_full(cur_job)-1, cur_job);        
                        end    
                        start_time_full(last_j_full(cur_job), cur_job) = end_time_full(last_j_full(cur_job)-1, cur_job);
                        end_time_full(last_j_full(cur_job) , cur_job) = end_time_full(last_j_full(cur_job)-1, cur_job) + db_data{id_n(jj),16}; %q_ws %*db_data{id_n(jj),20} умножать на количество как оказалось не надо!
                    end
                    last_j_full(cur_job) = last_j_full(cur_job) + 1;
                    ccount=ccount+1;
                    fl_en = 1;
                    if fl_crm_tx ==1 % тележки CRM:TX сами по себе (не надо склеивать)
                        break;
                    end
                else
                    if fl_en == 1
                        break;
                    end    
                end
            end
            %tmp1 = last_job_t + 60*p(cur_ws_ex, cur_job)*q_ws; %
            if last_j_full(cur_job)>1
                tmp2 = end_time_full( last_j_full(cur_job) - 1, cur_job);
            else
                tmp2 = end_time_full(1, cur_job);
            end    
            end_time(cur_ws, cur_job) = tmp2; %last_job_ws
        else
            end_time(cur_ws, cur_job) = end_time(prev_ws, cur_job);
        end
    end
    if q_ws==1
        if p(cur_ws_ex, cur_job)>0  || p_curs_ws_cur_job>0
            s_e_ws{cur_ws}(1,m_j_s_e_ws{cur_ws}) = start_time_ex(cur_ws_ex, cur_job);
            %заплатка на случай 1-3-1-3 
            if (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && fl_en_STI == 1
                s_e_ws{cur_ws}(2,m_j_s_e_ws{cur_ws}) = st_STI;
                %займем склад ... %всегда 1-я ячейка %потом перебросим
                s_e_ws{3}{1}(1,m_j_s_e_ws{3}(1)) = st_STI;
                s_e_ws{3}{1}(2,m_j_s_e_ws{3}(1)) = end_STI;
                s_e_ws{3}{1}(3,m_j_s_e_ws{3}(1)) = cur_job;
                m_j_s_e_ws{3}(1)=m_j_s_e_ws{3}(1)+1;
            else    
                s_e_ws{cur_ws}(2,m_j_s_e_ws{cur_ws}) = end_time(cur_ws, cur_job);
            end    
            s_e_ws{cur_ws}(3,m_j_s_e_ws{cur_ws}) = cur_job;
            m_j_s_e_ws{cur_ws}=m_j_s_e_ws{cur_ws}+1;
            %заплатка на случай 1-3-1-3-1
            if (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && fl_en_CRM_1 == 1
                s_e_ws{cur_ws}(1,m_j_s_e_ws{cur_ws}) = st_CRM_1;
                s_e_ws{cur_ws}(2,m_j_s_e_ws{cur_ws}) = end_CRM_1;
                s_e_ws{cur_ws}(3,m_j_s_e_ws{cur_ws}) = cur_job;
                m_j_s_e_ws{cur_ws}=m_j_s_e_ws{cur_ws}+1;
            end    
            if (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && fl_en_CRM_2 == 1
                s_e_ws{cur_ws}(1,m_j_s_e_ws{cur_ws}) = st_CRM_2;
                s_e_ws{cur_ws}(2,m_j_s_e_ws{cur_ws}) = end_CRM_2;
                s_e_ws{cur_ws}(3,m_j_s_e_ws{cur_ws}) = cur_job;
                m_j_s_e_ws{cur_ws}=m_j_s_e_ws{cur_ws}+1;
            end    
            if (cur_job_state(cur_job)==19 || cur_job_state(cur_job)==41 || cur_job_state(cur_job)==63) && fl_en_CRM_3 == 1
                s_e_ws{cur_ws}(1,m_j_s_e_ws{cur_ws}) = st_CRM_3;
                s_e_ws{cur_ws}(2,m_j_s_e_ws{cur_ws}) = end_CRM_3;
                s_e_ws{cur_ws}(3,m_j_s_e_ws{cur_ws}) = cur_job;
                m_j_s_e_ws{cur_ws}=m_j_s_e_ws{cur_ws}+1;
            end    
        end
    else
        if p(cur_ws_ex, cur_job)>0  || p_curs_ws_cur_job>0
            s_e_ws{cur_ws}{s_e_ind}(1,m_j_s_e_ws{cur_ws}(s_e_ind)) = start_time_ex(cur_ws_ex, cur_job);
            s_e_ws{cur_ws}{s_e_ind}(2,m_j_s_e_ws{cur_ws}(s_e_ind)) = end_time(cur_ws, cur_job);
            s_e_ws{cur_ws}{s_e_ind}(3,m_j_s_e_ws{cur_ws}(s_e_ind)) = cur_job;
            m_j_s_e_ws{cur_ws}(s_e_ind)=m_j_s_e_ws{cur_ws}(s_e_ind)+1;
        end
    end
    %if p(cur_ws_ex, cur_job)>0
    if p(cur_ws_ex, cur_job)==0 && p_curs_ws_cur_job==0
        start_time_full(cur_ws_ex, cur_job) = end_time(cur_ws, cur_job);
    end
    end_time_ex(cur_ws_ex, cur_job) = end_time(cur_ws, cur_job);
    %end
    prev_job(cur_ws) = cur_job;
    
    %stU=[];
    %last_j=0;
    
    j=j+1;
end
%end
%снизим размерность
[i,~] = find(start_time_full);
start_time_full=start_time_full(1:max(i),:);
[i,~] = find(end_time_full);
end_time_full=end_time_full(1:max(i),:);
% [i,j] = find(end_time_ex);
% end_time_ex = end_time_ex(:,1:max(j));
% end_time_ex = end_time_ex(1:max(i),:);
%
% [i,j] = find(end_time);
% end_time = end_time(:,1:max(j));
% end_time = end_time(1:max(i),:);

% d1=d(1:max(j));

%T=max(0,end_time(ma(numel(ma)),:)-d1);
%T=max(0,end_time(ma(numel(ma)),:)-d);
%sumT = sum(T);
%Cmax1 = end_time(prev_ws, cur_job);
%Cmax2 = end_time(cur_ws, cur_job);
end

