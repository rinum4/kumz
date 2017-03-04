function [sumT, Cmax1, Cmax2, end_time_ex] = KUMZTargFun19fga(p, Sga, n, m, ma, d, su, mq, r, mq_u)
% ������������� ������� ��� Fm | reentr | Sum{Tmax}
%�����! � ������ ������� ��������� ������ �������� su
%�����! � ������ ������� �������� ������������ ��������������� ��
%�����! � ������ ������� �������� ������������ r - ������ �� �����
%�����! ��� ������ ����� �� ��������
%�����! ��������� ������ ����� �� ������� ������
%�����! ����� ���������� �������� ����������� �������

%�����: ����� S - ������ ���������� M1-M2-M3-...M1-M2-M3-...

%�����: ������� ��� ga ���������, � ��� ������� ��������������� � ������
nk = numel(Sga)/19;
S=zeros(1,nk);
for k=1:nk
    indS=find(Sga(1+(k-1)*n:k*n),1);
    if ~isempty(indk)
        
    indk = find(Uk,S(k));
    if ~isempty(indk)
        
    end
    end    
end

cur_job_state=zeros(n,1); %��������� ������� ������
ma_ex = 1:numel(ma);

s_e_ws = cell(m,1);
m_j_s_e_ws = cell(m,1);
for i=1:m %parfor
    s_e_ws{i}=zeros(3,n*60);
    m_j_s_e_ws{i} = 1;
    if mq_u(i)>1 %��� ������������� �� �������� ������ ������ ���������
        s_e_ws{i}=cell(mq_u(i),1);
        m_j_s_e_ws{i} = ones(mq_u(i),1);
        for j=1:mq_u(i)
            s_e_ws{i}{j}=zeros(3,n*60);
        end
    end
end

%minq_job_ws = ones(m,1); %���� ������������� ��, �� �������� ������ ���� �� � ��� ����
%minq_job_ws_ex = ones(numel(ma),1); %���� ������������� ��, �� �������� ������ ���� �� � ��� ����

last_job_ws_q = cell(m,1);
for i=1:m %parfor
    if mq_u(i)>1 %��� ������������� �� �������� ������ ������ ���������
        last_job_ws_q{i}=zeros(mq_u(i),n); %��������� ��������� ������ �� ���� �� (�����������, � ��������)
    end    
end

%start_job_ws_q = cell(m,1);
% last_job_ws_q_ex = cell(numel(ma),1);
% for i=1:numel(ma) %parfor
%     if mq(i)>1 %��� ������������� �� �������� ������ ������ ���������
% %       %start_job_ws_q{i}=zeros(mq(i),n);
%         last_job_ws_q_ex{i}=zeros(mq(i),n); %��������� ��������� ������ �� ���� �� (�����������, � ��������)
%     end
% end

%��� ����������� ������ ����� �������� ������
end_time = zeros(max(ma),n);
start_time_ex = zeros(numel(ma),n);
end_time_ex = zeros(numel(ma),n);

[~,nc] = size(S);
prev_job = zeros(max(ma),1);

%for i=1:nr %������ ���� nr=1 ������
%for j=1:nc
j=1;
%last_j = 0;
while j<nc+1 %for ������ �� while ��� �������������
    cur_job = S(1,j);
    %     if cur_job==9
    %         stop=1;
    %     end
    cur_job_state(cur_job) =  cur_job_state(cur_job) + 1;
    
    cur_ws_ex = ma_ex(cur_job_state(cur_job));
    
    cur_ws = ma(cur_job_state(cur_job));
    %     if cur_ws==2 && p(cur_ws_ex, cur_job)>0 %&& cur_job==8
    %         stop=1;
    %     end
    
    if cur_ws==0
        cur_ws = m;
    end
    
    q_ws = mq(cur_ws_ex); %���������� ����� � ��
    
    if cur_job_state(cur_job)>1
        prev_ws = ma(cur_job_state(cur_job)-1);
        prev_ws_ex = cur_ws_ex - 1;
    else
        prev_ws = 0;
        prev_ws_ex = 0;
    end
    
    %�������� ���������� ����� ��������� ������ �� ������ ��
    if p(cur_ws_ex, cur_job)>0
        if prev_ws_ex>0
            s_t_b = end_time_ex(prev_ws_ex, cur_job); %end_time(prev_ws, cur_job);
        else
            s_t_b = r(cur_job);
        end
        
        s_t_e = s_t_b + p(cur_ws_ex, cur_job)*q_ws; %��� ����������� �� �����! + su{cur_ws_ex}(cur_job, prev_job(cur_ws))
        
        last_job_t = 0;
        %--------------------------------------------------------------------------------------------------------------------------------------
        %�����! ������ ���� ����� ��� ������������� ��!
        flag=1;
        cc = 1; % ������ 
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
                [sort_e_t, I_sort] = sort(s_e_ws{cur_ws}{cc}(2,:));
                sort_e_t_prev = s_e_ws{cur_ws}{cc}(1,I_sort);
                sort_e_t_ind = s_e_ws{cur_ws}{cc}(3,I_sort);
            end
            
            ind = find(sort_e_t>s_t_e,1);
            nn = numel(sort_e_t);
            if ~isempty(ind)
                fl=0;
                
                for jj = ind-1:nn-1 %ind - 1:n-1
                    %����������� ����������� ���! �.�. ������ � ���� ������
                    %�� ����� ����� ������ ����� ����
                    %�����! ���������� ������ ����� � 2�� �������������! �
                    %������ � ������
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
            
           cc=cc+1;
           if cc > q_ws %q_ws  mq_u(cur_ws) 
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
        cc=cc-1;
        %-----------------------------------------------------------------------------------------
    else
        last_job_t = 0;
    end
    
    if  cur_job_state(cur_job) ==1 %�������� ����� �� ������ ������
        if p(cur_ws_ex, cur_job)>0
            start_time_ex(cur_ws_ex, cur_job) = last_job_t;
            if prev_job(cur_ws) == 0
                end_time(cur_ws, cur_job) = last_job_t + p(cur_ws_ex, cur_job)*q_ws;
            else
                end_time(cur_ws, cur_job) = last_job_t  + p(cur_ws_ex, cur_job)*q_ws + ...
                    su{cur_ws_ex}(cur_job, prev_job(cur_ws));
            end
        end
    else
        if p(cur_ws_ex, cur_job)>0
            if end_time(prev_ws, cur_job) == last_job_t %��������
              if prev_job(cur_ws)>0
                if last_job_t - su{cur_ws_ex}(cur_job, prev_job(cur_ws)) < end_time_ex(cur_ws_ex, prev_job(cur_ws)) && last_job_t > end_time_ex(cur_ws_ex, prev_job(cur_ws))
                    last_job_t = end_time_ex(cur_ws_ex, prev_job(cur_ws));
                else
                    last_job_t = last_job_t - su{cur_ws_ex}(cur_job, prev_job(cur_ws));
                end
              end
            end
            if prev_job(cur_ws) > 0
                start_time_ex(cur_ws_ex, cur_job) = max(end_time(prev_ws, cur_job) - su{cur_ws_ex}(cur_job, prev_job(cur_ws)), last_job_t);
                end_time(cur_ws, cur_job) = max(end_time(prev_ws, cur_job) - su{cur_ws_ex}(cur_job, prev_job(cur_ws)), last_job_t) + ...
                    p(cur_ws_ex, cur_job)*q_ws + su{cur_ws_ex}(cur_job, prev_job(cur_ws));
            else
                start_time_ex(cur_ws_ex, cur_job) = max(end_time(prev_ws, cur_job) , last_job_t);
                end_time(cur_ws, cur_job) = max(end_time(prev_ws, cur_job) , last_job_t) + ...
                    p(cur_ws_ex, cur_job)*q_ws;
            end
        else
            end_time(cur_ws, cur_job) = end_time(prev_ws, cur_job);
        end
    end
    if q_ws==1
        if p(cur_ws_ex, cur_job)>0
            s_e_ws{cur_ws}(1,m_j_s_e_ws{cur_ws}) = start_time_ex(cur_ws_ex, cur_job);
            s_e_ws{cur_ws}(2,m_j_s_e_ws{cur_ws}) = end_time(cur_ws, cur_job);
            s_e_ws{cur_ws}(3,m_j_s_e_ws{cur_ws}) = cur_job;
            m_j_s_e_ws{cur_ws}=m_j_s_e_ws{cur_ws}+1;
        end
    else
        if p(cur_ws_ex, cur_job)>0
            s_e_ws{cur_ws}{cc}(1,m_j_s_e_ws{cur_ws}(cc)) = start_time_ex(cur_ws_ex, cur_job);
            s_e_ws{cur_ws}{cc}(2,m_j_s_e_ws{cur_ws}(cc)) = end_time(cur_ws, cur_job);
            s_e_ws{cur_ws}{cc}(3,m_j_s_e_ws{cur_ws}(cc)) = cur_job;
            m_j_s_e_ws{cur_ws}(cc)=m_j_s_e_ws{cur_ws}(cc)+1;
            
            %start_job_ws_q{cur_ws_ex}(minq_job_ws_ex(cur_ws_ex),cur_job) = start_time_ex(cur_ws_ex, cur_job);
            %last_job_ws_q{cur_ws}(minq_job_ws(cur_ws),cur_job) = end_time(cur_ws, cur_job);
            %last_job_ws_q_ex{cur_ws_ex}(minq_job_ws_ex(cur_ws_ex),cur_job) = end_time(cur_ws, cur_job);
            %             if cur_ws==2
            %                 stop=1;
            %             end
            
            %tmp_min = zeros(1,mq_u(cur_ws));
            
            %for jj=1:mq_u(cur_ws)
                % �����! ������� ������������ ������� ������ �� ����������
                % tmp_min(jj) = max(last_job_ws_q{cur_ws}(jj,:)); %
                % ����� ��������� - ����� ����������� �� ����
            %    tmp_min(jj) = max(last_job_ws_q{cur_ws}(jj,:)); %
            %end
            %             tmp_min = zeros(1,mq(cur_ws_ex));
            %             for jj=1:mq(cur_ws_ex)
            %                 tmp_min(jj) = max(last_job_ws_q_ex{cur_ws_ex}(jj,:));
            %             end
            %[~,minq_job_ws_ex(cur_ws_ex)] = min(tmp_min);
            %[~,minq_job_ws(cur_ws)] =  min(tmp_min); %minq_job_ws_ex(cur_ws_ex); - ��� ������!
        end
    end
    %if p(cur_ws_ex, cur_job)>0
    if p(cur_ws_ex, cur_job)==0
       start_time_ex(cur_ws_ex, cur_job) = end_time(cur_ws, cur_job);
    end    
    end_time_ex(cur_ws_ex, cur_job) = end_time(cur_ws, cur_job);
    %end
    prev_job(cur_ws) = cur_job;
    
    %stU=[];
    %last_j=0;
    
    j=j+1;
end
%end

%������ �����������
% [i,j] = find(end_time_ex);
% end_time_ex = end_time_ex(:,1:max(j));
% end_time_ex = end_time_ex(1:max(i),:);
%
% [i,j] = find(end_time);
% end_time = end_time(:,1:max(j));
% end_time = end_time(1:max(i),:);

% d1=d(1:max(j));

%T=max(0,end_time(ma(numel(ma)),:)-d1);
T=max(0,end_time(ma(numel(ma)),:)-d);
sumT = sum(T);
Cmax1 = end_time(prev_ws, cur_job);
Cmax2 = end_time(cur_ws, cur_job);
end

