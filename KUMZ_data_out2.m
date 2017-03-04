function [data_out,data_out2] = KUMZ_data_out2( start_time_full, end_time_full, db_data, session_id, id_t, id_in)

[~,nc]=size(end_time_full);
data_out = {};
%[~,I]=sort(end_time_full(1,:));
data_out2 = [];
st_b = now;

cc=1;

id_in_u = unique(id_in,'stable');
id_in_u = nonzeros(id_in_u);%find(id_in_u~=0);

for ii=1:nc
    cur_job = ii;
    f_srch = id_in_u(cur_job);
    id_n = find(id_t==f_srch); %cur_job
    % id_n = find(id_in==I(ii));
       
     %      if ii==1
     %         st_t = st_b;
     %      else
     %          st_t = end_t_nc;
     %      end
     for jj=1:numel(id_n)
         %          if jj==1
         %              st_t = st_t;
         %              end_t_nc = st_b + end_time_ex_full(jj,I(ii))/3600/24;
         %          else
         %              st_t = end_t;
         %          end
        st_t =  st_b + start_time_full(jj, ii)/3600/24;
        end_t = st_b + end_time_full(jj, ii)/3600/24;
        %         if I(ii) ==1
        %             end_time_full(jj,I(ii))
        %          end
        
        if end_t==0
            break;
        end    
        data_out{cc,1} = {session_id num2str(db_data{id_n(jj),1}) datestr(st_t,'yyyy-mm-dd HH:MM:SS') datestr(end_t,'yyyy-mm-dd HH:MM:SS')}; %datestr
        data_out2(cc,:) = [db_data{id_n(jj),1} st_t end_t];
        cc = cc + 1;
     end
end   

end

