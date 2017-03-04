s1=tic;
dpath = 'C:\Program Files\Microsoft JDBC Driver 6.0 for SQL Server\sqljdbc_6.0\rus\sqljdbc4.jar';

javaaddpath(dpath,'-end')

conn = database('APPPS','admin','adminkonsom',...
'Vendor','Microsoft SQL Server',...
'Server','172.18.20.101','PortNumber',1433);

sqlquery = '{call _auto_data_get_new(null, null, ''logic'',''oper'',''2017-03-03 15:50:00'')}';
curs = exec(conn,sqlquery);
curs = fetch(curs);
db_data = curs.Data;
e1=toc(s1);

% готовим данные
s2=tic;
%[db_data,id_t,id_in,r_in,d_in,p_in,su1_in,su2_in] = test_data_prep4c(db_data);
[db_data1,id_t_s,id_in_s,r_in_s,d_in_s,p_in_s,su1_in_s,su2_in_s] = test_data_prep5(db_data);
e2=toc(s2);

% снижаем размерность
s3=tic;
[id_t,id_in,r_in,d_in,p_in,su1_in,su2_in,S_out,S_umerge,S_merge] = low_dim(id_t_s,id_in_s,r_in_s,d_in_s,p_in_s,su1_in_s,su2_in_s);
e3=toc(s3);

s4=tic;
S2min = model9p(r_in,d_in,p_in,su1_in,su2_in);
e4=toc(s4);

s5=tic;
session_id = 1;
[data_out,data_out2] =  make_full_sched(S2min,S_umerge,S_merge,S_out,r_in_s,d_in_s,p_in_s,su1_in_s,su2_in_s,db_data1,session_id,id_t_s, id_in_s);
e5=toc(s5);

s6=tic;
[nr,~]=size(data_out);
colnames = {'date_st','date_fn'};
tablename = '_PO';

% sqlquery = 'select * from FK_PO_PO_MOUNT';
% curs = exec(conn,sqlquery);
% curs = fetch(curs);
% po_data = curs.Data;

%очистка
whereclause = '';
update(conn,tablename,colnames,[NaN NaN],whereclause);

for i=1:nr
    %whereclause = strcat('where OP_ID =' ,data_out{i}(2),' and InsTime=','''','2017-02-28 10:00:00',''''); 
    whereclause = strcat('where OP_ID =' ,data_out{i}(2)); 
    update(conn,tablename,colnames,data_out{i}(3:4),whereclause);
end 
e1=toc(s6);