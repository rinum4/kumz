function [data_out,data_out2] = make_full_sched(S2min,S_umerge,S_merge,S_out,r_in_s,d_in_s,p_in_s,su1_in_s,su2_in_s,db_data1,session_id,id_t_s, id_in_s)
% формируем полное расписание 
[~,n]=size(p_in_s);

%восстанавливаем полное расписание
%шаг1. замена 1-n на S_umerge
nsz = numel(S2min);
S3min = zeros(1,nsz);
for i=1:nsz
    S3min(i)=S_umerge(S2min(i));
end

%шаг2. добавление одинаковых работ из S_merge
S4min = [];
for i=1:nsz
    ind = find(S_merge(1,:)==S3min(i));
    if ~isempty(ind)
        S4min=[S4min S3min(i) S_merge(2,ind)];
    else
        S4min=[S4min S3min(i)];
    end
end

%шаг3. добавление вперед работ из S_out
S5min = [];
for i=1:numel(S_out)
    S5min = [S5min ones(1,19)*S_out(i)];
end
S5min = [S5min S4min];

%добавим тележки
ma3=[6 1 6 6 1 6 6 1 6 3 6 1 6 6 1 6 6 1 6 2 3 4 6 1 6 6 1 6 6 1 6 3 6 1 6 6 1 6 6 1 6 2 3 4 6 1 6 6 1 6 6 1 6 3 6 1 6 6 1 6 6 1 6 2 3 4 5]; % последовательность РЦ 31
mq3 =[14 1 14 14 1 14 14 1 14 473 14 1 14 14 1 14 14 1 14 18 473 1 14 1 14 14 1 14 14 1 14 473 14 1 14 14 1 14 14 1 14 18 473 1 14 1 14 14 1 14 14 1 14 473 14 1 14 14 1 14 14 1 14 18 473 1 2]; % кол-во машин в РЦ
mq_u3 = [1 18 473 1 2 14];
m3=6;

su3=cell(1,67);
%ast1= tic;
for ii=1:67 %parfor % бессмысленнен
    su3{ii}=zeros(n,n);
end
%aend1 = toc(ast1);

% валки
su3{2}=su1_in_s(1:n,1:n);
su3{5}=su1_in_s(1:n,1:n);
su3{8}=su1_in_s(1:n,1:n);
su3{12}=su1_in_s(1:n,1:n);
su3{15}=su1_in_s(1:n,1:n);
su3{18}=su1_in_s(1:n,1:n);
su3{24}=su1_in_s(1:n,1:n);
su3{27}=su1_in_s(1:n,1:n);
su3{30}=su1_in_s(1:n,1:n);
su3{34}=su1_in_s(1:n,1:n);
su3{37}=su1_in_s(1:n,1:n);
su3{40}=su1_in_s(1:n,1:n);
su3{46}=su1_in_s(1:n,1:n);
su3{49}=su1_in_s(1:n,1:n);
su3{52}=su1_in_s(1:n,1:n);
su3{56}=su1_in_s(1:n,1:n);
su3{59}=su1_in_s(1:n,1:n);
su3{62}=su1_in_s(1:n,1:n);

% ножи G03
su3{22}=su2_in_s{1}(1:n,1:n);
su3{44}=su2_in_s{2}(1:n,1:n);
su3{66}=su2_in_s{3}(1:n,1:n);

%clear p3
%p то же необходимо расширить
p3(1,:)=ones(1,n);
p3(2,:)=p_in_s(1,:)/3;
p3(3,:)=ones(1,n);
p3(4,:)=ones(1,n);
p3(5,:)=p_in_s(1,:)/3;
p3(6,:)=ones(1,n);
p3(7,:)=ones(1,n);
p3(8,:)=p_in_s(1,:)/3;
p3(9,:)=ones(1,n);
p3(10,:)=p_in_s(2,:);
p3(11,:)=ones(1,n);
p3(12,:)=p_in_s(3,:)/3;
p3(13,:)=ones(1,n);
p3(14,:)=ones(1,n);
p3(15,:)=p_in_s(3,:)/3   ;
p3(16,:)=ones(1,n);
p3(17,:)=ones(1,n);
p3(18,:)=p_in_s(3,:)/3;
p3(19,:)=ones(1,n);
p3(20,:)=p_in_s(4,:);
p3(21,:)=p_in_s(5,:);
p3(22,:)=p_in_s(6,:);
p3(23,:)=ones(1,n);
p3(24,:)=p_in_s(7,:)/3;
p3(25,:)=ones(1,n);
p3(26,:)=ones(1,n);
p3(27,:)=p_in_s(7,:)/3;
p3(28,:)=ones(1,n);
p3(29,:)=ones(1,n);
p3(30,:)=p_in_s(7,:)/3;
p3(31,:)=ones(1,n);
p3(32,:)=p_in_s(8,:);
p3(33,:)=ones(1,n);
p3(34,:)=p_in_s(9,:)/3;
p3(35,:)=ones(1,n);
p3(36,:)=ones(1,n);
p3(37,:)=p_in_s(9,:)/3;
p3(38,:)=ones(1,n);
p3(39,:)=ones(1,n);
p3(40,:)=p_in_s(9,:)/3;
p3(41,:)=ones(1,n);
p3(42,:)=p_in_s(10,:);
p3(43,:)=p_in_s(11,:);
p3(44,:)=p_in_s(12,:);
p3(45,:)=ones(1,n);
p3(46,:)=p_in_s(13,:)/3;
p3(47,:)=ones(1,n);
p3(48,:)=ones(1,n);
p3(49,:)=p_in_s(13,:)/3;
p3(50,:)=ones(1,n);
p3(51,:)=ones(1,n);
p3(52,:)=p_in_s(13,:)/3;
p3(53,:)=ones(1,n);
p3(54,:)=p_in_s(14,:);
p3(55,:)=ones(1,n);
p3(56,:)=p_in_s(15,:)/3;
p3(57,:)=ones(1,n);
p3(58,:)=ones(1,n);
p3(59,:)=p_in_s(15,:)/3;
p3(60,:)=ones(1,n);
p3(61,:)=ones(1,n);
p3(62,:)=p_in_s(15,:)/3;
p3(63,:)=ones(1,n);
p3(64,:)=p_in_s(16,:);
p3(65,:)=p_in_s(17,:);
p3(66,:)=p_in_s(18,:);
p3(67,:)=p_in_s(19,:);

%Smin_ex = KUMZ_19_67(S2min,n);
Smin_ex = KUMZ_19_67(S5min,n);

%[end_time_ex_12_SACO_12_full]= KUMZ_data_load(p,S, n, m, ma, d, su, mq, db_data,id_in);
%функция с учетом r
%[end_time_ex_12_SACO_12_full]= KUMZ_data_load2(p,S, n, m, ma, d, su, mq, db_data,id_in,r);
%[start_time_12_S_12_full, end_time_12_S_12_full]= KUMZ_data_load3(p,Smin, n, m, ma, d, su, mq, r, mq_u, db_data, id_t, id_in);
%[start_time_12_S_12_full, end_time_12_S_12_full]= KUMZ_data_load3c(p, Smin, n, m, ma, d, su, mq, r, mq_u, db_data, id_t, id_in);
[start_time_12_S_12_full, end_time_12_S_12_full]= KUMZ_data_load3c(p3, Smin_ex, n, m3, ma3, d_in_s, su3, mq3, r_in_s, mq_u3, db_data1, id_t_s, id_in_s);

%session_id = 7;
[data_out,data_out2] = KUMZ_data_out2( start_time_12_S_12_full, end_time_12_S_12_full, db_data1, session_id, id_t_s, id_in_s);


end

