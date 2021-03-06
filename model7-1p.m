%function [data_out]=model6p(id_t,id_in,r_in,d_in,p_in,su1_in,su2_in,db_data)
% F4 |(1,3,1),(2,3,4),(1,3,1),(2,3,4),(1,3,1),(2,3,4,5) : (1-2-(1-2)-1-2)

flag1 = 0;
if flag1==1
    clc
    clear
    
    load('db_Data5.mat')
    load('id_data5.mat')
    load('p_data5.mat');
    load('su1_data5.mat');
    load('su2_data5.mat');
    load('r_data5.mat');
    load('d_data5.mat')
end

%global record;
%global xmin;
%global ccount; %������� ������� ��� �����

% global F_max;
% global a_max;

%record = inf;
%xmin=[];
%ccount=0;

n=20; %10
flag_branch = 0;
n_max=30;
m=5;

%p = ceil(rand(19,5)*10);        %������ 1
%p=p_in;
% p = [54  53    52  50  50  76  50  51  51  53;  %1 - 1
%     9     9      9     0    0     0   0     0    0    9;      %3  2
%     27   26   26   0    0     0   0     0    0    53;   %1   3
%     17   17   34   36 36   34  35  35  35  35;  %2   4
%     9     9     9      9   9     9    9     9    9    9;   %3     5
%     6     5     5      5   4     0    0     0    0    5;     %4   6
%     28   26   26  50  49  76  49   50  50  53;   %1   7
%     9     9      9    0    0     0     0    0    0    9;     %3   8
%     58   54   54  0    0     0     0    0    0    55;     %1 9
%     17   17   34  36  36  35  35   35  35  35;  %2   10
%     9     9     9     9    9     9    9     9    9    9;   %3     11
%     0     0     0     0    0     0    0     0    5    0;   %4     12
%     63   56  56   52  50   51  33  33  52  0;   %1     13
%     9     9     9     9    9      9   9     9    9    0;      %3  14
%     0     0     0     18  17   52  33  34  0    0;     %1   15
%     17   17  34    36 36   34  38  35  35  0;  %2      16
%     9     9     9      9   9     9     9    9    9    0;   %3     17
%     33   0    30    0    0     0     0    0    0    0;   %4    18
%     0     41   0     48   39  36  35  38  46  40];   %5 19

%������ ����������� ����� ����� �������
%p=round(p_in(:,1:n)*100)/100;

p=p_in(:,1:n);

su=cell(1,19);

% ��������� ������ �� ��� �����
%if isempty(gcp('nocreate')) %���
%    parpool('local',4); % �������� 4 ����������
%end
    
%�������������
parfor ii=1:19 %parfor
    su{ii}=zeros(n_max,n_max);
end    
% �����

su{1}=su1_in;
su{3}=su1_in;
su{7}=su1_in;
su{9}=su1_in;
su{13}=su1_in;
su{15}=su1_in;
% ����
su{6}=su2_in{1};
su{12}=su2_in{2};
su{18}=su2_in{3};

%r=r_in(1:n);
%���� �� ������ � ����������� �������� r
r = zeros(1,n); %release dates - ��������� � ��������

%r = ceil(rand(1,n)*2300);
%r = [150 0 195 214 156 174 171 90 150 39];

%����
%r(2) = 600;
%r(5) = 600;

%d = ones(1,n)*4700; %due dates
%d = ceil(rand(1,n)*4700); %4700 9300
%d = [3082 168 3991 4390 3191 3562 3493 1844 3081 805];
%d = [5431 14387 14321 7229 10789 1442 1026 10086 14805 17747];
%d = [3881 462 8396 8787 4566 4551 3141 8371 3434 1035 7257 3625 2248 3757 898 1228 8762 8893 5350 556 2184 3285 7638 144 401 1572 6037 6806 6025 4194];

d = d_in(1:n);
%d= max(d,0);

ma=[1 3 1 2 3 4 1 3 1 2 3 4 1 3 1 2 3 4 5]; % ������������������ ��

mq =[1 473 1 18 473 1 1 473 1 18 473 1 1 473 1 18 473 1 2]; % ���-�� ����� � ��
mq_u = [1 18 473 1 2];

%������ � ����������� ����������
% 1. EDD
[~, I] = sort(d);
%[~, Ir] = sort(r);

% �����! KumzFullS �������������
%S = KumzFullS(I,ma); % (perm) ���������� - ������-� ������ �������:
S_EDD = KumzFullSp(I,ma); % (perm) ���������� - ������-� ������ �������:
% �����! XieTargFun7/12f �� ������������
%[T_EDD1_1,~,~, end_time_ex_EDD1_1]= XieTargFun7(p, S, n, m, ma, d, su);
%��������� ��������������� ��
%[T_EDD1,~,~, end_time_ex_EDD1]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
%��������� �������� r
%tmp
%[T_EDD1,~,~, end_time_ex_EDD1]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_EDD,~,~, end_time_ex_EDD]= KUMZTargFun19f(p, S_EDD, n, m, ma, d, su, mq, r, mq_u);
%[start_time_EDD1_full, end_time_EDD1_full]= KUMZ_data_load3(p,S, n, m, ma, d, su, mq, r, mq_u, db_data, id_in);
%������ ���������� �������� ������ ��� ����������� �������!
%[end_time_ex_EDD_full]= KUMZ_data_load(p,S, n, m, ma, d, su, mq, db_data,id_in);

% �����! �������� ��������� ��������������� �� ������� ��������
%S = KumzFullS2(I,ma); % (perm) ���������� - ������-� �� ������ ��:
%[T_EDD2_1,~,~, end_time_ex_EDD2_1] =  XieTargFun7(p, S, n, m, ma, d, su);
%[T_EDD2_2,~,~, end_time_ex_EDD2_2] = XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_EDD2_2,~,~, end_time_ex_EDD2_2] = KUMZTargFun16f(p, S, n, m, ma, d, su, mq, r);

% 2. ������������ ������������
p1=sum(p(:,1:n));
[~,I]=sort(p1,'ascend');
%S = KumzFullS(I,ma);
S_maxp = KumzFullSp(I,ma);
%S = KumzFullS2(I,ma); % (perm) ���������� - ������-� �� ������ ��:
%[T_maxp,~,~, end_time_ex_maxp] = XieTargFun7(p, S, n, m, ma, d, su);
%[T_maxp,~,~, end_time_ex_maxp]=XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_maxp,~,~, end_time_ex_maxp] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_maxp,~,~, end_time_ex_maxp] = KUMZTargFun19f(p, S_maxp, n, m, ma, d, su, mq, r, mq_u);
k1 = 0.3;
Fmaxp = (T_maxp - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_maxp)) - max(max(end_time_ex_EDD)));

%
% % ������� ������ ������ (1,3,1)-1
%a1=p(1,:);
% b1=p(2,:);
% c1=p(3,:);
%
% a=sum(a1);
% b=sum(b1);
% c=sum(c1);
%
% pj=[a1;b1]';
% S1= jhons(pj);
% Cab = WangTargFun(a1,b1,c1,S1);
% if Cab == a + c
%     xmin1 = S1;
%     record1 = Cab;
%     exit=0;
% else
%     F=WangBranchAndBoundFirst(a1,b1,c1,S1,n, a, b, c);
%     xmin1 = xmin;
%     record1 = record;
% end
%
% % ������� ������ ������ (1,3,1)-2
% a1=p(7,:);
% b1=p(8,:);
% c1=p(9,:);
%
% a=sum(a1);
% b=sum(b1);
% c=sum(c1);
%
% pj=[a1;b1]';
% S1= jhons(pj);
% Cab = WangTargFun(a1,b1,c1,S1);
% if Cab == a + c
%     xmin2 = S1;
%     record2 = Cab;
%     exit=0;
% else
%     F=WangBranchAndBoundFirst(a1,b1,c1,S1,n, a, b, c);
%     xmin2 = xmin;
%     record2 = record;
% end

%--------------------------------------------------------------------------
% 3. �������� NEH
S=1:n;
% �����! NawazHehCmax �� ������������
%[~,NEH3]=NawazHehCmax(p',m,ma,S);
%��������� �������� r
%�����! � ����� ��������� ��� ���� �� ���� ����, ��� �� ���-�� ����������
%����������� TargFun15f, ���� ������� ���
[~,NEH3]=NawazHehCmax2(p',m,ma,S,r);
%S = KumzFullS(NEH3,ma); % (perm) ���������� - ������-� ������ �������:
S_NEH_f = KumzFullSp(NEH3,ma); % (perm) ���������� - ������-� ������ �������:
%S = KumzFullS2(NEH3,ma); % (perm) ���������� - ������-� �� ������ ��:
%[T_NEH_f,~,~, end_time_ex_NEH_f]= XieTargFun7(p, S, n, m, ma, d, su);
%[T_NEH_f,~,~, end_time_ex_NEH_f]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_NEH_f,~,~, end_time_ex_NEH_f] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_NEH_f,~,~, end_time_ex_NEH_f] = KUMZTargFun19f(p, S_NEH_f, n, m, ma, d, su, mq, r, mq_u);
FNEH_f = (T_NEH_f - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_NEH_f)) - max(max(end_time_ex_EDD)));

%--------------------------------------------------------------------------
% ������� ������ ��� ������: (1-2-1-2)-1
% �. 1. ������ ������ (1-3-1)=1,(2,3,4)=2,(1-3-1)-2=1,(2,3,4)-2=2,
% ma=[1 3 1 2 3 4 1 3 1 2 3 4 1 3 1 2 3 4 5];
p1=zeros(4,n);
p1(1,:) = sum(p(1:3,:));
p1(2,:) = sum(p(4:6,:));
p1(3,:) = sum(p(7:9,:));
p1(4,:) = sum(p(10:12,:));
%
su1 = cell(1,4);
%for ii=1:4 %parfor
    su1{1}(:,:)=su{1} + su{2} + su{3}; %sum(su{ii}(1:3,:));
    su1{2}(:,:)=su{4} + su{5} + su{6}; %sum(su{ii}(4:6,:));
    su1{3}(:,:)=su{7} + su{8} + su{9}; %sum(su{ii}(7:9,:));
    su1{4}(:,:)=su{10} + su{11} + su{12}; %sum(su{ii}(10:12,:));
%end

%d1 = d - sum(p(13:19,:));
%�.�. d ��� � ������ �������� �������������
d1 = d - p(13,:) - p(14,:)*473 - p(15,:) - p(16,:)*18 - p(17,:)*473 - p(18,:) - p(14,:)*2; 

%�������� �������� NEH ��� ���������� �����������
%�����! ��� ���� �����������
S_NEH=1:n;
m1=2;
ma1=[1 2 1 2];
%��� NEH ����� ���������������
%[~,NEH1]=NawazHehCmax(p1',m1,ma1,S_NEH);
[~,NEH]=NawazHehCmax2(p1',m1,ma1,S_NEH,r);

S_NEH = [];
%S_NEH = zeros(1,2*n);
parfor ii=1:n %parfor
    %S_NEH(2*ii-1) = NEH(ii);
    %S_NEH(2*ii) = NEH(ii);
    tmp = [NEH(ii) NEH(ii)];
    S_NEH = [S_NEH tmp]; %���������
end

%S = JeongFullS(S_NEH);
S = JeongFullSp(S_NEH);
%[T_NEH_1212,~,~,~]= XieTargFun7(p1, S, n, m1, ma1, d1, su1); % ��� 1-2-1-2 ���������, ��� ������ ���������
% �����! ������� r
[T_NEH_1212,~,~,~]= XieTargFun7r(p1, S, n, m1, ma1, d1, su1, r); % ��� 1-2-1-2 ���������, ��� ������ ���������

xmin=S_NEH;
record = T_NEH_1212; %��� ��� ������������

%�������� �������� ���������������� FL - MFL ��� ���������� �����������
%MEDD
d_m = d1 - (p1(3,:) + p1(4,:));
[~, I] = sort(d_m);

S_MEDD = [];
%S_MEDD = zeros(1,2*n);
parfor ii=1:n %parfor
    %S_MEDD(2*ii-1) = I(ii);
    %S_MEDD(2*ii) = I(ii);
    tmp = [I(ii) I(ii)];
    S_MEDD = [S_MEDD tmp];
end

%S_new = [5 4 5 9 3 8 9 7 6 2 7 1 2 8 1 6 10 3 10 4];
%S_new_full = JeongFullS(S_new);
%[T_S_new,C1_S_new,C2_S_new, end_time_ex_S_new]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1);

S0 = S_MEDD;
S=[];
%S_new_par = [];
ri=1;
%� ������ ���� �� ��������������
%while ri<=2*n
for  ri = 1:2*n   
    sr = S0(ri);
    T_min = inf;
    S_min = [];
    C_min = inf;
    for l=1:numel(S) +1
        if l==1
            S_new_par = [sr S];
        elseif l== numel(S) +1
            S_new_par = [S sr];
        else
            S_new_par = [S(1:l-1) sr S(l:numel(S))];
        end
        %����� ������ �� ��������� ������� �������������, �� �������� � ��� �������
        S_new = [S_new_par S0(ri+1:numel(S0))];
        S_new_full = JeongFullS(S_new);
         %�����! ������ d - d1
        %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % ��� 1-2-1-2 ���������, ��� ������ ��������� %(p1, S, n, m1, ma1, d, su1);
        % �����! ������� r
        [T_S_new,~,C2_S_new, ~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % ��� 1-2-1-2 ���������, ��� ������ ��������� %(p1, S, n, m1, ma1, d, su1);
        if T_S_new<T_min  || (T_S_new==T_min && C2_S_new<C_min)
            T_min = T_S_new;
            S_min = S_new_par;
            C_min = C2_S_new;
        end
    end
    S=S_min;
    
    %����� MNEH ������ FL
    if ri>=3
        %T_min = inf;
        S_min = [];
        %C_min = inf;
        %clear S_new_par; 
        T_S_new=ones(numel(S)-1,numel(S))*inf;
        %clear T_S_new;
        %cc=1;
        parfor ii=1:numel(S)-1
            %clear T_S_new;
            %T_S_new(ii,:) = ones(1,2);
           T_S_new(ii,:) = S_xchange(ii,S,ri,S0, p1, n, m1, ma1, d1, su1, r);
        end
        [T_S_new1,ind1] = min(T_S_new);
        [T_min,ind2] = min(T_S_new1);
        i=ind1(ind2);
        j=ind2;
        if T_min<T_min_old
            %S=S_min;
            tempo = S(j);
            S(j)=S(i);
            S(i) = tempo;
        end    
    else    
    end
    
    %ri=ri+1;
end

if T_min<record
    xmin=S;
    record = T_min;
end

U=xmin; %
S_full=[];
end_time_ex = [];
%k2=1;
%k3=1;
%LocDom2=zeros(numel(U),2);
%LocDom3=zeros(numel(U),3);
% LocDom2 = [];
% LocDom3 =[];

Un=unique(U,'stable');
%stt = tic;%
%[LocDom2, LocDom3] = KUMZ_get_dom(S_full,Un, p,n,m,ma,d,su,end_time_ex);
%[LocDom2, LocDom3] = KUMZ_get_domP(S_full,Un, p,n,m,ma,d,su,end_time_ex);
[LocDom2, LocDom3] = KUMZ_get_domPsu(S_full,Un, p,n,m,ma,d,su,end_time_ex, r);
%endt = toc(stt);

LocDom2t = unique(LocDom2,'rows');
LocDom3t2 = unique(LocDom3,'rows');
%LocDom3t2 = [];

%U=1:n; %����� ������� �������
%U=[U U]; %����� ������� �������
S=[];
%S_full=[];
%end_time_ex = [];

%������� ������������� 0 ������ �� �� ������
LocDom2 = [];
LocDom3 = [];
LocDom3t1 = [];

ccount = 0;

Cyc=2;

%record1=record;

if record>0 && flag_branch ==1 %1111111
    %������� ��� ??? ���� ���
    %matlabpool close
    if ~isempty(gcp('nocreate')) %���
        delete(gcp('nocreate'));
    end    
    
    stt = tic;%
    %[~] = JeongBranchAndBound9(p1,su1,r,d1,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc);
    in_arg = {p1,su1,r,d1,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc,record,xmin,ccount};%,record,xmin
    %���������������� �������
    %[~,xmin] = JeongBranchAndBound13(in_arg); % %�����! � ������ ������� ��������� ������ �������� su
    
    %������������ �������
    [~,xmin] = JeongBranchAndBound11p(in_arg); % %�����! � ������ ������� ��������� ������ �������� su
    
    %����
    %[T_S_test,~,C2_S_new_test, end_t_test]= XieTargFun7r(p1, xmin, n, m1, ma1, d1, su1, r); 
    
    endt = toc(stt);
    S121_1 = xmin;
else
    %S121_1 = [2 8 7 6 7 6 10 9 3 10 1 3 8 1 5 4 5 4 9 2];
    S121_1 = xmin;
end


%--------------------------------------------------------------------------
% ������� ������ ��� ������: (1-2-1-2)-2
% �. 1. ������ ������ (1-3-1)=1,(2,3,4)=2,(1-3-1)-2=1,(2,3,4)-2=2,
% ma=[1 3 1 2 3 4 1 3 1 2 3 4 1 3 1 2 3 4 5];
p1=zeros(4,n);
p1(1,:) = sum(p(7:9,:));
p1(2,:) = sum(p(10:12,:));
p1(3,:) = sum(p(13:15,:));
p1(4,:) = sum(p(16:19,:));
%
su1 = cell(1,n);
%for ii=1:4 %parfor
    su1{1}(:,:)=su{7} + su{8} + su{9}; %sum(su{ii}(1:3,:));
    su1{2}(:,:)=su{10} + su{11} + su{12}; %sum(su{ii}(4:6,:));
    su1{3}(:,:)=su{13} + su{14} + su{15}; %sum(su{ii}(7:9,:));
    su1{4}(:,:)=su{16} + su{17} + su{18} + su{19}; %sum(su{ii}(10:12,:));
%end

%d1 = d - sum(p(1:6,:));
%�.�. d ��� � ������ �������� �������������
d1 = d - p(1,:) - p(2,:)*473 - p(3,:) - p(4,:)*18 - p(5,:)*473 - p(6,:); 
% ����� ��� 2� ������ ����������� �� r �� �������
r1 = zeros(1,n);

%�������� �������� NEH ��� ���������� �����������
S_NEH=1:n;
m1=2;
ma1=[1 2 1 2];
%��� NEH ����� ���������������
%[F_NEH,NEH]=NawazHehCmax(p1',m1,ma1,S_NEH); %������ ��� r, �.�. ������ ����� ���������� �� r �� �������
[F_NEH,NEH]=NawazHehCmax2(p1',m1,ma1,S_NEH,r1); %

S_NEH = [];
%S_NEH = zeros(1,2*n);
parfor ii=1:n %parfor
    %S_NEH(2*ii-1) = NEH(ii);
    %S_NEH(2*ii) = NEH(ii);
    tmp = [NEH(ii) NEH(ii)];
    S_NEH = [S_NEH tmp]; %���������
end
%S = JeongFullS(S_NEH);
S = JeongFullSp(S_NEH);
%[T_NEH_1212,~,~,~]= XieTargFun7(p1, S, n, m1, ma1, d1, su1); %������ ��� r, �.�. ������ ����� ���������� �� r �� �������
[T_NEH_1212,~,~,~]= XieTargFun7r(p1, S, n, m1, ma1, d1, su1, r1); %

%%if T_NEH<record
xmin=S_NEH;
record = T_NEH_1212;
%%end

%�������� �������� ���������������� FL - MFL ��� ���������� �����������
%MEDD
d_m = d1 - (p1(3,:) + p1(4,:));
[d_sort, I] = sort(d_m);

S_MEDD = [];
%S_MEDD = zeros(1,2*n);
parfor ii=1:n %parfor
    %S_MEDD(2*ii-1) = I(ii);
    %S_MEDD(2*ii) = I(ii);
    tmp = [I(ii) I(ii)];
    S_MEDD = [S_MEDD tmp];
end

%S_new = [5 4 5 9 3 8 9 7 6 2 7 1 2 8 1 6 10 3 10 4];
%S_new_full = JeongFullS(S_new);
%[T_S_new,C1_S_new,C2_S_new, end_time_ex_S_new]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1);

S0 = S_MEDD;
S=[];
ri=1;
while ri<=2*n
    sr = S0(ri);
    T_min = inf;
    S_min = [];
    C_min = inf;
    for l=1:numel(S) +1
        if l==1
            S_new_par = [sr S];
        elseif l== numel(S) +1
            S_new_par = [S sr];
        else
            S_new_par = [S(1:l-1) sr S(l:numel(S))];
        end
        %����� ������ �� ��������� ������� �������������, �� �������� � ��� �������
        S_new = [S_new_par S0(ri+1:numel(S0))];
        S_new_full = JeongFullS(S_new);
        %S_new_full = JeongFullSp(S_new);
        %�����! ������ d - d1
        %[T_S_new,C1_S_new,C2_S_new, end_time_ex_S_new]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); %(p1, S, n, m1, ma1, d, su1);
        % �����! ������� r - ��� ������ ������ r1 �������
        [T_S_new,C1_S_new,C2_S_new, end_time_ex_S_new]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r1); %(p1, S, n, m1, ma1, d, su1);
        if T_S_new<T_min || (T_S_new==T_min && C2_S_new<C_min)
            T_min = T_S_new;
            S_min = S_new_par;
            C_min = C2_S_new;
        end
    end
    S=S_min;
    %����� MNEH ������ FL
    if ri>=3
        T_min = inf;
        S_min = [];
        C_min = inf;
        for ii=1:numel(S)-1
            for jj=1:numel(S)
                S_new_par = S;
                tmp = S_new_par(ii);
                S_new_par(ii) = S_new_par(jj);
                S_new_par(jj) = tmp;
                S_new = [S_new_par S0(ri+1:numel(S0))];
                S_new_full = JeongFullS(S_new);
                %S_new_full = JeongFullSp(S_new);
                %�����! ������ d - d1
                %[T_S_new,C1_S_new,C2_S_new, end_time_ex_S_new]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); %(p1, S, n, m1, ma1, d, su1);
                % �����! ������� r - ��� ������ ������ r1 �������
                [T_S_new,C1_S_new,C2_S_new, end_time_ex_S_new]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r1); %(p1, S, n, m1, ma1, d, su1);
                if T_S_new<T_min || (T_S_new==T_min && C2_S_new<C_min)
                    T_min = T_S_new;
                    S_min = S_new_par;
                    C_min = C2_S_new;
                end
            end
        end
    end
    S=S_min;
    ri=ri+1;
end

if T_min<record
    xmin=S;
    record = T_min;
end

U=xmin; %
S_full=[];
end_time_ex = [];
%k2=1;
%$k3=1;
%LocDom2=zeros(numel(U),2);
%LocDom3=zeros(numel(U),3);

Un=unique(U,'stable');
%[LocDom2, LocDom3] = KUMZ_get_dom(S_full,Un, p,n,m,ma,d,su,end_time_ex);
%[LocDom2, LocDom3] = KUMZ_get_domP(S_full,Un, p,n,m,ma,d,su,end_time_ex);

%[LocDom2, LocDom3] = KUMZ_get_domPsu(S_full,Un, p,n,m,ma,d,su,end_time_ex);

LocDom2t = unique(LocDom2,'rows');
LocDom3t2 = unique(LocDom3,'rows');

%U=1:n; %����� ������� �������
%U=[U U]; %����� ������� �������
S=[];
%S_full=[];
%end_time_ex = [];

%������� ������������� 0 ������ �� �� ������
LocDom2 = [];
LocDom3 = [];
LocDom3t1 = [];

Cyc=2;

flag_branch=0;
if record>0 && flag_branch ==1 %22222222222
    %������� ����
    if ~isempty(gcp('nocreate')) %���
        delete(gcp('nocreate'));
    end    
    
    stt1 = tic;%
    %[~] = JeongBranchAndBound9(p1,su1,r,d1,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc);
    in_arg = {p1,su1,r1,d1,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc,record,xmin,ccount};%,record,xmin
    %���������������� �������
    %[~,xmin] = JeongBranchAndBound13(in_arg); % %�����! � ������ ������� ��������� ������ �������� su
    
    %������������ �������
    [record,xmin] = JeongBranchAndBound11p(in_arg); %�����! � ������ ������� ��������� ������ �������� su
    
    endt1 = toc(stt1);
    S121_2 = xmin;
else
    %S121_2 = [7 6 7 6 10 9 10 1 9 3 2 8 3 2 8 1 5 4 5 4];
    S121_2 = xmin;
end


%--------------------����� ������� ���������� S121_1-----------------------
%������� ������� �������
S_1212f1 = KumzFullS12(S121_1); % (perm) ���������� - ������-� ������ �������:
%S = KumzFullS12_2(S121_1); % (perm) ���������� - ������-� �� ��
% �����! XieTargFun6 �� ������������
%[T_1212f1,~,~, end_time_ex_1212f1]= XieTargFun7(p, S, n, m, ma, d, su);
%[T_1212f1,~,~, end_time_ex_1212f1]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_1212f1,~,~, end_time_ex_1212f1]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_1212f1,~,~, end_time_ex_1212f1]= KUMZTargFun19f(p, S_1212f1, n, m, ma, d, su, mq, r, mq_u);
FT_1212f1 = (T_1212f1 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_1212f1)) - max(max(end_time_ex_EDD)));

%��� ��������� ��������� ��� 1-2 ������� ������������� �� a+b ������
ab1 = sum(p(13:19,:));
[~,S_ab] = sort(ab1);
S_1212_ab=KumzFullS6_1(S121_1,S_ab);
%[T_1212_ab,~,~, end_time_ex_1212_ab]= XieTargFun7(p, S, n, m, ma, d, su);
%[T_1212_ab,~,~, end_time_ex_1212_ab]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_1212_ab,~,~, end_time_ex_1212_ab]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_1212_ab,~,~, end_time_ex_1212_ab]= KUMZTargFun19f(p, S_1212_ab, n, m, ma, d, su, mq, r, mq_u);
FT_1212_ab = (T_1212_ab - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_1212_ab)) - max(max(end_time_ex_EDD)));

%--------------------����� ������� ���������� S121_2-----------------------
%������� ������� �������
S_1212f2 = KumzFullS12(S121_2); % (perm) ���������� - ������-� ������ �������:
%S = KumzFullS12_2(S121_2); % (perm) ���������� - ������-� �� ��
% �����! XieTargFun6 �� ������������
%[T_1212f2,~,~, end_time_ex_1212f2] = XieTargFun7(p, S, n, m, ma, d, su);
%[T_1212f2,~,~, end_time_ex_1212f2] = XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_1212f2,~,~, end_time_ex_1212f2] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
%[T_1212f2,~,~, end_time_ex_1212f2] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
%[T_1212f2,~,~, end_time_ex_1212f2] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_1212f2,~,~, end_time_ex_1212f2] = KUMZTargFun19f(p, S_1212f2, n, m, ma, d, su, mq, r, mq_u);
FT_1212f2 = (T_1212f2 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_1212f2)) - max(max(end_time_ex_EDD)));

S_1212=KumzFullS6(S121_1,S121_2);
%[T_1212,~,~, end_time_ex_1212]= XieTargFun7(p, S, n, m, ma, d, su);
%[T_1212,~,~, end_time_ex_1212]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_1212,~,~, end_time_ex_1212]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_1212,~,~, end_time_ex_1212]= KUMZTargFun19f(p, S_1212, n, m, ma, d, su, mq, r, mq_u);
FT_1212 = (T_1212 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_1212)) - max(max(end_time_ex_EDD)));

% ��� ��������� ��������� ��� 1-2 ������� ������������� �� a+b ������
ab1 = sum(p(1:6,:));
[~,S_ab] = sort(ab1);
S_ab_1212=KumzFullS6(S_ab,S121_2);
%[T_ab_1212,~,~, end_time_ex_ab_1212]= XieTargFun7(p, S, n, m, ma, d, su);
%[T_ab_1212,~,~, end_time_ex_ab_1212]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_ab_1212,~,~, end_time_ex_ab_1212]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_ab_1212,~,~, end_time_ex_ab_1212]= KUMZTargFun19f(p, S_ab_1212, n, m, ma, d, su, mq, r, mq_u);
FT_ab_1212 = (T_ab_1212 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_ab_1212)) - max(max(end_time_ex_EDD)));

% ��� ��������� ��������� ������ 1-2 ������ �� �������� �������������
a1=sum(p(1:3,:));
b1=sum(p(4:6,:));
p1=[a1;b1]';
S_john= jhons(p1);
S_joh_1212=KumzFullS6(S_john,S121_2);
%[T_joh_1212,~,~, end_time_ex_joh_1212]= XieTargFun7(p, S, n, m, ma, d, su);
%[T_joh_1212,~,~, end_time_ex_joh_1212]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_joh_1212,~,~, end_time_ex_joh_1212]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_joh_1212,~,~, end_time_ex_joh_1212]= KUMZTargFun19f(p, S_joh_1212, n, m, ma, d, su, mq, r, mq_u);
FT_joh_1212 = (T_joh_1212 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_joh_1212)) - max(max(end_time_ex_EDD)));


%-------------����� ������� ���������� S121_1 - x - S121_2-----------------
% �1. ������� ��������������� 2n - �� ������ ������ ��������
%[max_u_SACO,S_opt_SACO]=SACO_KUMZ_12_x_12(p, n, m, ma, d, su,S121_1,S121_2);  %��� SACO ����� ������ r !!!
% �2. ��������� ��������������� n(���������)
%[max_u_SACO,S_opt_SACO,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_3(p, n, m, ma, d, su,S121_1,S121_2);  %
% stt3 = tic;
%�� ������������ ������
% [max_u_SACO,S_opt_SACO,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_4(p, n, m, ma, d, su,S121_1,S121_2, r);  %
% endt3 = toc(stt3);
% stt3 = tic;
%������������ ������
[max_u_SACO,S_opt_SACO,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_4p(p, n, m, ma, d, su,S121_1,S121_2, r);  %
%endt3 = toc(stt3);
%[max_u_SACO,S_opt_SACO,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_5(p, n, m, ma, d, su,S121_1,S121_2,mq);  %
S_12_SACO_12=KumzFullS11(S1_s,S_opt_SACO,S2_e,S2_U);
%[T_12_SACO_12,~,~, end_time_12_SACO_12]= XieTargFun7(p, S, n, m, ma, d, su);
%[T_12_SACO_12,~,~, end_time_12_SACO_12]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
%[T_12_SACO_12,~,~, end_time_12_SACO_12]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
[T_12_SACO_12,~,~, end_time_12_SACO_12]= KUMZTargFun19f(p, S_12_SACO_12, n, m, ma, d, su, mq, r, mq_u);
FT_12_SACO_12 = (T_12_SACO_12 - T_EDD)*k1 + (1-k1)*(max(max(end_time_12_SACO_12)) - max(max(end_time_ex_EDD)));

%----------------------------------------------------------------------------------------------------
%����� ������� ���������� � ������ �������
FF = [Fmaxp FNEH_f FT_1212_ab FT_1212f1 FT_1212f2 FT_12_SACO_12 FT_ab_1212 FT_joh_1212];
[Fmin,ind] = min(FF);
SS = [S_maxp; S_NEH_f; S_1212_ab; S_1212f1; S_1212f2; S_12_SACO_12; S_ab_1212; S_joh_1212];
Smin = SS(ind,:);

%[end_time_ex_12_SACO_12_full]= KUMZ_data_load(p,S, n, m, ma, d, su, mq, db_data,id_in);
%������� � ������ r
%[end_time_ex_12_SACO_12_full]= KUMZ_data_load2(p,S, n, m, ma, d, su, mq, db_data,id_in,r);
[start_time_12_S_12_full, end_time_12_S_12_full]= KUMZ_data_load3(p,Smin, n, m, ma, d, su, mq, r, mq_u, db_data, id_t, id_in);
 
session_id = 3;
data_out = KUMZ_data_out2( start_time_12_S_12_full, end_time_12_S_12_full, db_data, session_id, id_t, id_in);

delete(gcp('nocreate'));

%-------��� �����������------------------------------
% %
% % %�.2. ��������� ������� ����� (1,3,1)-1 � (1,3,1)-2
% % % ������������ d
% [~,I]=sort(xmin1);
% d1=d(xmin1);
% H=100;
%
% for ii=1:n
%     d1(ii)=d1(ii)-H/ii+H/n;
% end
% d2=d1(I);
%
% % ������������ d �� ������ ������
% [~,I]=sort(xmin2);
% d1=d2(xmin1);
% %H=200;
%
% for ii=1:n
%     d1(ii)=d1(ii)-H/ii+H/n;
% end
% d2=d1(I);
%
%%F = JeongBranchAndBound8(p1,su1,d2,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc);
%�������1 [4 5 2 3 3 2 5 4 1 1] ���������� ���� ����� �� ������ �����
%S121 = [4 5 2 3 3 2 5 4 1 1];
%S121 = [4 3 3 2 5 2 5 4 1 1];
%S121 = xmin;
%S = KumzFullS4(S121); % (perm) ���������� - ������-� �� ������ ��: ��
%����� � ������ ����������!!!
%S = KumzFullS5(S121); % (perm) ���������� - ������-� ������ �������:
%[T_6,~,~, end_time_ex_6]= XieTargFun8(p, S, n, m, ma, d, su, mq);

%�������2 [4 5 4 2 3 3 2 5 1 1] ���������� ���� ����� �� ������ �����
% S121 = [4 5 4 2 3 3 2 5 1 1];
% %S121 = xmin;
% S = KumzFullS5(S121); % (perm) ���������� - ������-� ������ �������:
% [T_7,~,~, end_time_ex_7]= XieTargFun8(p, S, n, m, ma, d, su, mq);

%������ ������ ������� ���������� 1-3-1-2-3-4-1-3-1-2-3-4
%[max_u_SACO,S_opt_SACO]=SACO_KUMZ(p, n, m, ma, d, su,S);  %��� SACO ����� ������ r !!!
% S_opt_SACO = [1 3 3 1 3 4 3 2 1 5 4 3 2 1 4 5 5 1 4 2 3 1 1 1 1 5 4 5 3 4 1 5 2 1 1 4 3 4 5 5 3 2 5 5 4 3 2 3 4 3 2 4 2 2 4 2 5 2 5 2];
%[T_8,~,~, end_time_ex_8]= XieTargFun7(p, S_opt_SACO, n, m, ma, d, su);

% ������ � ���������� Xie
% � ���� ������ ������ ������ 1,3,1,(2,3,4),1 ������
% ��� ������ ������� ���������
%  2-� ����� ���������� �������� �� Wang 1-3-1
%  ��� 1-2 ������� ������������� �� a+b ������
% c1 = p(3,:);
% d1 = sum(p(4:6,:));
% e1 = p(7,:);
%
% c = sum(c1);
% d = sum(d1);
% e = sum(e1);
%
% pj=[c1;d1]';
% Sj= jhons(pj);
% F=WangBranchAndBoundFirst(c1,d1,e1,Sj,n, c, d, e);
%
% ab1 = a1 + b1;
% [ab1_sort,S_ab] = sort(ab1);

%end