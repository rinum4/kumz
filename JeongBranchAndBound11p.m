function [record, xmin]=JeongBranchAndBound11p(in_arg) %
% F2, (1-2-1-2), reentr, sdst, T
% Jeong (2014)
% su - setup times
% d - due dates
%�����! � ������ ������� ��������� ������ �������� su

%global record;
%global xmin;
%global ccount; %������� ������� ��� �����

p = in_arg{1};
su  = in_arg{2};
r  = in_arg{3};
d = in_arg{4};
S = in_arg{5};
U = in_arg{6};
lvl = in_arg{7};
n = in_arg{8};
m = in_arg{9};
ma = in_arg{10};
end_time_ex = in_arg{11};
LocDom2_in = in_arg{12};
LocDom3_in = in_arg{13};
LocDom2t_in = in_arg{14};
LocDom3t1_in = in_arg{15};
LocDom3t2_in = in_arg{16};
Cyc = in_arg{17};
record = in_arg{18};
xmin = in_arg{19};
ccount = in_arg{20};

profileName = 'local';
clust = parcluster(profileName);

%F=0;

S1 = S;
lvl_1 = lvl;

% if numel(S1)==1
%          if S1==[4] 
%                 fl = 1;
%           end
% end

%i=1;
U2=unique(U, 'stable');
%while i<numel(U2)+1

%
max_cp = 3;%�����-�� ������������ ������������
max_iter = ceil(n/max_cp);
%flag = 0;

for cp=1:max_iter
    flag = 0;
    job = createJob(clust);
    
    start_i = (cp-1)*max_cp + 1;
    if cp==max_iter
        if rem(n,max_cp) ~=0
            stop_i = start_i + rem(n,max_cp) - 1;
        else
            stop_i = start_i + max_cp - 1;
        end    
    else    
        stop_i = start_i + max_cp - 1;
    end    
    for i=start_i:stop_i %parfor �� ����������
        if i>1
            ind = find(S==U2(i-1), 1);
            if isempty(ind)
                %���� �� ���� ����� �������� 1-2 ���
                end_time_ex(1:2,U2(i-1)) = [0;0];
            else
                %���� ��� - �������� 3-4 ���
                end_time_ex(3:4,U2(i-1)) = [0;0];
            end
        end
        
        %     if numel(S)==1
        %          if S==[4] && U2(i)==5
        %                 fl = 1;
        %           end
        %     end
        
        fl_dom=0;
        
        %������� �������� �� 1-� ������ � ���������� �� ������������ r
        %���������
        %         ind = find(S==U2(i), 1);
        %         if isempty(ind)
        %             if r(U2(i))>0
        %                 if max(end_time_ex(1,:))<r(U2(i))
        %                     fl_dom=1;
        %                 end
        %             end
        %         end
        
        %�������� ��������� ���������� �������������
        if numel(S)>=2
            test_ijl = [S(lvl - 2) S(lvl - 1) U2(i)];
            test_ij = [S(lvl - 1) U2(i)];
        elseif  numel(S)>=1
            test_ijl =[];
            test_ij = [S(lvl - 1) U2(i)];
        else
            test_ijl=[];
            test_ij = [];
        end
        
        [nrl2,~]=size(LocDom2_in); %2-�� ���������� �� ������� ����
        if numel(test_ij)>0
            if numel(LocDom2_in)>0
                for j=1:nrl2
                    if test_ij == LocDom2_in(j,:) %
                        fl_dom=1;
                    end
                end
            end
        end
        
        [nrl3,~]=size(LocDom3_in); %3-�� ���������� �� ��������� ����
        if numel(test_ijl)>0
            if numel(LocDom3_in)>0
                for j=1:nrl3
                    if test_ijl == LocDom3_in(j,:)
                        fl_dom=1;
                    end
                end
            end
        end
        
        if fl_dom==1
            %i=i+1;
            continue;
        end
        
        U1 = U;
        
        S1(lvl_1) = U2(i); %��������� ���������� (2� ������ ���������)
        %��� ����������� ����� ���������� �������� end_time_ex8 ���� ��������
        %     if numel(S1)==7
        %         if S1==[1 2 1 4 2 4 3];
        %             fl = 1;
        %         end
        %     end
        
        %� 1. S2 ��� XieTargFun6
        %S2_1 = JeongFullS(S1); %������ ���������� �� 2-� ������� ��� ������������� �������
        
        %� 2. S2 ��� XieTargFun9
        S2_2 = [U2(i) U2(i)];
        S_full = JeongFullS(S); %������ ���������� �� 2-� ������� ��� ������������� �������
        %U1(i) = [];
        ind = find(U1==U2(i),1);
        U1(ind)=[];
        
        %� 1.
        %[T,C1_S,C2_S,end_time_ex]= XieTargFun7(p, S2_1, n, m, ma, d, su); %
        %� 2.
        %[T,C1_S,C2_S,end_time_ex] = XieTargFun10(p, S2_2, n, m, ma, d, su,S_full, end_time_ex); % % ��� 1-2-1-2 ���������, ��� ������ ���������
        % �����! ������� r
        [T,C1_S,C2_S,end_time_ex] = XieTargFun10r(p, S2_2, n, m, ma, d, su,S_full, end_time_ex,r); % % ��� 1-2-1-2 ���������, ��� ������ ���������
        
        %���������� ��������� ���������� ������������� ��� ���� ������
        %� ������� �� ���������� ���������� � ������ �����������
        %��������� ��������� �� ������ � ������ �����
        k2=1;
        k3=1;
        %LocDom2 = [];
        LocDom2=zeros(numel(U1)+1,2);
        %LocDom3 = [];
        LocDom3=zeros(numel(U1)+1,3);
        
        S1_full = JeongFullS(S1);
        %������� 1: �� ����������
        Un=unique(U1,'stable');
        %[LocDom2, LocDom3] = KUMZ_get_domP(S1_full,Un, p,n,m,ma,d,su,end_time_ex); %����������
        for ii=1:numel(Un)
            for j=1:numel(Un) %ii+1 %������� 3: �� ���������� �� �� ���� ��������
                %[U(i)  U(j)]
                if JeongCheckPropAll2_o5(S1_full,Un,Un(ii),Un(j), p,n,m,ma,d,su,end_time_ex,r)==1
                    LocDom2(k2,:) = [Un(j) Un(ii)];
                    k2=k2+1;
                end
                for l=1:numel(Un) %j+1  %������� 3: �� ���������� �� �� ���� ��������
                    %[U(i)  U(j) ,U(l)]
                    if JeongCheckPropAll3_o5(S1_full,Un,Un(ii),Un(j),Un(l), p,n,m,ma,d,su, end_time_ex,C1_S,r)==1
                        LocDom3(k3,:) = [Un(j) Un(ii) Un(l)];
                        k3=k3+1;
                    end
                end
            end
        end
        %������� 2: �� ����
        %     Un=U1;
        %     for ii=1:numel(Un)
        %         for j=1:numel(Un)
        %             %[U(i)  U(j)]
        %             if JeongCheckPropAll2_o4(S1_full,Un,Un(ii),Un(j), p,n,m,ma,d,su,end_time_ex)==1
        %                 LocDom2(k2,:) = [Un(j) Un(ii)];
        %                 k2=k2+1;
        %             end
        %             for l=1:numel(Un)
        %                 %[U(i)  U(j) ,U(l)]
        %                 if JeongCheckPropAll3_o4(S1_full,Un,Un(ii),Un(j),U1(l), p,n,m,ma,d,su, end_time_ex,C1_S)==1
        %                     LocDom3(k3,:) = [Un(j) Un(ii) Un(l)];
        %                     k3=k3+1;
        %                 end
        %             end
        %         end
        %     end
        
        LocDom2 = unique(LocDom2,'rows');
        %LocDom2(1,:)=[];
        LocDom3 = unique(LocDom3,'rows');
        %LocDom3(1,:)=[];
        
        LB = JeongCheckLB2(S1,U1, p, d, su, C1_S,C2_S);    %
        
        if LB<record %��� min < %T+LB
            if lvl_1==Cyc*n
                if T< record
                    record = min(record,T);
                    xmin = S1;
                end
                
                %ccount=ccount+1;
                %return;
            else
                %[record, xmin]=JeongBranchAndBound9(p,su,r,d,S1,U1,lvl_1 + 1,n,m,ma,end_time_ex, LocDom2t_in, LocDom3t1_in, LocDom2, LocDom3t2_in, LocDom3,Cyc,record,xmin);
                in_arg = {p,su,r,d,S1,U1,lvl_1 + 1,n,m,ma,end_time_ex, LocDom2t_in, LocDom3t1_in, LocDom2, LocDom3t2_in, LocDom3,Cyc,record,xmin,ccount}; %record,xmin
                %[record, xmin]=JeongBranchAndBound9(in_arg);
                createTask(job, @JeongBranchAndBound12, 2, {{in_arg}});
                flag = 1;
            end
        else
            %F=record;
        end
        %i=i+1;
    end %for i
    if flag == 1
        submit(job);
        wait(job);
        y = fetchOutputs(job);
        [ny,~] = size(y);
        delete(job);
    end
    
    for i=1:ny %parfor
        if y{i,1} < record
            record = y{i,1};
            xmin = y{i,2};
            %ccount = y{i,3};
        end
    end
    
end %cp

%delete(job);

% record = record;
% xmin = xmin;

end
