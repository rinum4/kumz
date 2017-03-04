function [max_u,xmin,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_4(p, n, m, ma, d, su,S1_const,S2_const, r)
%Ant Colony Optimization for
%ma=[(1 3 1)-S1 (2 3 4)-S1 (1 3 1)-x (2 3 4)-x (1 3 1)-S2 (2 3 4 5)-S2];
%�����! � ������ ������� ��������� ������ �������� su

%������������ ���������� ACO
%nVar = n*numel(ma);
nVar=n;

%�������� � S2 ������� ���������� ����������
[~,nc] = size(S2_const);
%S=[];
S2_U=[];
S2_e=[];
%�������� - x ������ ������
%for i=nr:-1:1
    for j=nc:-1:(nc-n+1) %parfor
        ind = find(S2_U==S2_const(1,j), 1);
        if isempty(ind)
            cur_job = ones(1,7)*S2_const(1,j);
        else
            cur_job = ones(1,6)*S2_const(1,j);
        end
        S2_U = [S2_U S2_const(1,j)];
        S2_e=[cur_job S2_e];
    end
%end

%S1 - ���� ������ � S2_U �������
%S1_un=unique(S1,'stable');
[~,nc] = size(S1_const);
S1_s=[];
S1_U=[];
count=0;
%�������� - S ������ ������
%for i=1:nr
    for j=1:nc %parfor �� ������
        if count<n
            ind = find(S2_U==S1_const(1,j));
            if numel(ind)==2
                ind2 = find(S1_U==S1_const(1,j),1);
                if isempty(ind2)
                    cur_job = ones(1,6)*S1_const(1,j);
                    S1_U = [S1_U S1_const(1,j)];
                    S1_s=[S1_s cur_job];
                    count = count+1;
                else
                    cur_job = [];
                end
            else
                cur_job = ones(1,6)*S1_const(1,j);
                S1_U = [S1_U S1_const(1,j)];
                S1_s=[S1_s cur_job];
                count = count+1;
            end
        else
            break;
        end    
    end
%end

T1 = unique(S2_const,'stable'); % ������ ������� ���������� ������������
%T1 = [T1 T1]; % 

[~,nc] = size(T1);
T=[];
%for i=1:nr
    for j=1:nc %parfor 
            ind1 = find(S1_U==T1(1,j));
            ind2 = find(S2_U==T1(1,j));
            if (numel(ind2)==2 && numel(ind1)==1) || (numel(ind2)==1 && numel(ind1)==2)
                cur_job = [];
            elseif (numel(ind2)==2 && numel(ind1)==0) || (numel(ind2)==1 && numel(ind1)==1) || (numel(ind2)==0 && numel(ind1)==2)
                cur_job = T1(1,j);
            elseif (numel(ind2)==1 && numel(ind1)==0) || (numel(ind2)==0 && numel(ind1)==1)
                cur_job = [T1(1,j) T1(1,j)];
            else
                cur_job = [T1(1,j) T1(1,j) T1(1,j)];
            end
            T=[T cur_job];
    end
%end

N = 300; %������������ ���������� �������� - ������������ ��������
N_max = 10; %������������ ���������� �������� ��� ���������
M = 10;   %���������� ��������

Q = 1;          %�������� �������� 1
%Q2 = 100000;
rho=0.9;       %�������� �������� 2
evap = 0.05;%�������� �������� 3
tau_max = 1*Q/(1-rho); %��������� ����������� ����������
tau_min = tau_max/5;

tau=tau_max*ones(nVar,nVar);   % ������� ����������
%S_opt = T;
max_u=inf;
max_u_old=inf;
%max_u_i=zeros(1,N);
iteration=1;
N_un = 0;
fl_st=0;
%for iteration=1:N % ������������ ��������
while N_un < N_max && iteration<N % �������� ���� ������� �����������������
    %max_u_i(iteration) = 0;
    p0=log(iteration)/log(N);
    %max_u_old=max_u;
    for Ant = 1:M
        %L=T; �� L � ���������� ������ �� �������
        S(Ant,:)=zeros(1,nVar);
        %S1=[];
        %Cost(Ant)=0;
        if Ant>1 && fl_st==0
            tau_i=tau;
            position = 1;
            %max_ant=inf;
            %k_old=0;
            while position <= nVar
                pr = rand();
                if pr>=p0
                    %��������������
                    %�������� ��������
                    P=tau_i(:,position);
                    P=P/sum(P); %�������� ����������� � �������� [0 1]
                    
                    r=rand;
                    C=cumsum(P);
                    k=find(r<=C,1,'first');
                else
                    %��������������
                    %�������� �� ������������ ������������
                    [~,k]=max(tau_i(:,position));
                end
                
                if isempty(k)
                    %������ �������, �������� �������
                    position = 1;
                    tau_i=tau;
                    %max_ant=inf;
                    continue;
                end
                
                %S1(position)=k;
                
                S(Ant,position)=k;
                tau_i(k,:) = 0;
                %max_ant=LB;
                %L(k)=0; % �� L � ���������� ������ �� �������
                position = position + 1;
            end %pos
        else
            S(Ant,:)=1:nVar;
            fl_st=1;
        end
        %API - Adjacent Pairwise Interchange - �������� ������������
        %��������� ����� - ������������ 1-API,2-API,...,(n-j)-API
        
        %�.�. � ������� ���������� �� 1 �� nVar ����� ���������� ���
        %��������� � ���� ���
        max_S = S(Ant,:);
        full_S = KumzFullS11(S1_s,T(S(Ant,:)),S2_e,S2_U);
        max_loc_s = XieTargFun7r(p, full_S, n, m, ma, d, su, r);
        for j=1:nVar-1
            for i=j+1:nVar
                cur_S = S(Ant,:);
                tempo = cur_S(j);
                cur_S(j)=cur_S(i);
                cur_S(i) = tempo;
                full_S = KumzFullS11(S1_s,T(cur_S),S2_e,S2_U);
                Cost_i_j = XieTargFun7r(p,full_S, n, m, ma, d, su, r);
                if Cost_i_j<max_loc_s
                    max_S = cur_S;
                    max_loc_s = Cost_i_j;
                end
            end
        end
        %��������� � ������� ������ ��������� ��������������� ����������
        S(Ant,:) =  max_S;
        %Cost(Ant) = CroceTargFun(p,S(Ant,:));
        Cost(Ant) = max_loc_s;
        
        if Cost(Ant)<max_u
            max_u=Cost(Ant);
            %max_u_i(iteration)=Cost(Ant); %��� �������
            xmin = T(S(Ant,:));
            %N_un = 0; % ������� ���������
        %else
        %    N_un = N_un+1; % ������� �� ���������
        end
    end %Ant
    if max_u<max_u_old
        N_un = 0; % ������� ���������
        max_u_old = max_u;
    else
        N_un = N_un+1; % ������� �� ���������
    end
    %��������� ������������ ����������
    for k=1:M
        tour=S(k,:);
        tour=[tour tour(1)]; %#ok
        for l=1:nVar
            i=tour(l);
            j=tour(l+1);
            tau(i,j)=tau(i,j)+Q/Cost(k); %Q2
        end
    end
    % ���������
    tau=(1-evap)*tau;
    
    tau=max(tau,tau_min); % ������������ �������� ����������
    
    %��� �������
    %if max_u_i(iteration)==0
    %    max_u_i(iteration) = max_u_i(iteration-1);
        %N_un = N_un+1; % ������� �� ���������
    %end
    iteration = iteration + 1;
end %iteration
%F=max_u;
%plot(max_u_i(:))
xmin=xmin;
end %fcnt




