function [S, Cost] = SACO_KUMZ_12_x_12_42pp(tau,S1_s,T2,S2_e,S2_U,nVar,p0,p, n, m, ma, d, su, r, flag_fl,flag_KUMZ, mq, mq_u, T_EDD, Cmax_EDD,k1)
%L=T; �� L � ���������� ������ �� �������
S=zeros(1,nVar);
%S1=[];
%Cost(Ant)=0;

%if Ant>1 && fl_st==0
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
            
            rr=rand;
            C=cumsum(P);
            k=find(rr<=C,1,'first');
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
        
        S(position)=k;
        tau_i(k,:) = 0;
        %max_ant=LB;
        %L(k)=0; % �� L � ���������� ������ �� �������
        position = position + 1;
    end %pos
%else %������� ��� �����������? - �� � ������������ ������
%   S=1:nVar;
%   fl_st=1;
%end

%API - Adjacent Pairwise Interchange - �������� ������������
%��������� ����� - ������������ 1-API,2-API,...,(n-j)-API

%�.�. � ������� ���������� �� 1 �� nVar ����� ���������� ���
%��������� � ���� ���
%max_S = S;
if ~isempty(S2_e)
    full_S = KumzFullS11(S1_s,T2(S),S2_e,S2_U);
else
    full_S = KumzFullS(T2(S),ma);
end    
if flag_KUMZ == 0
    %������� ������� ��� ����� ������������ ��
    [min_loc_s,~,Cmax2,~] = XieTargFun7r(p, full_S, n, m, ma, d, su, r);
    min_loc_s = (min_loc_s - T_EDD)*k1 + (1-k1)*(Cmax2 - Cmax_EDD);
else
    %������ �������
    [min_loc_s,~,Cmax2,~] = KUMZTargFun19f(p, full_S, n, m, ma, d, su, mq, r, mq_u); 
    min_loc_s = (min_loc_s - T_EDD)*k1 + (1-k1)*(Cmax2 - Cmax_EDD);
end    
if flag_fl == 1
    Costs=ones(nVar,nVar)*inf;
    parfor i=1:nVar-1 %parfor ����������
    %for i=1:1
        Costs(i,:) =  SACO_KUMZ_12_x_12_42ppp(i,S,S1_s,T2,S2_e,S2_U,nVar,p, n, m, ma, d, su, r,flag_KUMZ, mq, mq_u,T_EDD, Cmax_EDD,k1);
    end
    
    %��������� � ������� ������ ��������� ��������������� ����������
    [Costs1,ind1] = min(Costs);
    [Cost,ind2] = min(Costs1);
    i=ind1(ind2);
    j =ind2; 
    if Cost<min_loc_s
        tempo = S(j);
        S(j)=S(i);
        S(i) = tempo;
    else
        Cost = min_loc_s;
    end
else
    Cost = min_loc_s;
end

%if Cost < max_u %� ������������ ������ max_u �� ����� ������������ �����
%���������� !
%     max_u=Cost;
%     %max_u_i(iteration)=Cost(Ant); %��� �������
%     xmin = T(S);
%     %N_un = 0; % ������� ���������
%     %else
%     %    N_un = N_un+1; % ������� �� ���������
%end

end

