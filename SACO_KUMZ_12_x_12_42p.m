function [S, Cost] = SACO_KUMZ_12_x_12_42p(tau,S1_s,T2,S2_e,S2_U,nVar,p0,p, n, m, ma, d, su, r)
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
max_S = S;
full_S = KumzFullS11(S1_s,T2(S),S2_e,S2_U);
max_loc_s = XieTargFun7r(p, full_S, n, m, ma, d, su, r);

for j=1:nVar-1
    for i=j+1:nVar
        cur_S = S;
        tempo = cur_S(j);
        cur_S(j)=cur_S(i);
        cur_S(i) = tempo;
        full_S = KumzFullS11(S1_s,T2(cur_S),S2_e,S2_U);
        Cost_i_j = XieTargFun7r(p,full_S, n, m, ma, d, su, r);
        if Cost_i_j<max_loc_s
            max_S = cur_S;
            max_loc_s = Cost_i_j;
        end
    end
end

%��������� � ������� ������ ��������� ��������������� ����������
S =  max_S;
%Cost(Ant) = CroceTargFun(p,S);
Cost = max_loc_s;

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

