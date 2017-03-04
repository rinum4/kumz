function [ min_u_glob, xmin, S1_s,S2_e,S2_U ] = HGA_KUMZ( p, n, m, ma, d, su,S1_const,S2_const, r, mq, mq_u, flag_KUMZ,T_EDD, Cmax_EDD,k1)
% HGA  alg - 
%Chen (2009)

N = 300; %������������ ���������� ��������� - ������������ ��������
N_max = N/2; %������������ ���������� �������� ��� ���������
M = ceil(n/2);   %������ ��������� ceil(n/2) 50

Kross_p = 0.8; %����������� ����������
Mut_p = 0.3; %����������� �������
Other_p = 0.01; %����������� ������ �� ��������

%����������� - ������������ �� ������  = ������������������ �����

%�������� � S2 ������� ���������� ����������
[~,nc] = size(S2_const);
%S=[];
%S2_U=[];
S2_U=zeros(1,n);
S2_e=[]; %��������� ������ �.�. ����������� �����������
%�������� - x ������ ������
%for i=nr:-1:1
cc=1;
if ~isempty(S2_const)
    for j=nc:-1:(nc-n+1) %parfor �� ������ �.�. S2_U
        ind = find(S2_U==S2_const(1,j), 1);
        if isempty(ind)
            cur_job = ones(1,7)*S2_const(1,j);
        else
            cur_job = ones(1,6)*S2_const(1,j);
        end
        %S2_U = [S2_U S2_const(1,j)];
        S2_U(cc) = S2_const(1,j);
        cc=cc+1;
        S2_e=[cur_job S2_e];
    end
end    
%end

%S1 - ���� ������ � S2_U �������
%S1_un=unique(S1,'stable');
[~,nc] = size(S1_const);
S1_s=[]; %��������� ������ �.�. ����������� �����������
%S1_U=[];
S1_U=zeros(1,n);
count=0;
cc=1;
%�������� - S ������ ������
%for i=1:nr
if ~isempty(S1_const)
    for j=1:nc %parfor �� ������ �.�. S1_U
        if count<n
            ind = find(S2_U==S1_const(1,j));
            if numel(ind)==2
                ind2 = find(S1_U==S1_const(1,j),1);
                if isempty(ind2)
                    cur_job = ones(1,6)*S1_const(1,j);
                    %S1_U = [S1_U S1_const(1,j)];
                    S1_U(cc) = S1_const(1,j);
                    cc=cc+1;
                    S1_s=[S1_s cur_job];
                    count = count+1;
                else
                    %cur_job = [];
                end
            else
                cur_job = ones(1,6)*S1_const(1,j);
                %S1_U = [S1_U S1_const(1,j)];
                S1_U(cc) = S1_const(1,j);
                cc=cc+1;
                S1_s=[S1_s cur_job];
                count = count+1;
            end
        else
            break;
        end    
    end
end    
%end

if ~isempty(S2_const)
T1 = unique(S2_const,'stable'); % ������ ������� ���������� ������������
%T1 = [T1 T1]; % 

[~,nc] = size(T1);
%T2=[];
T2 = zeros(1,nc);
%for i=1:nr
ne=1;
    for j=1:nc %parfor ������������
            ind1 = find(S1_U==T1(1,j));
            ind2 = find(S2_U==T1(1,j));
            if (numel(ind2)==2 && numel(ind1)==1) || (numel(ind2)==1 && numel(ind1)==2)
                %cur_job = [];
            elseif (numel(ind2)==2 && numel(ind1)==0) || (numel(ind2)==1 && numel(ind1)==1) || (numel(ind2)==0 && numel(ind1)==2)
                cur_job = T1(1,j);
                T2(ne) = cur_job;
                ne=ne+1;
            elseif (numel(ind2)==1 && numel(ind1)==0) || (numel(ind2)==0 && numel(ind1)==1)
                cur_job = [T1(1,j) T1(1,j)];
                T2(ne:ne+1) = cur_job;
                ne=ne+2;
            else
                cur_job = [T1(1,j) T1(1,j) T1(1,j)];
                T2(ne:ne+2) = cur_job;
                ne=ne+3;
            end
            %T2=[T2 cur_job];
            %T2() = cur_job;
    end
else
    %T2=1:n;
    [~, I] = sort(d);
    T2 = I;
end    
%end

%������ � �2 ����������� ������

Par = zeros(M,n); %��������
%Un = zeros(1,n);
%�������� ��������� ��������� - �����������, �.�. � ��������� ������ ������
%� ��������� �������
for Chrom = 1:M % parfor �� ����������
    position = 1;
    Un = zeros(1,n);
    while position <= n
        k = ceil(rand()*n);
        ind = find(Un==k, 1);
        if ~isempty(ind)
            continue;
        end    
        Par(Chrom,position) = k;
        Un(position)=k;
        position = position + 1;
    end %pos
end    

iteration=1;
N_un = 0;
Off = zeros(M,n);
Cost = zeros(1,M);
min_u_glob = inf;
while N_un < N_max && iteration<N % �������� ���� ������� �����������������
    
    parfor Chrom = 1:M-1 % parfor ����������
        %������� ��������
        [Off(Chrom,:), Cost(Chrom)] = HGA_KUMZ_p(Par(Chrom,:),Par(Chrom+1,:),Kross_p,Mut_p,Other_p,S1_s,T2,S2_e,S2_U,p, n, m, ma, d, su, r, mq, mq_u, flag_KUMZ,T_EDD, Cmax_EDD,k1);
    end
    [Off(M,:), Cost(M)] = HGA_KUMZ_p(Par(M,:),Par(1,:),Kross_p,Mut_p,Other_p,S1_s,T2,S2_e,S2_U,p, n, m, ma, d, su, r, mq, mq_u,flag_KUMZ,T_EDD, Cmax_EDD,k1);
    
    [min_u,ind] = min(Cost(:));
    S_opt = Off(ind,:); %����������� ���������� ����� ��������
    
    %������ ������-�������
    %Gilles (1985) ��������� ������� ��������������
    F=zeros(1,M);
    alf = 1.005; %������� �� ������
    for i = 1:M % parfor ����������
        F(i) = (Cost(i) - min_u)^alf;
    end
    
    %��������� �������� - ������������ �������� ������� ���������
    %Goldberg (1985)
    %����� � ���������� ���� ������
    if sum(F)>0
        P=F/sum(F);
    else
        P=ones(1,M)/M;
    end    
    C=cumsum(P);
    for i=1:M
        rr=rand;
        k=find(rr<=C,1,'first');
        Par(i,:)=Off(k,:);
    end
    
    %�.�. ��� �������� ������ ������ ��������� �������� ������������ ��
    %������� ����� (����� ��������� ����� ��������� ������)
    k = ceil(rand()*M);
    Par(k,:) = S_opt;
    
    if min_u<min_u_glob % ������� ���������
        N_un = 0; 
        min_u_glob = min_u;
        xmin = T2(S_opt);
    else % ������� �� ���������
        N_un = N_un+1; 
    end
    
    iteration = iteration + 1;
end

end

