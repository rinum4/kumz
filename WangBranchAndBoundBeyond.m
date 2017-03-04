function [record,xmin]=WangBranchAndBoundBeyond(a1,b1,c1,N,J1,Jk,J2,n, a, b, c,record,xmin)
% F2|(1,2,1)-reentrant|Cmax
% Wang et al., 1997

%global record;
%global xmax;
%global Jk_max;
%global xcount;
%global ccount;

%record=record;

%i=1;
if numel(N)>0
    for i=1:2
        N1 = N;
        %�������� ������ ��������� ������� ���� ���� � J1 ���� � J2
        if i==1
            J1_ = [J1 N(1)];
            %����� ����� ������� - ������������� �� �������� �� a,b
            p1=[a1(J1_);b1(J1_)]';
            J1_ = jhons2(p1,J1_);
            
            J2_=J2;
            %����� ����� ������� - ������������� �� �������� �� b,c
            p2=[b1(J2_);c1(J2_)]';
            J2_ = jhons2(p2,J2_);
            
            N1(1) = [];
            %��� ����� J1
            J1_1=[N1 Jk J2_];
            %J = [J1_ Jk J2_];
            
            %�������� �� inadmissible C2(J1) > a
            [~, C2, C3] = WangC2Fun(a1, b1, c1, J1_);
            if C2>a
                continue;
            end
            
            %������ �������
            LB1 = C2 + WangLB1Fun(b1, c1, J1_1);  % Wang(1997)
            %LB2 = C1 + WangLB2Fun(c1, J1_1);  % Wang(1997) %������
            %����� ������ � �������
            LB2 = C3 + WangLB2Fun(c1, J1_1);  % Wang(1997) %������ �������
            LB = [LB1 LB2];
            [LB,~]=max(LB); %������ 1
            %[LB,ind]=min(LB); %������ 2
            
            %LB=8;
            
            if LB<record %��� min <
                [record,xmin] = WangBranchAndBoundBeyond(a1,b1,c1,N1,J1_,Jk,J2_,n, a, b, c,record,xmin);
            else
                %record=record;
            end
        else
            J1_ = J1;
            %����� ����� ������� - ������������� �� �������� �� a,b
            p1=[a1(J1_);b1(J1_)]';
            J1_ = jhons2(p1,J1_);
            
            J2_ = [J2 N(1)];
            %����� ����� ������� - ������������� �� �������� �� b,c
            p2=[b1(J2_);c1(J2_)]';
            J2_ = jhons2(p2,J2_);
            
            N1(1) = [];
            %��� ����� J1
            J1_1=[N1 Jk J2_];
            %J = [J1_ Jk J2_];
            
            %�������� �� inadmissible C2(J1) > a
            [~, C2, C3] = WangC2Fun(a1, b1, c1, J1_);
            if C2>a
                continue;
            end
            
            %������ �������
            LB1 = C2 + WangLB1Fun(b1, c1, J1_1);  % Wang(1997)
            %LB2 = C1 + WangLB2Fun(c1, J1_1);  % Wang(1997) %������
            %����� ������ � �������
            LB2 = C3 + WangLB2Fun(c1, J1_1);  % Wang(1997) %������ �������
            LB = [LB1 LB2];
            [LB,~]=max(LB);
            
            %LB=8;
            
            if LB<record %��� min <
                [record,xmin]=WangBranchAndBoundBeyond(a1,b1,c1,N1,J1_,Jk,J2_,n, a, b, c,record,xmin);
            else
                %record=record;
            end
        end
        
     end %for i
else
    %     �������
    %     ccount=ccount+1;
    %     xcount(ccount,:)=[J1 Jk J2];
    
    J = [J1 Jk J2];
    FF = WangTargFun(a1,b1,c1,J);
    if FF<record
        xmin=J;
        record = FF;
        %Jk_max = Jk;
    end
end

end
