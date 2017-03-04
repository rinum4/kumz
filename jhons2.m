function [n] = jhons2(p,N)
% алгоритм Джонсона
%N - доступные значения заданы
% F2|(perm), (pmtn)|Cmax

%для теста
%p=[4.1 5 6 4; 2 9 8 5]'; %[4 3.1 3 1 8; 8 3 3 4 7]
[nr,nc]=size(p);
%N=1:nr;

n1=[];
n2=[];

N1=N;
p1=p;
k=0;
while ~isempty(N1)
    if length(N1)>1
        [C1,I1]=min(p1);
        [C2,I2]=min(C1);
    else
        I1=[1 1];
        [C1,I2]=min(p1);
    end
    if I2==1
        j=I1(1);
        %n1=horzcat(n1,N1(j));
        n1=[n1 N1(j)];
    else
        j=I1(2);
        %n2=horzcat(n2,N1(j));
        n2=[N1(j) n2];
    end
    N1(j)=[];
    p1(j,:)=[  ];
    k=k+1;
    %if length(N1)==1
    %   n1=horzcat(n1,N1(1));
    %   N1(j)=[];
    %end
end
n=horzcat(n1,n2);
end