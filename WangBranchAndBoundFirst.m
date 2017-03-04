function [record,xmin]=WangBranchAndBoundFirst(a1,b1,c1,S1,n,a,b,c)
% F2|(1,2,1)-reentrant|Cmax
% Wang et al., 1997

%F=0;
%N = 1:n; %тест
N = S1;
record = inf;
xmin=[];

i=1; %1
while i<=n%n
    %револьвер
%     if i>1
%         i_=N(1);
%         for j=1:n-1
%             N(j)=N(j+1);
%         end    
%         N(n)=i_;
%     end
    
    N1 = N;
    %J1=[];
    Jk= N(i); %N(1) при револьвере
    %J2=[];
    N1(i) = []; %N(1) при револьвере
    %lvl = 1;
    
    [record,xmin]=WangBranchAndBoundBeyond(a1,b1,c1,N1,[],Jk,[],n, a, b, c,record,xmin);
    %F=WangBranchAndBoundBeyond3(a1,b1,c1,N1,Jk,2,n); %альтернативный
    %перебор
    i=i+1;
    
end %for i

end
