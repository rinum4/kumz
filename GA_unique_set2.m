function [c, ceq] = GA_unique_set2(Sga,n)
% 
% 
nk = numel(Sga)/n;
for k=1:nk
    xk=Sga(1+(k-1)*n:k*n);
    if sum(xk) == 1
        ceq(1+(k-1)*n:k*n)=0;
    else
        ceq(1+(k-1)*n:k*n)=100;
    end    
end
c=[;];
%ceq = [];
end

