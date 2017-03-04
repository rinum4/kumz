function [LB]=WangLB2Fun(c1, S)

ns=numel(S);
SumC=0;

for i=1:ns
    j=S(i);
    
    SumC=SumC+c1(j);
    
end

LB=SumC;

end