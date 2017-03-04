function [LB]=WangLB1Fun(b1,c1, S)

ns=numel(S);
SumB=0;
MinC=inf;

for i=1:ns
    j=S(i);
    
    SumB=SumB+b1(j);
    if c1(j)<MinC
        MinC = c1(j);
    end
end

LB=SumB + MinC;

end