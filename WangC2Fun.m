function [C1, C2, C3]=WangC2Fun(a1, b1, c1, S)

ns=numel(S);
C2=0;
a_s=sum(a1);
%b_s=sum(b1);
%c_s=sum(c1);
constr=a_s; %под вопросом
%e1=zeros(1,ns);
e1=0;
%e2=zeros(1,ns);
e2=0;
e3=0;
%e3=zeros(1,ns);
%phi=zeros(1,ns);
%phi2=zeros(1,ns);

for i=1:ns
    j=S(i);
    if i==1
        %e1(i)=a1(j);
        %e2(i)=e1(i)+b1(j);
        e1=a1(j);
        e2=e1+b1(j);
        if ns==1
            if e2>constr
                e3=e2;
            elseif i==ns
                e3=constr;
            end
            if e3>0
                e3=e3+c1(j);
           end
        end    
    else
        %e1(i)=e1(i-1)+a1(j);
        %phi(i)=max(e1(i),e2(i-1));
        %e2(i)=phi(i)+b1(j);
        e1=e1+a1(j);
        phi=max(e1,e2);
        e2=phi+b1(j);
        constr=constr+c1(S(i-1));
        
        if e3 == 0
            if e2>constr
                e3=e2;
            elseif i==ns
                e3=constr;
            end
            if e3>0
                e3=e3+c1(j);
           end
        else
            phi2=max(e3,e2);
            e3=phi2+c1(j);
        end
    end
end
C1 = e1;
C2 = e2;
C3 = e3;
end