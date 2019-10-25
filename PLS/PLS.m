function [p,cost]=PLS(h,C,plower,pupper,pd)
n=length(plower);
p=zeros(n,7);
for i=1:n
    cost(i,:)=h(i,:)*C(1,i);
end
for i=1:n
    fla(i,1)=(cost(i,1)+cost(i,2)*pupper(1,i)+cost(i,3)*pupper(1,i)^2)/pupper(1,i);
end
p(:,7)=1:n;
p(:,1)=fla;
p(:,2)=plower';
p(:,3)=pupper';
p=sortrows(p);
p(1,4)=p(1,2);
p(1,5)=p(1,3);
for i=2:n
    p(i,4)=p(i,2)+p(i-1,4);
    p(i,5)=p(i,3)+p(i-1,5);
end
for i=1:n
    if pd>=p(i,5)
        p(i,6)=p(i,3);
    else
        g=0;
        for j=1:i-1
            g=g+p(j,3);
        end
        p(i,6)=pd-g;
    end
end
for i=1:n
    if p(i,6) < 0
        p(i,6)=0;
    end
end
cost=sum( [sum(cost(i,1)) sum(cost(i,2).*p(:,6)) sum(cost(i,3).*p(:,6).*p(:,6))]);
end