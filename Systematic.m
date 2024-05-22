function Q = Systematic(H)
[r,~]=size(H);
Q = H;
for y=1:r
    if Q(y,y)~=1
        n =find(Q(:,y)==1);
        if length(n)>1
            n = n(end);
        end
        a = Q(y,:);
        Q(y,:)=Q(n,:);
        Q(n,:) = a;
    end
    for x=1:r
        if Q(x,y)==1 && x~=y
            Q(x,:)=mod(Q(x,:)+Q(y,:),2);
        end
    end
end
end