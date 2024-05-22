function extrinsic_LLRs = check_node(priori_LLRs)


N= length(priori_LLRs);
f=zeros(size(priori_LLRs)-1);
for i=1:N-1
    if i==1
        f(i)=priori_LLRs(i);
    elseif i > 1
        f(i)=boxplus(priori_LLRs(i),f(i-1));
    end
end

b=zeros(size(priori_LLRs)-1);
for i=N:-1:2
    if i==N
        b(i-1)=priori_LLRs(N);
    elseif i < N
        b(i-1)=boxplus(priori_LLRs(i),b(i));    
    end
end

extrinsic_LLRs = zeros(size(priori_LLRs));
%extrinsic_LLRs(i) = the BOXPLUS sum of all priori_LLRs EXCEPT priori_LLRs(i)
   for i=1:N
       if i==1
           extrinsic_LLRs(i)=b(1);
      
       elseif i==N
           extrinsic_LLRs(i)=f(N-1);
       else
           extrinsic_LLRs(i)=boxplus(f(i-1),b(i));  
       end
   end
end
   