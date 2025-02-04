function extrinsic_LLRs = check_node1(priori_LLRs)


N= length(priori_LLRs);
f=zeros(size(priori_LLRs));
for i=1:N
    if i==1
        f(i)=priori_LLRs(i);
    elseif i > 1
        f(i)=boxplus1(priori_LLRs(i),f(i-1));
    end
end

b=zeros(size(priori_LLRs));
for i=N:-1:1
    if i==N
        b(i)=priori_LLRs(N);
    elseif i < N
        b(i)=boxplus1(priori_LLRs(i),b(i+1));    
    end
end

extrinsic_LLRs = zeros(size(priori_LLRs));
%extrinsic_LLRs(i) = the BOXPLUS sum of all priori_LLRs EXCEPT priori_LLRs(i)
   for i=1:N
       if i==1
           extrinsic_LLRs(i)=b(2);
      
       elseif i==N
           extrinsic_LLRs(i)=f(N-1);
       else
           extrinsic_LLRs(i)=boxplus1(f(i-1),b(i+1));  
       end
   end
end
   