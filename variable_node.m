function extrinsic_LLRs = variable_node(priori_LLRs)
N= length(priori_LLRs);
extrinsic_LLRs = zeros(N);
for i = 2:N
    extrinsic_LLRs(i)=sum(priori_LLRs)-priori_LLRs(i);
end
end