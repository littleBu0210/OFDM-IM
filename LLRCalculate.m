function LLR = LLRCalculate(yF,NsubCarry,NactiveCarry,h,s_ref,N_F)

OutpuSum = zeros(size(yF),'double');
for i =1:length(s_ref)
    OutpuSum = OutpuSum+exp(-1*conj(yF-h.*s_ref(i)).*(yF-h.*s_ref(i))/N_F);
end

LLR = log(NactiveCarry)-log(NsubCarry-NactiveCarry)+yF.*conj(yF)/N_F+log(OutpuSum)-log(length(s_ref));

% LLR = log(NactiveCarry)-log(NsubCarry-NactiveCarry)+yF.*conj(yF)/N_F+log(OutpuSum);
end