function [Aki] = cal_Ak_to_lambdai(CCMk,Lambdai,Nr)
%Lambdai是单用户即Lambda(:,:,k);
%CCMk也是单用户的CCM即CCM(:,:,k)
%已验证
Aki_ele = zeros(Nr,1);

for i = 1:Nr
    Aki_ele(i) = trace(diag(CCMk(i,:)*Lambdai));
end


Aki = diag(Aki_ele);


end

