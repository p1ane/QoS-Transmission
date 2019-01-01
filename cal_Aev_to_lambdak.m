function [Aev] = cal_Aev_to_lambdak(CCMev,Lambda,Nr)
%Lambda输入单个用户的即Lambda(:,:,k)
Aev_ele = zeros(Nr,1);


for i = 1:Nr
    Aev_ele(i) = trace(diag(CCMev(i,:)*Lambda));
end


Aev = diag(Aev_ele);


end

