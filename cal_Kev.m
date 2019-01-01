function [Kev] = cal_Kev(CCMev,Lambda,Nr,NumUsers)


Kev = zeros(Nr,Nr,NumUsers);
for k = 1:NumUsers
    Kev(:,:,k) = eye(Nr) + cal_Aev_to_lambdak(CCMev,Lambda(:,:,k),Nr);
end

end

