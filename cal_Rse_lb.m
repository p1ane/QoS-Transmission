function [Rse_lb] = cal_Rse_lb(CCM,CCMev,Lambda,Nt,Nr,NumUsers)

R = zeros(NumUsers,1);
[~,~,~,~,K] = cal_para_of_R1de_exp(CCM,Lambda,Nt,Nr,NumUsers,1e-5);
[R1_de,~] = cal_R1_de(CCM,Lambda,K,Nt,Nr,NumUsers,1e-5);
[Kev] = cal_Kev(CCMev,Lambda,Nr,NumUsers);

for k = 1:NumUsers
    R(k) = R1_de(k) - ...
        log(det(K(:,:,k))) - log(det(Kev(:,:,k)));
end

Rse_lb = min(R);




end

