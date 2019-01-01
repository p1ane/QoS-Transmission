function [diff_of_logdet_Keve,diff_of_logdet_Kk] = cal_diffpara_of_cvx(CCM,CCMev,Lambda_l,Nt,Nr,NumUsers)


diff_of_logdet_Keve = zeros(Nt,Nt,NumUsers);
diff_of_logdet_Kk = zeros(Nt,Nt,NumUsers);


for k = 1:NumUsers
[diff_of_logdet_Keve(:,:,k)] = cal_diff_of_Kevek(CCMev,Lambda_l(:,:,k),Nt,Nr);
[diff_of_logdet_Kk(:,:,k)] = cal_diff_of_Kk_to_lambdai(CCM,Lambda_l,NumUsers,Nt,Nr,k);
end


end

