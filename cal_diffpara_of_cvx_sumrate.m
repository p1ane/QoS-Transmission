function [diff_of_R2] = cal_diffpara_of_cvx_sumrate(CCM,CCMev,Lambda_l,Nt,Nr,NumUsers)

diff_of_logdet_Keve = zeros(Nt,Nt,NumUsers);
diff_of_logdet_Kk = zeros(Nt,Nt,NumUsers);
diff_of_R2 = zeros(Nt,Nt,NumUsers);
temp = zeros(Nt,Nt,NumUsers);

for k = 1:NumUsers
[diff_of_logdet_Keve(:,:,k)] = cal_diff_of_Kevek(CCMev,Lambda_l(:,:,k),Nt,Nr);
[diff_of_logdet_Kk(:,:,k)] = cal_diff_of_Kk_to_lambdai(CCM,Lambda_l,NumUsers,Nt,Nr,k);
end

for i = 1:NumUsers
    
    for k = 1:NumUsers
        
        if(k == i)
            temp(:,:,i) = temp(:,:,i);
        end
        
        if(k ~= i)
            temp(:,:,i) = temp(:,:,i) + diff_of_logdet_Kk(:,:,k);
        end
        
    end
    
    diff_of_R2(:,:,i) = diff_of_logdet_Keve(:,:,i) + temp(:,:,i);
    
end

      
end

