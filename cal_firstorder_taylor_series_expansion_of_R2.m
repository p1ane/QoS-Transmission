function [firstorder_taylor_series_expansion_of_R2] = cal_firstorder_taylor_series_expansion_of_R2(CCMev,CCM,H_freq_beam,Hev_freq_beam,Nr,Nt,Lambda,Lambda_l,NumUsers,NumSamples)
diff_of_Kevek = zeros(Nt,Nt,NumUsers);
temp_RB = zeros(NumUsers,1);
temp_Rk = zeros(NumUsers,1);
firstorder_taylor_series_expansion_of_R2 = zeros(NumUsers,1);


R2_l = cal_Rk2(H_freq_beam,Hev_freq_beam,Nr,Lambda,NumUsers,NumSamples);


for k = 1:NumUsers
    
    
    diff_of_Kevek(:,:,k) = cal_diff_of_Kevek(CCMev,Lambda(:,:,k),Nt,Nr);
    temp_RB(k) = trace(diff_of_Kevek(:,:,k)*(Lambda(:,:,k) - Lambda_l(:,:,k)));
    %注意转置操作未进行因为没有实际变化
    
    
end


for k = 1:NumUsers
    temp1 = 0;
    
    
    for i = 1:NumUsers
        if(i == k)
            temp1 = temp1 + 0;
        end
        if(i ~= k)
            diff_of_Kk_to_lambdai = cal_diff_of_Kk_to_lambdai(CCM(:,:,k),Lambda,NumUsers,Nt,Nr,k);
            temp1 = temp1 + trace(diff_of_Kk_to_lambdai*(Lambda(:,:,i) - Lambda_l(:,:,i)));
            %未进行转置操作因为没有实际意义
        end
    end
    
    
    temp_Rk(k) = temp1;
end


for k = 1:NumUsers
    firstorder_taylor_series_expansion_of_R2(k) = R2_l(k) + ...
        temp_RB(k) + temp_Rk(k);
end


end

