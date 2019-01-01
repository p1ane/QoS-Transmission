lambda_ele = zeros(Nt,NumUsers);

for k = 1:NumUsers
    lambda_ele(:,k) = diag(Lambda_temp_new(:,:,k));
end

P_test_sum = sum(sum(lambda_ele));

% [R1_de,count] = cal_Rk1_de(H_freq_beam,CCM,Lambda_temp_new,Nt,Nr,NumUsers,NumSamples,1e-4);
% [R2] = cal_Rk2(H_freq_beam,Hev_freq_beam,Nr,Lambda_temp_new,NumUsers,NumSamples);
% 
% R_se = R1_de - R2;