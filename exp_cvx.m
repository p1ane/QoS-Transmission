



clear all;
clc;
[Nr,Nt,NumUsers,NumSamples,NumLinks,H_freq,H_freq_beam,CCM] = gen_of_channel;
[Hev_freq,Hev_freq_beam,CCMev] = gen_of_evchannel;
[H_freq_beam_test,CCM_test] = reform_H_freq_beam(CCM,Nt,Nr,NumLinks,NumSamples);
[Hev_freq_beam_test,CCMev_test] = reform_Hev_freq_beam(CCMev,Nt,Nr,NumSamples);


P = 100;
Lambda_l = zeros(Nt,Nt,NumUsers);
for i = 1:NumUsers
    Lambda_l(:,:,i) = (P/Nt/NumUsers).*eye(Nt);
end



[PHi_b,PHi,Gamma,Gamma_b,K] = cal_para_of_R1de(CCM,Lambda_l,Nt,Nr,NumUsers,1e-4);
[Kev] = cal_Kev(CCMev,Lambda_l,Nr,NumUsers);
[R2_l] = cal_R2(K,Kev,NumUsers);
[diff_of_logdet_Keve,diff_of_logdet_Kk] = cal_diffpara_of_cvx(CCM,CCMev,Lambda_l,Nt,Nr,NumUsers);


cvx_begin SDP quiet

variable Lambda(Nt,Nt,NumUsers) hermitian semidefinite diagonal;


expression R1_de(NumUsers)
for k = 1:NumUsers
    R1_de(k) = log_det(eye(Nt) + Gamma(:,:,k)*Lambda(:,:,k)) + log_det(Gamma_b(:,:,k) + K(:,:,k)) - trace(eye(Nr) - eye(Nr)/PHi_b(:,:,k));
end


expression temp1(NumUsers)
for k = 1:NumUsers
    temp1(k) = trace(diff_of_logdet_Keve(k)*(Lambda(:,:,k) - Lambda_l(:,:,k)));
end


expression temp2(NumUsers)
for k = 1:NumUsers
    
    
    for i = 1:NumUsers
        
        
        if(i == k)
            temp2(NumUsers) = temp2(NumUsers);
        end
        
        
        if(i ~= k)
            temp2(NumUsers) = temp2(NumUsers) + ...
                trace(diff_of_logdet_Kk(:,:,k)*(Lambda(:,:,i) - Lambda_l(:,:,i)));
        end
        
        
    end
    
    
end


expression Rse(NumUsers)
for k = 1:NumUsers
    Rse(k) = R1_de(k) - R2_l(k) - ...
        temp1(k) - temp2(k);
end


expression R
R = min(Rse);


expression sum_Lambda
for k = 1:NumUsers
    sum_Lambda = sum_Lambda + ...
        trace(Lambda(:,:,k));
end


maximize(R)


subject to


sum_Lambda <= P;


cvx_end;

[K_test] = cal_K_final(CCM,Lambda,Nr,NumUsers);
[R1_de_test,count] = cal_R1_de(CCM,Lambda,K,Nt,Nr,NumUsers,1e-4);
[Kev_test] = cal_Kev(CCMev,Lambda,Nr,NumUsers);
[R2_test] = cal_R2(K_test,Kev_test,NumUsers);
R_test = min(R1_de_test - R2_test);












