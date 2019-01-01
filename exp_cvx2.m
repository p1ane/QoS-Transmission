clear all;
clc;
[Nr,Nt,NumUsers,NumSamples,NumLinks,H_freq,H_freq_beam,CCM,R_b_ele] = gen_of_channel;
[Hev_freq,Hev_freq_beam,CCMev,R_b_ele_ev] = gen_of_evchannel;
[H_freq_beam_test,CCM_test] = reform_H_freq_beam(CCM,Nt,Nr,NumLinks,NumSamples);
[Hev_freq_beam_test,CCMev_test] = reform_Hev_freq_beam(CCMev,Nt,Nr,NumSamples);


P = 100;
Lambda_l = zeros(Nt,Nt,NumUsers);
for i = 1:NumUsers
    Lambda_l(:,:,i) = (P/Nt/NumUsers).*eye(Nt);
end

Lambda_l = zeros(Nt,Nt,NumUsers);
for i = 1:NumUsers
    Lambda_l(:,:,i) = (P/Nt).*eye(Nt);
    Lambda_l(12,12,i) = Lambda_l(12,12,i) - 0.5;
    Lambda_l(64,64,i) = Lambda_l(64,64,i) + 0.2;
    Lambda_l(100,100,i) = Lambda_l(100,100,i) + 0.3;
end



[PHi_b,PHi,Gamma,Gamma_b,K] = cal_para_of_R1de(H_freq_beam,CCM,Lambda_l,Nt,Nr,NumUsers,NumSamples,1e-4);


[R2_l] = cal_Rk2(H_freq_beam,Hev_freq_beam,Nr,Lambda_l,NumUsers,NumSamples);
[diff_of_logdet_Keve,diff_of_logdet_Kk] = cal_diffpara_of_cvx(CCM,CCMev,Lambda_l,Nt,Nr,NumUsers);


[Lambda_mat,Rse_lb] = cal_Lambda_cvx(PHi_b,Gamma,Gamma_b,K,Lambda_l,R2_l,diff_of_logdet_Keve,diff_of_logdet_Kk,Nt,Nr,NumUsers,P);



