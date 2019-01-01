clear all;
close all;
clc;



addpath('.\data');
SNRdB = 15;
a = length(SNRdB);
Rse_lb = zeros(a,1);
Rse_lb_test = zeros(a,1);
Rse_lb_temp_new = zeros(a,1);
Rse_lb_temp_old = zeros(a,1);
Rse_lb_diff = zeros(a,1);
count = 0;



%生成信道
[Nr,Nt,NumUsers,NumSamples,NumLinks,H_freq,H_freq_beam,CCM] = gen_of_channel;
[Hev_freq,Hev_freq_beam,CCMev] = gen_of_evchannel;
Lambda_final_ele = zeros(Nt,NumUsers,a);
% load('CCM_Nt_256.mat');
% load('CCMev_Nt_256.mat');
% load('Lambda_10.mat');
[H_freq_beam_test,CCM_test] = reform_H_freq_beam(CCM,Nt,Nr,NumLinks,NumSamples);
[Hev_freq_beam_test,CCMev_test] = reform_Hev_freq_beam(CCMev,Nt,Nr,NumSamples);
Lambda_temp_old_ele = zeros(Nt,NumUsers);

for snr_num = 1:a
    count = 0;
    %功率矩阵初始化
    P = 10^(SNRdB(snr_num)/10);
    Lambda_temp_old = zeros(Nt,Nt,NumUsers);
    Lambda_temp_new = zeros(Nt,Nt,NumUsers);
    
    
    %     Lambda_temp_new = Lambda_temp_old;
    
    for k = 1:NumUsers
        Lambda_temp_old(:,:,k) = ((P/Nt/NumUsers).*eye(Nt));
        Lambda_temp_new(:,:,k) = ((P/Nt/NumUsers).*eye(Nt));
    end
    
    %解cccp
    %         while 1
    while(count <= 200)
        %计算参量
        [para1_R1,para2_R1,Gamma,K] = cal_para_of_R1de(CCM,Lambda_temp_old,Nt,Nr,NumUsers,1e-4);
        [Kev] = cal_Kev(CCMev,Lambda_temp_old,Nr,NumUsers);
        [R2_temp_old] = cal_R2(K,Kev,NumUsers);
        [diff_of_logdet_Keve,diff_of_logdet_Kk] = cal_diffpara_of_cvx(CCM,CCMev,Lambda_temp_old,Nt,Nr,NumUsers);
        
        for k = 1:NumUsers
            Lambda_temp_old_ele(:,k) = diag(Lambda_temp_old(:,:,k));
        end
        
        %迭代赋值
        Rse_lb_temp_old(snr_num) = Rse_lb_temp_new(snr_num);
        Lambda_temp_old = Lambda_temp_new;
        
        %cvx解凸问题
        [Lambda_temp_new,Rse_lb_diff(snr_num),P_sum] = cal_Lambda_para_cvx(Gamma,para1_R1,para2_R1,Lambda_temp_old_ele,R2_temp_old,diff_of_logdet_Keve,diff_of_logdet_Kk,Nt,NumUsers,P);
        
        %测试参数
        Rse_lb(snr_num) = Rse_lb_temp_new(snr_num);
        [K_test] = cal_K_final(CCM,Lambda_temp_new,Nr,NumUsers);
        [R1_de_test,count2] = cal_R1_de(CCM,Lambda_temp_new,K,Nt,Nr,NumUsers,1e-4);
        [Kev_test] = cal_Kev(CCMev,Lambda_temp_new,Nr,NumUsers);
        [R2_test] = cal_R2(K_test,Kev_test,NumUsers);
        Rse_lb_test(snr_num) = min(R1_de_test - R2_test);
        Rse_lb_temp_new(snr_num) = Rse_lb_test(snr_num);
%         判断收敛
%                         if (abs(Rse_lb_temp_old(snr_num) - Rse_lb_temp_new(snr_num)) <= 1e-5 && count >= 20)
%                             break;
%                         end
%         
        count = count + 1;
        Rse_lb_test(snr_num)
        P_sum;
        Lambda_final_ele(:,:,snr_num) = Lambda_temp_old_ele;
    end
    snr_num
    Rse_lb_test(snr_num)
    %     Rse_lb(snr_num) = Rse_lb_temp_new(snr_num);
    %     [K_test] = cal_K_final(CCM,Lambda_temp_new,Nr,NumUsers);
    %     [R1_de_test,count2] = cal_R1_de(CCM,Lambda_temp_new,K,Nt,Nr,NumUsers,1e-4);
    %     [Kev_test] = cal_Kev(CCMev,Lambda_temp_new,Nr,NumUsers);
    %     [R2_test] = cal_R2(K_test,Kev_test,NumUsers);
    %     Rse_lb_test(snr_num) = min(R1_de_test - R2_test);
    
    
    
end