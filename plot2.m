clear all;
clc;
close all;


addpath('.\data');
load('CCM2.mat');
load('CCMev2.mat');
load('H_freq_beam2');
load('Hev_freq_beam2');
load('Lambda_-10_202');
load('R_b2.mat');


%参量获取
[Nr,Nt,x,NumSamples,NumUsers] = size(H_freq_beam_test);
K = zeros(Nr,Nr,NumUsers,7);
Kev = zeros(Nr,Nr,NumUsers,7);
Lambda = zeros(Nt,Nt,NumUsers,7);
R = zeros(NumSamples,NumUsers,7);
R_b = zeros(NumSamples*NumUsers,7);


for i = 1:7
    
    for k = 1:NumUsers
        
        Lambda(:,:,k,i) = diag(Lambda_final_ele(:,k,i));
        
    end
    
    
    [K(:,:,:,i)] = cal_K_final(CCM, Lambda(:,:,:,i),Nr,NumUsers);
    [Kev(:,:,:,i)] = cal_Kev(CCMev,Lambda(:,:,:,i),Nr,NumUsers);
    
end


for i = 1:7
    
    for k = 1:NumUsers
        
        for sam = 1:NumSamples
            
            R(sam,k,i) = real(log(det(K(:,:,k,i) + H_freq_beam_test(:,:,1,sam,k)*Lambda(:,:,k,i)*H_freq_beam_test(:,:,1,sam,k)')) - ...
                log(det(K(:,:,k,i))) - log(det(Kev(:,:,k,i))));
            
        end
        
    end
    
end


for i = 1:7
    num = 1;
    for sam = 1:NumSamples
        
        for k = 1:NumUsers
            
            R_b(num,i) = R(sam,k,i);
            num = num + 1;
        end
        
    end
    
end

figure;
hold on;
%标记是marker

hfig(:,1) = cdfplot(R_b(:,1));
set(hfig(:,1),'Color', 'r','LineStyle','-','LineWidth',2);
hfig(:,2) = cdfplot(R_b1(:,1));

hfig(:,3) = cdfplot(R_b(:,3));
set(hfig(:,3),'Color', 'r','LineStyle','-','LineWidth',2);
hfig(:,4) = cdfplot(R_b1(:,3));

hfig(:,5) = cdfplot(R_b(:,5));
set(hfig(:,5),'Color', 'r','LineStyle','-','LineWidth',2);
hfig(:,6) = cdfplot(R_b1(:,5));

fig_legend_str{1,1}=('-10 dB');
fig_legend_str{1,2}=('-10 dB');
fig_legend_str{1,3}=('0 dB');
fig_legend_str{1,4}=('0 dB');
fig_legend_str{1,5}=('10 dB');
fig_legend_str{1,6}=('10 dB');

legend(hfig(1,:),fig_legend_str,'Location','SouthEast');


