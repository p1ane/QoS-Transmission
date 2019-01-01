% clear all;
% clc;
% close all;
% 
% 
% addpath('.\data');
% load('CCM1.mat');
% load('CCMev1.mat');
% load('H_freq_beam3');
% load('Hev_freq_beam3');
% load('Lambda_-10_201');
% 
% 
% %参量获取
% [Nr,Nt,x,NumSamples,NumUsers] = size(H_freq_beam_test);
% Lambda = zeros(Nt,Nt,NumUsers,7);
% K = zeros(Nr,Nr,NumUsers,7);
% Kev = zeros(Nr,Nr,NumUsers,7);
% R = zeros(NumSamples,NumUsers,7);
% R_b = zeros(NumSamples*NumUsers,7);
% 
% 
% for i = 1:7
%     
%     for k = 1:NumUsers
%         
%         Lambda(:,:,k,i) = diag(lambda_iter_par(:,k,i));
%         
%     end
%     
%     [K(:,:,:,i)] = cal_K_final(CCM, Lambda(:,:,:,i),Nr,NumUsers);
%     [Kev(:,:,:,i)] = cal_Kev(CCMev,Lambda(:,:,:,i),Nr,NumUsers);
%     
% end
% 
% 
% for i = 1:7
%     
%     for k = 1:NumUsers
%         
%         for sam = 1:NumSamples
%             
%             R(sam,k,i) = real(log(det(K(:,:,k,i) + H_freq_beam_test(:,:,1,sam,k)*Lambda(:,:,k,i)*H_freq_beam_test(:,:,1,sam,k)')) - ...
%                 log(det(K(:,:,k,i))) - log(det(Kev(:,:,k,i))));
%             
%         end
%         
%     end
%     
% end
% 
% 
% for i = 1:7
%     num = 1;
%     for sam = 1:NumSamples
%         
%         for k = 1:NumUsers
%             
%             R_b(num,i) = R(sam,k,i);
%             num = num + 1;
%         end
%         
%     end
%     
% end

cdfplot(R_b(:,4));







% log(det(K(:,:,k,i))) - log(det(Kev(:,:,k,i)));
%[Kev] = cal_Kev(CCMev,Lambda,Nr,NumUsers)












