CCM = ones(4,128,8);
CCMev = ones(4,128);
Lambda = zeros(128,128,4);
for i = 1:8
    Lambda(:,:,i) = eye(128);
end
% 
[diff_of_R2] = cal_diffpara_of_cvx_sumrate(CCM,CCMev,Lambda,128,4,8);

% CCMev = ones(4,128);
% Lambda = eye(128);
% [diff_of_Kevek] = cal_diff_of_Kevek(CCMev,Lambda,128,4);



