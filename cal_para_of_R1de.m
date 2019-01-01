function [para1_R1,para2_R1,Gamma,K] = cal_para_of_R1de(CCM,Lambda_l,Nt,Nr,NumUsers,it_con)
PHi_b = zeros(Nr,Nr,NumUsers);
PHi_b_old = zeros(Nr,Nr,NumUsers);
PHi = zeros(Nt,Nt,NumUsers);
Gamma = zeros(Nt,Nt,NumUsers);
Gamma_b = zeros(Nr,Nr,NumUsers);
test = zeros(NumUsers,1);
para1_R1 = zeros(NumUsers,1);
para2_R1 = zeros(NumUsers,1);

for k = 1:NumUsers
    PHi_b(:,:,k) = eye(Nr);
end


[K] = cal_K_final(CCM,Lambda_l,Nr,NumUsers);


while 1 
    
    
    for k = 1:NumUsers
        temp = inv(PHi_b(:,:,k))*inv(K(:,:,k));
        Gamma(:,:,k) = cal_eta(Nt,Nr,k,CCM,temp);
        PHi(:,:,k) = eye(Nt) + Gamma(:,:,k)*Lambda_l(:,:,k);
        %temp2 = inv(PHi(:,:,k))*Lambda(:,:,k);%
        temp2 = PHi(:,:,k)\Lambda_l(:,:,k);
        Gamma_b(:,:,k) = cal_eta_b(Nt,Nr,k,CCM,temp2);
        PHi_b_old(:,:,k) =  PHi_b(:,:,k);
        PHi_b(:,:,k) = eye(Nr) + Gamma_b(:,:,k)/K(:,:,k);
        test(k) = norm(PHi_b_old(:,:,k) - PHi_b(:,:,k),2);
    end
    
    if(max(test) <= it_con)
        break;
    end
    
    
end

Gamma = real(Gamma);
PHi_b = real(PHi_b);
Gamma_b = real(Gamma_b);
K = real(K);

for k = 1:NumUsers
%     para1_R1(k) = log(det(Gamma_b(:,:,k) + K(:,:,k)));
    para1_R1(k) = sum(log(diag(Gamma_b(:,:,k) + K(:,:,k))));
    para2_R1(k) = trace(eye(Nr) - eye(Nr)/PHi_b(:,:,k));
end
 
end

