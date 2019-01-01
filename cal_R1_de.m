function [R1_de,count] = cal_R1_de(CCM,Lambda,K,Nt,Nr,NumUsers,it_con)
PHi_b = zeros(Nr,Nr,NumUsers);
PHi_b_old = zeros(Nr,Nr,NumUsers);
PHi = zeros(Nt,Nt,NumUsers);
Gamma = zeros(Nt,Nt,NumUsers);
Gamma_b = zeros(Nr,Nr,NumUsers);
test = zeros(NumUsers,1);
R1_de = zeros(NumUsers,1);
count = 0;

for k = 1:NumUsers
    PHi_b(:,:,k) = eye(Nr);
end


while 1 
    
    
    for k = 1:NumUsers
        temp = inv(PHi_b(:,:,k))*inv(K(:,:,k));
        Gamma(:,:,k) = cal_eta(Nt,Nr,k,CCM,temp);
        PHi(:,:,k) = eye(Nt) + Gamma(:,:,k)*Lambda(:,:,k);
        %temp2 = inv(PHi(:,:,k))*Lambda(:,:,k);%
        temp2 = PHi(:,:,k)\Lambda(:,:,k);
        Gamma_b(:,:,k) = cal_eta_b(Nt,Nr,k,CCM,temp2);
        PHi_b_old(:,:,k) =  PHi_b(:,:,k);
        PHi_b(:,:,k) = eye(Nr) + Gamma_b(:,:,k)/K(:,:,k);
        test(k) = norm(PHi_b_old(:,:,k) - PHi_b(:,:,k),2);
    end
    count = count + 1;
    
    if(max(test) <= it_con)
        break;
    end
    
    
end


for k = 1:NumUsers
%     R1_de(k) = log(det(eye(Nt) + Gamma(:,:,k)*Lambda(:,:,k))) + ...
%         log(det(Gamma_b(:,:,k) + K(:,:,k))) - trace(eye(Nr) - eye(Nr)/PHi_b(:,:,k));
    
    R1_de(k) = sum(log(diag(eye(Nt) + Gamma(:,:,k)*Lambda(:,:,k)))) + ...
        sum(log(diag(Gamma_b(:,:,k) + K(:,:,k)))) - trace(eye(Nr) - eye(Nr)/PHi_b(:,:,k));
  
    
end


R1_de = real(R1_de);

end

