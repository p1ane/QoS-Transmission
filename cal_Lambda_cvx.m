function [Lambda_mat,Rse_lb] = cal_Lambda_cvx(PHi_b,Gamma,Gamma_b,K,Lambda_l,R2_l,diff_of_logdet_Keve,diff_of_logdet_Kk,Nt,Nr,NumUsers,P)


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

Lambda_mat = Lambda;

Rse_lb = R;

end

