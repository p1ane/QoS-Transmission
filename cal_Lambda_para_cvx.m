function [Lambda_mat,Rse_lb,P_sum] = cal_Lambda_para_cvx(Gamma,para1_R1,para2_R1,Lambda_l_ele,R2_l,diff_of_logdet_Keve,diff_of_logdet_Kk,Nt,NumUsers,P)
digits(50);
cvx_begin  quiet

variable Lambda_ele(Nt,NumUsers);

expression R1_de(NumUsers)
for k = 1:NumUsers
    R1_de(k) = sum(log(1 + diag(Gamma(:,:,k)) .* (Lambda_ele(:,k)))) + ...
        para1_R1(k) - para2_R1(k);
end

expression temp1(NumUsers)
for k = 1:NumUsers
    temp1(k) = sum(diag(diff_of_logdet_Keve(k)) .* (Lambda_ele(:,k) - Lambda_l_ele(:,k)));
end

expression temp2(NumUsers)
for k = 1:NumUsers
    
    for i = 1:NumUsers
        
        
        if(i == k)
            temp2(NumUsers) = temp2(NumUsers);
        end
        
        
        if(i ~= k)
            temp2(NumUsers) = temp2(NumUsers) + ...
                sum(diag(diff_of_logdet_Kk(:,:,k)) .* (Lambda_ele(:,i) - Lambda_l_ele(:,i)));
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
        sum(Lambda_ele(:,k));
end

maximize(R)

subject to

sum_Lambda <= P;
Lambda_ele >= 0;

cvx_end;


for k = 1:NumUsers
    Lambda_mat(:,:,k) = diag(Lambda_ele(:,k));
end

Rse_lb = R;

P_sum = sum_Lambda;
end

