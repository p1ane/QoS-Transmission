function [Lambda_mat,sumRate,P_sum] = cal_Lambda2_para_cvx(Gamma,para1_R1,para2_R1,diff_of_R2,Nt,NumUsers,P)
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
    temp1(k) = sum(diag(diff_of_R2(:,:,k)) .* (Lambda_ele(:,k)));
end

expression Rse(NumUsers)
for k = 1:NumUsers
    Rse(k) = R1_de(k) - temp1(k);
end

expression sumR
sumR = sum(Rse);

expression sum_Lambda
for k = 1:NumUsers
    sum_Lambda = sum_Lambda + ...
        sum(Lambda_ele(:,k));
end

maximize(sumR)

subject to

sum_Lambda <= P;
Lambda_ele >= 0;

cvx_end;


for k = 1:NumUsers
    Lambda_mat(:,:,k) = diag(Lambda_ele(:,k));
end

sumRate = sumR;

P_sum = sum_Lambda;
end

