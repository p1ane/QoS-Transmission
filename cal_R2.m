function [R2] = cal_R2(K,Kev,NumUsers)

R2 = zeros(NumUsers,1);

for k = 1:NumUsers
%     R2(k) = log(det(K(:,:,k))) + log(det(Kev(:,:,k)));
    R2(k) = sum(log(diag(K(:,:,k)))) + sum(log(diag(Kev(:,:,k))));
end

end

