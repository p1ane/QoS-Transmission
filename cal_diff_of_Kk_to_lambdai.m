function [diff_of_Kk] = cal_diff_of_Kk_to_lambdai(CCM,Lambda,NumUsers,Nt,Nr,k)
%CCM是矩阵，即CCM;
%Lambda是包含所有功率矩阵的三维矩阵
%输出是一个二维矩阵
% temp_a = zeros(NumUsers,1);
temp_b = zeros(Nr,1);
temp_c = zeros(Nt,1);
diff_of_Kk_ele = zeros(Nt,1);


for m = 1:Nt
    
    temp_c = zeros(Nt,1);
    
    for n = 1:Nr
        
        temp_a = zeros(NumUsers,1);
        temp_b = zeros(Nr,1);
        
        for j = 1:NumUsers
            
            
            for q = 1:Nt
                temp_a(j) = temp_a(j) + CCM(n,q,k)*Lambda(q,q,j);
            end
            
            
            if(j==k)
                temp_b(n) = temp_b(n) + 0;
            end
            if(j~=k)
                temp_b(n) = temp_b(n) + temp_a(j);
            end
            
            
        end
        temp_c(m) = temp_c(m) + CCM(n,m,k)/(1 + temp_b(n));
    end
    
    
    diff_of_Kk_ele(m) = temp_c(m);
end


diff_of_Kk = diag(diff_of_Kk_ele);
end

