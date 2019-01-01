function [K] = cal_K_final(CCM,Lambda,Nr,NumUsers)
%checked

K = zeros(Nr,Nr,NumUsers);

for k = 1:NumUsers
    
    for i = 1:NumUsers
        
        if(i ~= k)
            K(:,:,k) = K(:,:,k) + cal_Ak_to_lambdai(CCM(:,:,k),Lambda(:,:,i),Nr);
        end
        
        if(i == k)
            K(:,:,k) = K(:,:,k) + 0;
        end
        
    end
    
    K(:,:,k) = eye(Nr) + K(:,:,k);
    
end

end

