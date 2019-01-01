function [eta_b] = cal_eta_b(Nt,Nr,k,CCM,X)
eta_b = zeros(Nr,Nr);
eta_b_ele = zeros(Nr,1);
temp = diag(X);



for n = 1:Nr
    
    
    for m = 1:Nt
        eta_b_ele(n) = eta_b_ele(n) +  CCM(n,m,k)*temp(m);
    end
    
    
    
end


% for  n = 1:Nr
%     eta_b_ele(n) = trace(diag(CCM(n,:,k)*X));
% end





eta_b(:,:) = diag(eta_b_ele(:,1));






end

