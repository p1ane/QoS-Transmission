function [eta] = cal_eta(Nt,Nr,k,CCM,X)
eta = zeros(Nt,Nt);
eta_ele = zeros(Nt,1);
temp = diag(X);



for m = 1:Nt
    
    for n = 1:Nr
        eta_ele(m) = eta_ele(m) + ...
            CCM(n,m,k)*temp(n);
    end
    
    
end

% for m = 1:Nt
%     eta_ele(m) = trace(diag(CCM(:,m,k)*X));
% end

    
eta(:,:) = diag(eta_ele(:,1));



end

