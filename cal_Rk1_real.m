function [Rk1_real] = cal_Rk1_real(H_freq_beam,Nr,Lambda,NumUsers,NumSamples)
%CAL_RK1_REAL calculate the real value of the Rk
Rk1_real = zeros(NumUsers,1);
K = zeros(Nr,Nr,NumUsers);

[K] = cal_K_new(H_freq_beam,Nr,Lambda,NumUsers,NumSamples);



for k = 1:NumUsers
    R_temp = 0;
    
    
    for Sample_n = 1:NumSamples
        R_temp = R_temp + ...
            log(det(K(:,:,k) + H_freq_beam(:,:,1,Sample_n,k)*Lambda(:,:,k)*H_freq_beam(:,:,1,Sample_n,k)'));
    end
    
    
    Rk1_real(k) = real(R_temp/NumSamples);
end



end







