function [H_freq_beam_test,CCM_test] = reform_H_freq_beam(CCM,Nt,Nr,NumLinks,NumSamples)
%REFORM_H_FREQ_BEAM 依据已有的统计特性，生成方差均值恰好满足已有统计特性高斯信道
j = sqrt(-1);
H_freq_beam_test = zeros(Nr,Nt,1,NumSamples,NumLinks);
CCM_test = zeros(Nr,Nt,NumLinks); 
for link_n = 1:NumLinks
    for sample_n = 1:NumSamples
        H_freq_beam_test(:,:,1,sample_n,link_n) = ...
        sqrt(0.5*CCM(:,:,link_n)).*(randn(Nr,Nt) + j*randn(Nr,Nt));
    end
end
for link_n = 1:NumLinks
    for sample_n = 1:NumSamples
        CCM_test(:,:,link_n) = ... 
        CCM_test(:,:,link_n) + H_freq_beam_test(:,:,:,sample_n,link_n).*conj(H_freq_beam_test(:,:,:,sample_n,link_n));
    end
    CCM_test(:,:,link_n) = CCM_test(:,:,link_n)/NumSamples;
end
end

