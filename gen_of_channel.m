function [Nr,Nt,NumUsers,NumSamples,NumLinks,H_freq,H_freq_beam,CCM] = gen_of_channel(~)
addpath('.\SCM');
j = sqrt(-1);
% Basic Parameters
% system parameters
Scenario = 'suburban_macro';
Nt = 256;                   % number of transmit antennas
Nr = 4;                     % number of receive antennas
NumUsers = 8;               % number of users in the cell and one passive eavesdrropper
NumSamples = 2000;         % number of samples in each channel realization
NumPaths = 4;               % number of paths
NumSubpaths = 20;           % number of subpaths in each path
CenterFrequency = 2e9;      % central frequency - 2G Hz
Ts = 1/(16*3.84e6);         % delay sampling interval                  % total transmit power unit:dBW
gamma = ones(NumUsers,1);   % rate weight
% OFDM parameters
NumSubcarriers = 256;
SubcarrierLength = 1;
SubcarrierStart = 0;
% SCM parameters
scmpar = scmparset;
scmpar.Scenario = Scenario;
scmpar.NumBsElements = Nt;
scmpar.NumMsElements = Nr;
scmpar.NumTimeSamples = NumSamples;
scmpar.NumPaths = NumPaths;
scmpar.NumSubPathsPerPath = NumSubpaths;
scmpar.CenterFrequency = CenterFrequency;
scmpar.DelaySamplingInterval = Ts;
scmpar.PathLossModelUsed = 'no';
scmpar.ShadowingModelUsed = 'no';
scmpar.AnsiC_core = 'yes';
% link parameters
NumLinks = (NumUsers);
linkpar = linkparset(NumLinks);
linkpar.ThetaMs = zeros(1,NumLinks);
% antenna parameters
antpar = antparset;
antpar.BsElementPosition = 1/2;
antpar.MsElementPosition = 1/2;
FM = dftmtx(Nt)/sqrt(Nt);
% fix point iteration prarameter
miu = 1e1;                      % approximation parameter in entropy regularization
epsilon = 1e-5;                 % terminate accuracy
NumIterMin = 100;                % minimum number of iteration
NumIterMax = 100;                % maximum number of iteration
% Generate channel matrix
[H_time, delay, out] = scm(scmpar,linkpar,antpar); % time-domain channel matrix - wideband
taps = delay/Ts;                                   % delay taps
CCM = zeros(Nr,Nt,NumLinks);                       % channel coupling matrix - average in all subcarriers
Covmat_t = zeros(Nt,Nt,NumLinks);                  % transmit channel covariance
Covmat_r = zeros(Nr,Nr,NumLinks);                  % receive channel covariance
Ur = zeros(Nr,Nr,NumLinks);
Ut = zeros(Nt,Nt,NumLinks);
H_freq = zeros(Nr,Nt,SubcarrierLength,NumSamples,NumLinks);% frequency-domain channel matrix (on first subcarrier)
H_freq_beam = zeros(Nr,Nt,SubcarrierLength,NumSamples,NumLinks);% frequency-domain beamspace channel matrix (on first subcarrier)
for Link_n = 1:NumLinks
    for Sample_n = 1:NumSamples
        for Path_n = 1:NumPaths
            for Subcarrier_n = 1:SubcarrierLength
                H_freq(:,:,Subcarrier_n,Sample_n,Link_n) = H_freq(:,:,Subcarrier_n,Sample_n,Link_n)+...
                    H_time(:,:,Path_n,Sample_n,Link_n)*exp(-j*2*pi*(Subcarrier_n-1+SubcarrierStart)*taps(Link_n,Path_n)/NumSubcarriers);
            end
        end
        Covmat_r(:,:,Link_n) = Covmat_r(:,:,Link_n) + H_freq(:,:,:,Sample_n,Link_n)*H_freq(:,:,:,Sample_n,Link_n)';
        Covmat_t(:,:,Link_n) = Covmat_t(:,:,Link_n) + H_freq(:,:,:,Sample_n,Link_n)'*H_freq(:,:,:,Sample_n,Link_n);    
    end
    Covmat_r(:,:,Link_n) = Covmat_r(:,:,Link_n)/NumSamples;
    Covmat_t(:,:,Link_n) = Covmat_t(:,:,Link_n)/NumSamples;
    [Ur(:,:,Link_n),~] = eig(Covmat_r(:,:,Link_n));
    [Ut(:,:,Link_n),~] = eig(Covmat_t(:,:,Link_n));
end
for Link_n = 1:NumLinks
   for Sample_n = 1:NumSamples
       for Subcarrier_n = 1:SubcarrierLength
           H_freq_beam(:,:,Subcarrier_n,Sample_n,Link_n) = Ur(:,:,Link_n)'*H_freq(:,:,Subcarrier_n,Sample_n,Link_n)*FM;
       end
       CCM(:,:,Link_n) = CCM(:,:,Link_n) + H_freq_beam(:,:,:,Sample_n,Link_n).*conj(H_freq_beam(:,:,:,Sample_n,Link_n));
   end
   CCM(:,:,Link_n) = CCM(:,:,Link_n)/NumSamples;
end
% initialization of Lambda
R_b = zeros(Nt,Nt,NumLinks);
R_b_ele = zeros(Nt,1,NumLinks);
for k = 1:NumLinks
    for m = 1:Nt
    R_b_ele(m,1,k) = sum(CCM(:,m,k));
    end
    R_b(:,:,k) = diag(R_b_ele(:,:,k));
end


%CCMπÈ“ªªØ
% for i = 1:NumUsers
%     sum_of_CCM = sum(sum(CCM(:,:,i)));
%     CCM(:,:,i) = CCM(:,:,i)./sum_of_CCM;
% end

end
    

