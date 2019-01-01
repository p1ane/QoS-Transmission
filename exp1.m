clear all;
clc;

addpath('.\data');
[Nr,Nt,NumUsers,~,NumLinks,H_freq,H_freq_beam,CCM] = gen_of_channel;
[Hev_freq,Hev_freq_beam,CCMev] = gen_of_evchannel;

load('CCM1.mat');
load('CCMev1.mat');

NumSamples = 5000;

[H_freq_beam_test,CCM_test] = reform_H_freq_beam(CCM,Nt,Nr,NumLinks,NumSamples);
[Hev_freq_beam_test,CCMev_test] = reform_Hev_freq_beam(CCMev,Nt,Nr,NumSamples);
