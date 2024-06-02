clc;clear;close;
load 2024-06-02.mat
addpath ..\
%Theoretical OFDM Curve Plotting
ber_Rayleigh = ber_QAM(EbN0dB,M,'Rayleigh');
semilogy(EbN0dB,ber_Rayleigh,'r-');
hold on
%OFDM_IM
semilogy(EbN0dB,BER,'r--o'); 
hold on
legend({'OFDM\_Theoretical','OFDM\_IM'})
title("OFDM\_IM")
xlabel("Eb/N0")
ylabel("BER")
