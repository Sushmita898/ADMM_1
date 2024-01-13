

clc
clear all
close all


%% OTFS parameters%%%%%%%%%%
% number of symbol
N = 8;
% number of subcarriers
M = 8;
% size of constellation
taps=3;
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;
n_irs = 4;

 SNR_dB = 0:5:20;
%    SNR_dB =10;
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
noise_var = 1./SNR;
% noise_var = abs((noise_var_sqrt).^2)/eng_sqrt;
 eng_sqrt = 1;
  sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;
 % sigma_2 = abs(noise_var_sqrt).^2;
%% ADMM Para0meters

 eps_par = 0.005;
alpha=1;
%%
rng(1)
N_fram = 10000;
n_ite = 40;
err_ber = zeros(length(SNR_dB),1);
  max_error=[ 1000, 800,100,50,10];

% max_error=[ 10];
theta = zeros(1,n_irs);
phaseFactor = zeros(1,n_irs);

theta_list = [0,pi*2/8,pi*4/8,pi*6/8,pi*8/8,pi*10/8,pi*12/8,pi*14/8];
for iesn0 = 1:length(SNR_dB)
%  alpha = alpha*sigma_2(iesn0);
%   rho = rho.^2
   rho = eps_par*sigma_2(iesn0);
%  rho=eps_par;
 ifram = 0;
   while err_ber(iesn0)<max_error(iesn0)
% for ifram = 1:N_fram
        %% random input bits generation%%%%%
        data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));        
        x = qammod(data_temp,M_mod);
        x = reshape(x,N,M);
       
        x = reshape(x,N*M,1);
     
        %% OTFS modulation%%%%
        s = OTFS_modulation(N,M,x);

         %% OTFS channel generation%%%%


        R = zeros(N*M,1);
        H = [];
       for k = 1:n_irs
%        theta(k) = theta_list(randi(length(theta_list)));%从theta list里面随机取一个
%     phaseFactor(k) = cos(theta(k))+1i*sin(theta(k));
  %% generating OTFS-IRS Channel parameter
      taps=3;
        delay_taps1=0:taps-1;
        Doppler_taps1=0:taps-1;
       %[chan_coef1] = OTFS_channel_gen(N,M,taps);    % BS to IRS channel
       pow_prof = (1/taps) * (ones(1,taps));
chan_coef1 = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
        delay_taps2=0:taps-1;
        Doppler_taps2=0:taps-1;
        %[chan_coef2] = OTFS_channel_gen(N,M,taps);   % IRS to User channel
pow_prof = (1/taps) * (ones(1,taps));
chan_coef2 = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
        [delaytaps_irs,dopplertaps_irs,chan_coef_irs] = generate_channel_for_irs(taps,delay_taps1,delay_taps2,Doppler_taps1,Doppler_taps2,chan_coef1,chan_coef2);  % Cascaded channel from BS to User via IRS
       
        Taps = length(delaytaps_irs);
 chan_coef =  (randn(1,Taps)+1i*randn(1,Taps));
sum_ch_path_coeff = zeros(1,Taps);
%h_bkp=zeros(B,K,P,BS_antennas,User_antennas);
% for b=1:B
%     for k=1:K
%         for p=1:P
for t = 1:Taps      %% Number of channel path
    %for r = 1:N_ris*R   %% Number of IRS element
  sum_ch_path_coeff(t) = sum_ch_path_coeff(t)+ abs(chan_coef_irs(1,t)) ;  %% summing the channel coefficient of all the irs element for each taps

    %end
end
 sum_ch_path_coeff = sum_ch_path_coeff +abs(chan_coef);
%% finding the channel path with maximum power
t_hat = find(max(sum_ch_path_coeff));  %% Tap which gives maximum power

%% Obtaining optimum angle
% theta_hat = zeros(N_ris*R,1);
% for r = 1:N_ris*R
   theta_hat(k) =   angle(chan_coef(t_hat)) - angle(chan_coef_irs(1,t_hat));
% end
% Theta = diag(theta_hat);
    phaseFactor(k) = cos(theta_hat(k))+1i*sin(theta_hat(k));     
     theta = rand(1,n_irs)+1j*rand(1,n_irs) ;  
        
        %% Generating effective channel matrix for cascaded channel
        H_rect_d=find_effective_H_rect_new(delay_taps1,Doppler_taps1,chan_coef,M,N);
      H_rect=find_effective_H_rect_new(delaytaps_irs,dopplertaps_irs,(chan_coef_irs),M,N);
%      
%          H_rect=find_effective_H_rect_new(delaytaps_irs,dopplertaps_irs,sum_ch_path_coeff,M,N);
     H_rect = H_rect_d+H_rect;
    
        H = [H H_rect];  %% cascaded channel of all n_irs IRS element
     end
       Theta = kron(phaseFactor,eye(N*M));
%  Theta = kron(theta,eye(N*M));
%  Theta= kron(ones(1,n_irs),eye(N*M));
        H_rect1 = H*Theta';  
%          H_rect1 = H_rect_d+H_rect1;
       H = H_rect1;
   % noise = sqrt(sigma_2(iesn0)/2)*(randn(size(s)) + 1i*randn(size(s)));
    noise = 1/sqrt(2)*(sigma_2(iesn0))*(randn(size(s)) + 1i*randn(size(s)));
%      noise = sqrt(noise_var(iesn0)/2)*(randn(size(s)) + 1i*randn(size(s)));
     x = reshape(x,N*M,1);
       y = H*x+noise;
 

  x_est =OTFS_ADMM_detector(H,N,M,M_mod,y, rho, alpha,n_ite);
       
        
        %% output bits and errors count%%%%%
        data_demapping = qamdemod(x_est,M_mod);
        data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
        errors = sum(xor(data_info_est,data_info_bit));
        err_ber(iesn0) = errors + err_ber(iesn0)
       % fprintf('\n finished for %d frame', ifram);
        ifram = ifram+1;
         end
    fprintf('\n finished for %d SNR', SNR_dB(iesn0));
    %end
    
end
err_ber_fram = err_ber/N_bits_perfram./ifram
semilogy(SNR_dB, err_ber_fram,'-*','LineWidth',2);
title(sprintf('OTFS-IRS STM'))
ylabel('BER'); xlabel('SNR in dB');grid on

hold on