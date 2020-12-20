% CPFM-Full Response Rectangular Filter using Weakly Orthogonal Signals
% M = 2 (Minimum shift keying); h (modulation index) = 1/2
clear all
close all
clc
EbNo_dB = 9; % Eb/No in (dB)
num_bit = 1000; % number of data bits in each frame
sym_dur = 1; % symbol duration (seconds)
mod_index = 1/2; % modulation index
band_width = 100; % bandwidth (Hz)
fs = 4*band_width; % sampling frequency
decoding_delay = 20; % decoding delay of the Viterbi Algorithm
sim_runs = 10^0; % simulation runs

% SNR parameters
EbNo = 10^(0.1*EbNo_dB);
noise_var_1D = sym_dur*fs/(2*EbNo); % awgn variance

% time vector
time_vec = 0:1/fs:num_bit*sym_dur;
time_vec(end)=[];
reshape_time_vec = reshape(time_vec,sym_dur*fs,num_bit);
%--------------------------------------------------------------------------
% Precomputations required for branch metrics for VA
temp1=(1/sqrt(sym_dur))*exp(1i*(mod_index*pi/sym_dur)*(1)*time_vec(1:4*sym_dur*fs));
temp2 =(1/sqrt(sym_dur))*exp(1i*(mod_index*pi/sym_dur)*(-1)*time_vec(1:4*sym_dur*fs));
%              matched filter
mf_sig1 = conj(fliplr(temp1(1:sym_dur*fs)));
mf_sig2 = conj(fliplr(temp1(sym_dur*fs+1:2*sym_dur*fs)));
mf_sig3 = conj(fliplr(temp1(2*sym_dur*fs+1:3*sym_dur*fs)));
mf_sig4 = conj(fliplr(temp1(3*sym_dur*fs+1:4*sym_dur*fs)));
mf_sig5 = conj(fliplr(temp2(1:sym_dur*fs)));
mf_sig6 = conj(fliplr(temp2(sym_dur*fs+1:2*sym_dur*fs)));
mf_sig7 = conj(fliplr(temp2(2*sym_dur*fs+1:3*sym_dur*fs)));
mf_sig8 = conj(fliplr(temp2(3*sym_dur*fs+1:4*sym_dur*fs)));
%--------------------------------------------------------------------------
C_Ber = 0; % total number of bit errors
tic()
for i1 = 1:sim_runs
%--------------------------------------------------------------------------
%                       Transmitter
% source
a = randi([0 1],1,num_bit);

% 2-PAM mapping
pam_sig = 1-2*a;

% phase modulation
init_phase = 0; % initial phase
phase_mod_sig = zeros(sym_dur*fs,num_bit); % initialization
for i2 = 1:num_bit
 phase_mod_sig(:,i2)=init_phase+(mod_index*pi/sym_dur)*pam_sig(i2)*(reshape_time_vec(:,i2)-(i2-1)*sym_dur);
 init_phase =  phase_mod_sig(end,i2);  
end
phase_mod_sig = transpose(phase_mod_sig(:)); % now a row vector

% complex low-pass equivalent signal
trans_sig = (1/sqrt(sym_dur))*exp(1i*phase_mod_sig);
%--------------------------------------------------------------------------
%                           channel
awgn =normrnd(0,sqrt(noise_var_1D),1,num_bit*sym_dur*fs)+1i*normrnd(0,sqrt(noise_var_1D),1,num_bit*sym_dur*fs);
chan_op = trans_sig + awgn; % channel output
%--------------------------------------------------------------------------
%                            Receiver
reshape_chan_op = reshape(chan_op,sym_dur*fs,num_bit);
chan_op1 = reshape_chan_op.'; % transpose
mf_op1 = real(conv2(chan_op1,mf_sig1));
mf_op2 = real(conv2(chan_op1,mf_sig2));
mf_op3 = real(conv2(chan_op1,mf_sig3));
mf_op4 = real(conv2(chan_op1,mf_sig4));
mf_op5 = real(conv2(chan_op1,mf_sig5));
mf_op6 = real(conv2(chan_op1,mf_sig6));
mf_op7 = real(conv2(chan_op1,mf_sig7));
mf_op8 = real(conv2(chan_op1,mf_sig8));

% Branch metrics for the Viterbi algorithm
branch_metric =  zeros(8,num_bit);
branch_metric(1,:) = mf_op1(:,sym_dur*fs).';
branch_metric(2,:) = mf_op2(:,sym_dur*fs).';
branch_metric(3,:) = mf_op3(:,sym_dur*fs).';
branch_metric(4,:) = mf_op4(:,sym_dur*fs).';
branch_metric(5,:) = mf_op5(:,sym_dur*fs).';
branch_metric(6,:) = mf_op6(:,sym_dur*fs).';
branch_metric(7,:) = mf_op7(:,sym_dur*fs).';
branch_metric(8,:) = mf_op8(:,sym_dur*fs).';

dec_a = Viterbi_alg(branch_metric,num_bit,20);

C_Ber = C_Ber + nnz(dec_a - a(1:num_bit-decoding_delay));
end
toc()
%bit error rate
BER = C_Ber/(sim_runs*(num_bit-decoding_delay)) % ignoring last transient samples
