clc
clear;
close all;
%% Estimation Params:

ENABLE_CFO = 1;                      % flag to enable car freq offset
DETECTION_OFFSET = 100;              % to add packet detection error 
SNR = 10; 
snr_lin = 10^(SNR/20);
TRIGGER_OFFSET_TOL_NS   = 3000;      % Trigger time offset toleration between Tx and Rx that can be accomodated

%% OFDM params
PILOT_SUB_CAR_IND = [8 22 44 58];                           % Pilot subcarrier indices
DATA_SUB_CAR_IND = [2:7 9:21 23:27 39:43 45:57 59:64];      % Data subcarrier indices
NUM_SUB_CAR = 64;                                           % Number of subcarriers
CP_LEN = 16;                                                % Cyclic prefix length
channel_coding = .5;                                        %  R = k/n k=1,n=2
trel_len = 8;                                               % bits for trellis to end(horizontal) 
NUM_OFDM_SYMS = 500;                                        % Number of OFDM symbols
MOD_ORDER = 4;                                              % {1,2,4,6} = {BSPK, QPSK, 16-QAM, 64-QAM}
NUM_DATA_SYMS = NUM_OFDM_SYMS * length(DATA_SUB_CAR_IND);   % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
	

%% Preamble

% STS
short_ts_f = zeros(1,64);
short_ts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
short_ts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
figure
plot(abs(short_ts_f));
xlabel('sub-car index');
ylabel('magnitude');
title('Short Training Symbol');
short_ts_t = ifft(sqrt(13/6).*short_ts_f, 64);  % scalar is for normalizing the power
short_ts_t = short_ts_t(1:16); 

% LTS
long_ts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
long_ts_t = ifft(long_ts_f, 64);
figure
plot(abs(long_ts_f));
xlabel('sub car index');
ylabel('magnitude');
title('Long Training Symbol');
ylim([-1,2]);

preamble = [repmat(short_ts_t, 1, 10)  long_ts_t(33:64) long_ts_t long_ts_t];
figure
plot(abs(preamble));
xlabel('sample index');
ylabel('magnitude');
title('Preamble time domain');

%% Generate a payload of random integers
num_tx_bits= (NUM_DATA_SYMS * MOD_ORDER - 2*trel_len) * channel_coding;
tx_data = randi(2, 1, num_tx_bits) - 1; % converts 1,2 into 0,1

%% Forward Error Correction 
tx_data = double([tx_data zeros(1,trel_len) ]);    % 8 bits padding
trel = poly2trellis(7, [141 133]);      % Define trellis with constraint length 7 and octal polnomial
tx_code = convenc(tx_data,trel);        % convultional encoder 

%% modulation
tx_vec = modulator(tx_code', MOD_ORDER, 1);
figure
scatter(real(tx_vec), imag(tx_vec),'filled');
title('Constellation of transmitted bits');
xlabel('In-phase'); ylabel('Quadrature-phase');
% reshape gives a matrix
tx_syms_mat = reshape(tx_vec, length(DATA_SUB_CAR_IND), NUM_OFDM_SYMS); %one column per OFDM symbol
% pilot symbols always BPSK mod
BPSK_pilots = [1 1 -1 1].';
% append the pilots to symbols
pilots_mat = repmat(BPSK_pilots, 1, NUM_OFDM_SYMS);


%% IFFT
ifft_input_mat = zeros(NUM_SUB_CAR, NUM_OFDM_SYMS);
% insert data and pilot values; force other sub car to 0
ifft_input_mat(DATA_SUB_CAR_IND, :)   = tx_syms_mat;
ifft_input_mat(PILOT_SUB_CAR_IND, :) = pilots_mat;
tx_data_mat_t = ifft(ifft_input_mat, NUM_SUB_CAR, 1); %ifft across dimension 1(each column)
%% append  cyclic prefix
if(CP_LEN > 0)
    tx_cp_t = tx_data_mat_t((end-CP_LEN+1 : end), :);
    tx_data_mat_t = [tx_cp_t; tx_data_mat_t];
end

% convert back to vector
tx_payload_vec = reshape(tx_data_mat_t, 1, numel(tx_data_mat_t));
% time-domain OFDM 
tx_vec_t = [preamble tx_payload_vec];


%% Interpolation
% DAC frequency response will be a sinc function because of its sample and
% hold nature
% Define a half-band 2x interpolation filter response (sinc function)
interp_filt = sinc(-10.5:0.5:10);
% put zeros at alternate indices to give a vec of twice the original length
% apply low pass filter(sinc) to properly interpolate
tx_vec_zero_pad = [tx_vec_t, zeros(1, ceil(length(interp_filt)/2))];
tx_vec_2x_interp = zeros(1, 2*numel(tx_vec_zero_pad)); %numel is mXn for a matrix
tx_vec_2x_interp(1:2:end) = tx_vec_zero_pad;
tx_vec_dac = filter(interp_filt, 1, tx_vec_2x_interp); %after conversion from digital into analog

figure
plot(db(abs(fftshift(fft(tx_vec_dac)))));
title('frequency domain tx vector');
xlabel('frequency');
xlim([8000,80000]); ylim([0,65]);
figure
subplot(2,1,1);
plot(real(tx_vec_dac));
title('real part of tx vector');
subplot(2,1,2);
plot(imag(tx_vec_dac));
title('imaginary part of tx vector');
xlabel('index n');

%% simulating the wireless channel.
if(ENABLE_CFO)
    tx_vec_dac = tx_vec_dac .* exp(-1i*2*pi*1e-4*[0:length(tx_vec_dac)-1]);
end

% AWGN:
SAMP_FREQ = 40e6;        % Sampling frequency
rx_vec_air = [zeros(1,80),tx_vec_dac, zeros(1,ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)))-80]; % simulating delay
rx_sig_power = sum(abs(rx_vec_air).^2);
%rx_vec_air = rx_vec_air  +  1*0.009*complex(randn(1,length(rx_vec_air)), randn(1,length(rx_vec_air)));

% Decimate
raw_rx_data = filter(interp_filt, 1, rx_vec_air);
raw_rx_data = [zeros(1,DETECTION_OFFSET) raw_rx_data(1:2:end)];


%% %%%%%%%%%%%%% RECEIVER %%%%%%%%%%%%%


%% Rx processing params
mod_type = 4;
LTS_CORR_THRESH = 0.8;         % Normalized threshold for LTS correlation

%% Packet Detection, Time Estimation
lts_xcorr = abs(xcorr(long_ts_t,sign(raw_rx_data)));
lts_xcorr = fliplr(lts_xcorr(32:(length(raw_rx_data)+32)));
max_lts_xcorr =  max(lts_xcorr)*0.8;
lts_xcorr_peaks = find(lts_xcorr(1:end) > max_lts_xcorr);
payload_ind = max(lts_xcorr_peaks)+32; % The "+32" corresponds to the 32-sample cyclic prefix on the preamble LTS
lts_index = payload_ind-160; % The "-160" corresponds to the length of the preamble LTS (2.5 copies of 64-sample LTS)
figure
plot(lts_xcorr);
xlim([1 2000]);
title('cross correlation of LTS symbols');

Rxy = abs(xcorr(raw_rx_data(1:length(preamble)),preamble)); %becomes auto corr in absence of noise
figure
plot(Rxy);
title('auto correlation of rx vector');
xlabel('samples');
ylabel('Rxy');


%% CFO estimation and correction
SAMP_FREQ_OFFSET = 4; % Number of CP samples to use in FFT(FFT OFFSET)
N=64;
rx_lts_t = raw_rx_data(lts_index : lts_index+159);
rx_lts_1_t = rx_lts_t(-64+-SAMP_FREQ_OFFSET + [97:160]);
rx_lts_2_t = rx_lts_t(-SAMP_FREQ_OFFSET + [97:160]);
rx_cfo_est_lts = mean(unwrap(angle(rx_lts_2_t .* conj(rx_lts_1_t))));
rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*N);    
rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[1:length(raw_rx_data)]);
rx_cfo_corrected = raw_rx_data .* rx_cfo_corr_t;


%% FFT and CP Removal
% payload is data
payload_data_rx = rx_cfo_corrected(payload_ind : payload_ind+NUM_OFDM_SYMS*(NUM_SUB_CAR+CP_LEN)-1); % obtaining data from CFO corrected vec
payload_data_matrix_CP = reshape(payload_data_rx, (NUM_SUB_CAR+CP_LEN), NUM_OFDM_SYMS);
% Remove the cyclic prefix, keep FFT_OFFSET samples of CP intact
payload_mat_noCP = payload_data_matrix_CP(CP_LEN - SAMP_FREQ_OFFSET + [1:NUM_SUB_CAR], :);% gives matrix of size N_SC X N_OFDM_SYMS
syms_mat_f = fft(payload_mat_noCP, NUM_SUB_CAR, 1); % 64 point fft of the cols

%% Channel estimation and correction
% two copies of LTS to estimate the channel
rx_lts_t = rx_cfo_corrected(lts_index : lts_index+159); % #160
rx_lts_1_t = raw_rx_data(lts_index+96 : lts_index+159); % #64
rx_lts_2_t = raw_rx_data(lts_index+32 : lts_index+95);  % #64
rx_lts1_f = fft(rx_lts_1_t); 
rx_lts2_f = fft(rx_lts_2_t);
rx_lts_f = fft(rx_lts_t(-SAMP_FREQ_OFFSET + [97:160])); 
rx_H_est = long_ts_f .* (rx_lts_f);
channel_matrix = repmat(transpose(rx_H_est),1,NUM_OFDM_SYMS);
syms_eq_mat = syms_mat_f ./ channel_matrix; 
x = (20/NUM_SUB_CAR) * (-(NUM_SUB_CAR/2):(NUM_SUB_CAR/2 - 1));

figure
stairs(x - (20/(2*NUM_SUB_CAR)), fftshift(real(rx_H_est)));
hold on
stairs(x - (20/(2*NUM_SUB_CAR)), fftshift(imag(rx_H_est)));
hold off
axis([min(x) max(x) -1.1*max(abs(rx_H_est)) 1.1*max(abs(rx_H_est))])
grid on;
xlabel('subcarrier index(DC at 0)');
title('Channel Estimates Plot (I and Q)');
legend('I-phase','Q-phase');

%% SFO estimation and correction using pilots
pilots_mat_f = syms_eq_mat(PILOT_SUB_CAR_IND, :); %obtain pilots from the matrix
pilots_mat_f_comp = pilots_mat_f.*pilots_mat; %equalizing pilots by their tx values
% phase of every Rx pilot symbol
pilot_phase_angles = unwrap(angle(fftshift(pilots_mat_f_comp,1)), [], 1);
% slope of pilot tone phases wrt frequency
pilot_spacing_mat = repmat(mod(diff(fftshift(PILOT_SUB_CAR_IND)),64).', 1, NUM_OFDM_SYMS);                        
pilot_slope_mat = mean(diff(pilot_phase_angles) ./ pilot_spacing_mat); %mean of derivatives

%% Calculate the SFO correction phases for each OFDM symbol
 pilot_phase_sfo_corr = fftshift((-32:31)' .* pilot_slope_mat, 1);
pilot_phase_corr = exp(-1i*(pilot_phase_sfo_corr));
syms_eq_mat = syms_eq_mat .* pilot_phase_corr;% pilot phase corrected equalized mat 

%% Phase Error Correction using pilots
% Extract the pilots and calculate per-symbol phase error
pilots_mat_f = syms_eq_mat(PILOT_SUB_CAR_IND, :);
pilots_mat_f_comp = pilots_mat_f.*pilots_mat;
pilot_phase_err = angle(mean(pilots_mat_f_comp));  
pilot_phase_err_corr = repmat(pilot_phase_err, NUM_SUB_CAR, 1);
pilot_phase_corr = exp(-1i*(pilot_phase_err_corr));
% Apply the pilot phase correction per symbol
sym_pc_corrected_mat = syms_eq_mat .* pilot_phase_corr;
payload_data_sym_mat = sym_pc_corrected_mat(DATA_SUB_CAR_IND, :);
final_sym_vec = reshape(payload_data_sym_mat, 1, NUM_DATA_SYMS);

%% Demodulation and  viterbi decoder
figure
scatter(real(final_sym_vec), imag(final_sym_vec),'filled');
title(' Received Constellation');
xlabel('In-phase'); ylabel('Qudrature-phase');
Demod_out = demodulator(final_sym_vec,mod_type,1);
rx_data_viterbi_dec = vitdec(Demod_out,trel,7,'trunc','hard'); %truncated operating mode and hard decoding
[num_err_bits,ber] = biterr(tx_data,rx_data_viterbi_dec);
disp(ber);
disp(num_err_bits)
