% Parameters
num_symbols = 100; % Number of symbols
num_subcarriers = 64; % Number of subcarriers
SNR_dB = 20; % Signal-to-Noise Ratio in dB

% Generate random QPSK symbols for each transmit antenna
tx_symbols_antenna1 = qpskmod(randi([0 3], num_symbols, 1));
tx_symbols_antenna2 = qpskmod(randi([0 3], num_symbols, 1));

% Combine the symbols for both antennas
tx_symbols = [tx_symbols_antenna1, tx_symbols_antenna2];

% OTFS modulation
modulated_symbols = otfs_modulate(tx_symbols, num_subcarriers);

% Channel (you can modify this based on your scenario)
H = sqrt(0.5) * (randn(num_subcarriers, 2, 2) + 1i * randn(num_subcarriers, 2, 2));

% Transmit through the channel
received_symbols = modulated_symbols * H;

% Add AWGN noise
SNR = 10^(SNR_dB / 10);
noise_var = 1 / (2 * SNR); % Two antennas
noise = sqrt(noise_var) * (randn(size(received_symbols)) + 1i * randn(size(received_symbols)));
received_symbols_noisy = received_symbols + noise;

% OTFS demodulation
demodulated_symbols = otfs_demodulate(received_symbols_noisy, num_subcarriers);

% QPSK demodulation
rx_symbols_antenna1 = qpskdemod(demodulated_symbols(:, 1));
rx_symbols_antenna2 = qpskdemod(demodulated_symbols(:, 2));

% Calculate symbol error rate
SER_antenna1 = sum(abs(rx_symbols_antenna1 - tx_symbols_antenna1) > 0) / num_symbols;
SER_antenna2 = sum(abs(rx_symbols_antenna2 - tx_symbols_antenna2) > 0) / num_symbols;

disp(['Symbol Error Rate (Antenna 1): ', num2str(SER_antenna1)]);
disp(['Symbol Error Rate (Antenna 2): ', num2str(SER_antenna2)]);
