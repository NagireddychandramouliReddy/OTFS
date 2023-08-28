% Parameters
num_symbols = 100; % Number of symbols

% Generate random bits (0s and 1s)
bits = randi([0 1], num_symbols * 2, 1);

% Reshape bits into symbols
symbols = reshape(bits, 2, num_symbols).';

% Map bits to QPSK symbols
qpsk_symbols = 1/sqrt(2) * (1 - 2*symbols(:, 1)) + 1i/sqrt(2) * (1 - 2*symbols(:, 2));

disp('Generated QPSK symbols:');
plot(qpsk_symbols);
title ("QPSK signal")
