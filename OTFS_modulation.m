function s = OTFS_modulation(N,M,x)
%% OTFS Modulation: 1. ISFFT, 2. Heisenberg transform
X = fft(ifft(x).').'/sqrt(M/N); %%%ISFFT
s_mat = ifft(X.')*sqrt(M); % Heisenberg transform
s = s_mat(:);
end