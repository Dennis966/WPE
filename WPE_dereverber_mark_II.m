clc
clear
close all

%% Weighted Prediction Error dereverber (WPE) Mark II

% ***********************************************************************************************
%% REQUIREMENT: 
% 1. Specify the path of the input & output audio files.
% 2. Specify the angle of target arrival.

%% CONTENT: 
% An algorithm that performs dereverberation via a linear prediction
% model for reverberant environments.

%% Important Parameters 
% n_iters = 1;                % Alternating Optimization iterations 
% epsilon = 10^(-5);           % Minimum energy for each stft value 
% delta = 10^(-10);            % Diagonal loading (numerical purposes)
% b = 8;                       % Prediction delay
% Lw_1 = 40;                   % order of AR model for 0 ~ 0.8kHz (1 ~ 14) (rounded) (the higher the better ???) 
% Lw_2 = 30;                   % order of AR model for 0.8k+ ~ 1.5kHz (15 ~ 25) (the higher the better ???)
% Lw_3 = 25;                   % order of AR model for 1.5k+ ~ 8kHz (26 ~ 129) (the higher the better ???)

%% Note 1: (4/13/2022) 
% If Lw_1, Lw_2, Lw_3 are too high, the algorithm will be extremely 
% slow. This is because the matrices will very large, and the
% pseudo-inverse will be computationally expensive! So you can increase
% these values but not too much.

%% Note 2: (4/13/2022) 
% delta only exists to ensure the condition number of the matrices
% are not too large (i.e., matrices are invertible).

%% Note 3: (4/13/2022) 
% The Alternating Optimization seem to converge very quickly.
% However, one can still try increasing n_iters. Maybe the dereverberation
% will be even better regardless of the fact that the algorithm already
% converged ?

%% Note 3: (4/13/2022) 
% The larger the value of epsilon, the larger the effect of the 
% dereverberation (not sure why), but we get a really large noise floor. 

%% Note 4: (4/13/2022) 
% I haven't tried the exact parameters as in the original paper.

%% Note 5: (4/14/2022) (My personal perceptual observations)
% For the current filter parameters, let's vary the n_iters
% n_iters = 1: medium dereverberation already !!
% n_iters = 2: slightly better.
% n_iters = 5: audibly better than n_iters = 1.                     
% n_iters = 10: maybe just a tiny bit better than n_iters = 5 (Really not obvious). 
% n_iters = 25: a bit unnecessary, audibly not much different than n_iters = 10.

% (For n_iters = [1, 25], waveforms seem to always looks better as n_iters 
% increases, but I cannot hear too much of a difference after n_iters = 5.)

%% Reference:
% T. Nakatani, T. Yoshioka, K. Kinoshita, M. Miyoshi and B. Juang, 
% "Speech Dereverberation Based on Variance-Normalized Delayed Linear 
% Prediction," in IEEE Transactions on Audio, Speech, and Language 
% Processing, vol. 18, no. 7, pp. 1717-1731, Sept. 2010, 
% doi: 10.1109/TASL.2010.2052251.
% ***********************************************************************************************

%% Initialize Default Parameters
input_dir = "./input/t_0_e_0.0000_T60_0point6/";
% input_dir = "./input/t_0_i_45_e_0.0000/";
t_i_e_specs = "target_0_error_0.0000";                   % target_interferer_error_specifications
% t_i_e_specs = "target_0_int_45_error_0.0000";

mic_1_pth = input_dir + "mic_1_" + t_i_e_specs + ".wav";
mic_2_pth = input_dir + "mic_2_" + t_i_e_specs + ".wav";
mic_3_pth = input_dir + "mic_3_" + t_i_e_specs + ".wav";
mic_4_pth = input_dir + "mic_4_" + t_i_e_specs + ".wav";
mic_5_pth = input_dir + "mic_5_" + t_i_e_specs + ".wav";
mic_6_pth = input_dir + "mic_6_" + t_i_e_specs + ".wav";
mic_7_pth = input_dir + "mic_7_" + t_i_e_specs + ".wav";
mic_8_pth = input_dir + "mic_8_" + t_i_e_specs + ".wav";
mic_9_pth = input_dir + "mic_9_" + t_i_e_specs + ".wav";
mic_10_pth = input_dir + "mic_10_" + t_i_e_specs + ".wav";

output_dir = "./output/t_0_e_0.0000_T60_0point6/";
% output_dir = "./output/t_0_i_45_e_0.0000/";

%% Read Wav Files
n_mics = 10;

[mic_1_wav, sr] = audioread(mic_1_pth);
[mic_2_wav, ~] = audioread(mic_2_pth);
[mic_3_wav, ~] = audioread(mic_3_pth);
[mic_4_wav, ~] = audioread(mic_4_pth);
[mic_5_wav, ~] = audioread(mic_5_pth);
[mic_6_wav, ~] = audioread(mic_6_pth);
[mic_7_wav, ~] = audioread(mic_7_pth);
[mic_8_wav, ~] = audioread(mic_8_pth);
[mic_9_wav, ~] = audioread(mic_9_pth);
[mic_10_wav, ~] = audioread(mic_10_pth);

%% STFT Conversion
% Parameter settings (NOTE: COLA needs to be satisfied)
win_length = 128;
win = hann(win_length, "periodic");                       % "periodic" to ensure COLA (Now, COLA seems to be a bit unnecessary?)
win_name = "hann";
n_fft = win_length*2;

%% obtain stft 
% MATLAB's stft(.)
overlap_perc = 75;                                                  % default
stft_1 = stft(mic_1_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_2 = stft(mic_2_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_3 = stft(mic_3_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_4 = stft(mic_4_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_5 = stft(mic_5_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_6 = stft(mic_6_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_7 = stft(mic_7_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_8 = stft(mic_8_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_9 = stft(mic_9_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 
stft_10 = stft(mic_10_wav, sr, "Window", win, "FFTLength", n_fft);    % 75% Overlap 

[stft_height, n_frames] = size(stft_1);
half_stfts = zeros(n_mics, n_fft - n_fft/2 + 1, n_frames);

half_stfts(1, :, :) = stft_1(n_fft/2:end, :);
half_stfts(2, :, :) = stft_2(n_fft/2:end, :);
half_stfts(3, :, :) = stft_3(n_fft/2:end, :);
half_stfts(4, :, :) = stft_4(n_fft/2:end, :);
half_stfts(5, :, :) = stft_5(n_fft/2:end, :);
half_stfts(6, :, :) = stft_6(n_fft/2:end, :);
half_stfts(7, :, :) = stft_7(n_fft/2:end, :);
half_stfts(8, :, :) = stft_8(n_fft/2:end, :);
half_stfts(9, :, :) = stft_9(n_fft/2:end, :);
half_stfts(10, :, :) = stft_10(n_fft/2:end, :);

d = 0.02;
c = 340;                                                      

%% Plot input signals 
% figure
% image(400*log10(abs(flip(half_stfts(1,:,:)))))
% xlabel("frames")
% title("Input spectrogram")
% 
% figure
% plot(10*log10(abs(half_stfts(1,:,50))))
% ylabel("magnitude (dB)")
% xlabel("normalized frequency (Hz)")
% title("Spectrum of frame 50")

%% Dereverb (WPE)
[half_y_stft, L_function, n_iters, epsilon, delta, b, Lw_1, Lw_2, Lw_3] = WPE(half_stfts, n_mics, n_fft, n_frames);

%% plot Likelihood function
%[FF, NN] = size(half_y_stft);
iter = 1:1:n_iters;
figure
for f_sections = 1:3
    hold on
    if f_sections == 1                    % plot mean likelihood function for low frequencies
        plot(iter, mean(L_function(:, 1:14), 2))
    elseif f_sections == 2                % plot mean likelihood function for mid frequencies
        plot(iter, mean(L_function(:, 15:25), 2))
    else                                  % plot mean of likelihood function for mid frequencies
        plot(iter, mean(L_function(:, 26:129), 2))
    end
    hold off
end
grid on
legend(["low", "mid", "high"])
xlabel("Iterations")
ylabel("$\mathcal{L}$", Interpreter="latex")
xticks(iter)
title("Likelihood Function Score", Interpreter="latex")

%% ISTFT conversion

upper_half_y_stft = conj(flipud(half_y_stft(2:n_fft/2, :)));
y_stft = [upper_half_y_stft; half_y_stft];
y = istft(y_stft, sr, "Window", win, "FFTLength", n_fft);       % 75% Overlap 
y = real(y);
y = y/max(abs(y));

%% Plot results (wave files before and after beamforming)
ref = mic_1_wav/max(abs(mic_1_wav));
figure
hold on
plot(ref)
plot(y)
xlim([1, sr*3])       
xlabel("n")
title("Dereverbed Results")
legend("original", "Dereverbed audio")
hold off
savefig("sr_" + num2str(sr) ...
      + "_win_length_" + num2str(win_length) ...
      + "_Overlap_perc_" + num2str(overlap_perc) ...
      + "_win_" + win_name ...
      + "_n_mics_" + num2str(n_mics) ...
      + "_n_fft_" + num2str(n_fft) ...
      + "_n_iters_" + num2str(n_iters) ...
      + "_epsilon_" + num2str(epsilon) ... 
      + "_delta_" + num2str(delta) ...
      + "_b_" + num2str(b) ... 
      + "_Lw_1_" + num2str(Lw_1)...
      + "_Lw_2_" + num2str(Lw_2) ...
      + "_Lw_3_" + num2str(Lw_3))

%% Store Output Wave File
% audiowrite(output_dir + t_i_e_specs + "_original.wav", ref, sr)
audiowrite(output_dir + t_i_e_specs + "_WPE" ...
                                    + "_sr_" + num2str(sr) ...
                                    + "_win_length_" + num2str(win_length) ...
                                    + "_Overlap_perc_" + num2str(overlap_perc) ...
                                    + "_win_" + win_name ...
                                    + "_n_mics_" + num2str(n_mics) ...
                                    + "_n_fft_" + num2str(n_fft) ...
                                    + "_n_iters_" + num2str(n_iters) ...
                                    + "_epsilon_" + num2str(epsilon) ... 
                                    + "_delta_" + num2str(delta) ...
                                    + "_b_" + num2str(b) ... 
                                    + "_Lw_1_" + num2str(Lw_1)...
                                    + "_Lw_2_" + num2str(Lw_2) ...
                                    + "_Lw_3_" + num2str(Lw_3) + ".wav", y, sr)


