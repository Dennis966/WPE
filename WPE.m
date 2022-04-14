function [half_y_stft, L_function, n_iters, epsilon, delta, b, Lw_1, Lw_2, Lw_3] = WPE(half_stfts, n_mics, n_fft, n_frames)
% half_y_stft: upper half of output

% Hyperparameters
% epsilon = 10^(-7);
n_iters = 5;  
epsilon = 10^(-5);
delta = 10^(-10);
b = 8;                                             % prediction delay
Lw_1 = 40;                                         % order of AR model for 0 ~ 0.8kHz (1 ~ 14) (rounded)
Lw_2 = 30;                                         % order of AR model for 0.8k+ ~ 1.5kHz (15 ~ 25)
Lw_3 = 25;                                         % order of AR model for 1.5k+ ~ 8kHz (26 ~ 129)

len_1 = n_mics * (Lw_1);
len_2 = n_mics * (Lw_2);
len_3 = n_mics * (Lw_3);

% 0 ~ 0.8kHz (1 ~ 14) (rounded)
%Rxx_1 = zeros(len_1, len_1, 14); phi_1 = zeros(len_1, 14); c_1 = zeros(len_1, 14);        

% 0.8k+ ~ 1.5kHz (15 ~ 25)
%Rxx_2 = zeros(len_2, len_2, 11); phi_2 = zeros(len_2, 11); c_2 = zeros(len_2, 11);                    

% 1.5k+ ~ 8kHz (26 ~ 129)
%Rxx_3 = zeros(len_3, len_3, 104); phi_3 = zeros(len_3, 104); c_3 = zeros(len_3, 104);               

% initialize estimates
half_y_stft = zeros(n_fft/2 + 1, n_frames);                   % Let's first optimize channel 1                                             % around 25, likelihood function converges
L_function = zeros(n_iters, n_fft/2 + 1);

%% Alternating Optimization
fprintf("************** WPE Alternating Optimization **************\n");
tic
for iter = 1:n_iters
    fprintf("Iteration = %d \n", iter);
    %% Optimize rhos
    if iter == 1
        rhos = max(abs( squeeze( half_stfts(1, :, :) ) ).^2, epsilon);
    else
        rhos = max(abs(half_y_stft).^2, epsilon);
    end

    %% Optimize c_bars

    for f_bin = 1:n_fft/2 + 1
        %% Parameter settings 
        if (1 <= f_bin) && (f_bin <= 14)       % 0 ~ 0.8kHz (14 subbands)
            Lw = Lw_1; len = len_1; 
        elseif (15 <= f_bin) && (f_bin <= 25)  % 0.8k+ ~ 1.5kHz (11 subbands)
            Lw = Lw_2; len = len_2; 
        else                                   % 1.5k+ ~ 8kHz (26 ~ 129)
            Lw = Lw_3; len = len_3; 
        end

        Rxx = zeros(len, len); phi = zeros(len, 1);
        for frame = 1:n_frames
            x_bar = x_bar_gen(half_stfts, len, Lw, f_bin, frame, b, n_mics);
            Rxx = Rxx + (x_bar * x_bar') / rhos(f_bin, frame);
            phi = phi + (x_bar * conj(half_stfts(1, f_bin, frame))) / rhos(f_bin, frame);    
        end
        c = (Rxx + delta*eye(len)) \ phi;
        
        for frame = 1:n_frames
            x_bar = x_bar_gen(half_stfts, len, Lw, f_bin, frame, b, n_mics);
            half_y_stft(f_bin, frame) = half_stfts(1, f_bin, frame) - c'*x_bar;

            % calculate likelihood function
            L_function(iter, f_bin) = L_function(iter, f_bin) - (abs(half_y_stft(f_bin, frame))^2/rhos(f_bin, frame) + log(rhos(f_bin, frame)));
        end
    end
    fprintf("Likelihood Score (Low freq) = %d \n", mean( L_function(iter, 1:14), 2) );              % no log() yet
    fprintf("Likelihood Score (Mid freq) = %d \n", mean( L_function(iter, 15:25), 2) );             % no log() yet
    fprintf("Likelihood Score (High freq) = %d \n\n", mean( L_function(iter, 26:129), 2) );         % no log() yet
end
toc
fprintf("************** End **************\n");
end