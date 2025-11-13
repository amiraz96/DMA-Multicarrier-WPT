clear all
close all
clc
%% Optimization

phi_lorentz = linspace(0, 1, 10e3).*2.*pi;
lorentz_vec = (1i + exp(1i.*phi_lorentz))./2;
it_main_max = 15;
p = 1000e-3;
G = 1; gamma = 1;
As = 1;
M_set = [1 2];
L_set = [0.25 0.25];
N_set = [1 2 4 8];

x_min = 0; x_max = 5; y_min = 0; y_max = 5; z_min = 0; z_max = 3;

max_it = 1;

%% General parameters

f_0 = 5.18e9; %% Operating Frequency
boresight_gain = 2;
BW = 10e6;
DC_Threshold_main = 0.2e-4; 
eps_thr = 1e-3;

%% Rectifier Parameters
i_s = 5*10^(-6); eta_0 = 1.05; V_o = 25.86 * 10^(-3); R_ant = 50; R_L = 10e3;  

%% Sampling
TT = 1e-3;
nfac = 1;
samp_num = 10;
load('user_location.mat');
user_iters = size(user_loc1, 2);

P_Cons_Results_DMA_overN = zeros(length(L_set), length(M_set), length(N_set), user_iters);
Pin_Cons_Results_DMA_overN = zeros(length(L_set), length(M_set), length(N_set), user_iters);
PA_Cons_Results_DMA_overN = zeros(length(L_set), length(M_set), length(N_set), user_iters);

start_L = 1; start_M = 1; start_N = 1; 



for ll = start_L:size(L_set, 1)
    antenna_length_vec = L_set(ll, :);
    % DMA Location
    [DMA_element_loc, H, RFC_num, passive_num] = ...
        DMA_deploy(f_0, antenna_length_vec, x_max, y_max, z_max);
    for mm = start_M:length(M_set)
        user_num = M_set(mm);
        if user_num == 1
            user_locc = user_loc1;
        elseif user_num == 2
            user_locc = user_loc2;
        else
            user_locc = user_loc4;
        end
        for nn = start_N:length(N_set)
            sub_number = N_set(nn);
            delta_f = BW/sub_number;
            freq_mat = zeros(1, sub_number); %% Sub-Carriers
            for ii = 1:sub_number 
                freq_mat(ii) = f_0 + (ii - 1)*delta_f;
            end
            samp_rate = 2*max(freq_mat); % Nyquist
            samp_int = 1/(nfac*samp_rate);
            samp_vec = 0:samp_int:samp_num*samp_int; % Sample vector
            Nt = length(samp_vec);
            Q_dma = zeros(RFC_num, passive_num);
            H_dma = zeros(RFC_num, passive_num);
            phi_vec = rand(RFC_num, passive_num);
            Q = zeros(passive_num*RFC_num, RFC_num);
            for i = 1:RFC_num
                for l = 1:passive_num
                    idx = (i-1)*passive_num + l;
                    H_dma(i,l) = H(idx,idx);
                    Q((i - 1)*passive_num + l, i) = Q_DMA(phi_vec(i,l) * 2 * pi);
                    Q_dma(i, l) = Q_DMA(phi_vec(i,l) * 2 * pi);
                end
            end
            for user_it = 1:10
                user_loc = reshape(user_locc(:, user_it, :), [user_num, 3]);
                [channel_vec] = Channel_comp(freq_mat, DMA_element_loc, user_loc, boresight_gain, sub_number);
    
                fval0 = 0;
                w_beam = rand(RFC_num, sub_number).*exp(1i.*rand(RFC_num, sub_number));
                DMA_scenario;
                ym0 = Y_Rcv_Matrix';
                zdc0 = Z_DC;
                xtr0 = X_TR;
                it_main = 1;
                fval00 = 0;
                while true
                    it_in = 0;
                    DC_Threshold = DC_Threshold_main;
                    DC_Threshold_old = DC_Threshold;
                    fval0 = 100;
                    while true 
                        [ym, zdc, xtr, qvec, fval] = ...
                            fixed_W_SCA(w_beam, ym0, channel_vec, ...
                            As, H_dma, freq_mat, samp_vec);
                        
                        if norm(fval - fval0) < eps_thr
                            break
                        end
                        fval0 = fval;
                        Q = zeros(passive_num*RFC_num, RFC_num);
                        for i = 1:RFC_num
                            for l = 1:passive_num
                                idx = (i-1)*passive_num + l;
                                temp_val = qvec(i, l);
                                [a, b] = min(lorentz_vec - temp_val);
                                Q((i - 1)*passive_num + l, i) = lorentz_vec(b(1));
                                Q_dma(i, l) = lorentz_vec(b(1));
                            end
                        end
                        
                        DMA_scenario;
                        ym0 = Y_Rcv_Matrix';
                        zdc0 = Z_DC;
                        xtr0 = X_TR;
                        it_in = it_in + 1;
                    end
    
                    while true
                        fval0 = 0;
                        while true
                            [w_beam, ym, zdc, xtr, fval] = phase_init_func(ym0, zdc0, channel_vec, As,...
                                Q_dma, H_dma, freq_mat, samp_vec, DC_Threshold, xtr0);
                            
                            if norm(fval - fval0) < eps_thr
                                break
                            end
                            fval0 = fval;                    
                            ym0 = Y_Rcv_Matrix';
                            zdc0 = Z_DC;
                            xtr0 = X_TR;
                        end
                        DMA_scenario;
                        if P_RCV_Users >= DC_Threshold_old
                            break
                        else
                            DC_Threshold = (DC_Threshold./P_RCV_Users).*DC_Threshold;
                        end
                    end
                    DMA_scenario;
                    ym0 = Y_Rcv_Matrix';
                    zdc0 = Z_DC;
                    xtr0 = X_TR;
                    it_main = it_main + 1;
                    if norm(fval - fval00) < 10*eps_thr || it_main == it_main_max
                        break
                    end
                    fval00 = fval;
                end
                DMA_scenario;
                P_Cons_Results_DMA_overN(ll, mm, nn, user_it) = P_Consumption;
                Pin_Cons_Results_DMA_overN(ll, mm, nn, user_it) = P_in;
                PA_Cons_Results_DMA_overN(ll, mm, nn, user_it) = PA_consumption;

                disp(strcat('User_', string(user_it), '_Length_', string(ll), '_User_num_', string(user_num), ...
                    '_sub_carrier_', string(sub_number)))
                filename = strcat('DMA_OverN_User_', string(user_it), '_Length_', string(ll), '_User_num_', string(user_num), ...
                    '_sub_carrier_', string(sub_number), '.mat');
                save(filename);
            end
        end
    end
end

