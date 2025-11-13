

%% Rectifier Parameters
i_s = 5*10^(-6); eta_0 = 1.05; V_o = 25.86 * 10^(-3); R_ant = 50; R_L = 10e3;  
K2 = R_ant*Rect_K(2, i_s, V_o, eta_0); K4 = (R_ant^2)*Rect_K(4, i_s, V_o, eta_0);

%% Parameters

N = RFC_num*passive_num;
sig_weight = w_beam;

QH_dma = H_dma.*Q_dma;

%% HPA Parameters
eff_max = pi/4;

BW = 10e6;
delta_f = BW/sub_number;
delta_f = 125e3;
freq_mat = zeros(1, sub_number); %% Sub-Carriers
for ii = 1:sub_number 
    freq_mat(ii) = f_0 + (ii - 1)*delta_f;
end

nfac = 20000;
samp_total = nfac/max(freq_mat);
samp_int = samp_total/(2*nfac*sub_number);
samp_vec = 0:samp_int:samp_total;
Nt = length(samp_vec);

   
%% Multi-Sine Waveform
Y_Rcv_Matrix = zeros(Nt, user_num);
X_TR = zeros(Nt, RFC_num, passive_num);

P_in = sum((abs(w_beam)).^2, 'all');


for tt = 1:Nt
    for kk = 1:RFC_num
        for lll = 1:passive_num
            for nnn = 1:sub_number
                X_TR(tt, kk, lll) = X_TR(tt, kk, lll) + real(Q_dma(kk, lll)*H_dma(kk, lll)*w_beam(kk, nnn)*exp(1i*2*pi*freq_mat(nnn)*samp_vec(tt)));
            end
        end
    end
end

P_out_HPA = reshape(sum(abs(X_TR).^2, 3), [length(samp_vec), RFC_num]);
HPA_eff_Matrix = eff_max.*sqrt(P_out_HPA./(As^2));
HPA_Cons_Matrix = P_out_HPA./HPA_eff_Matrix;

for tt = 1:Nt
    for mm = 1:user_num
        for kk = 1:RFC_num
            for lll = 1:passive_num
                for nnn = 1:sub_number
                    Y_Rcv_Matrix(tt, mm) = Y_Rcv_Matrix(tt, mm) + real(channel_vec(mm, nnn, kk, lll)*...
                        Q_dma(kk, lll)*H_dma(kk, lll)*w_beam(kk, nnn)*exp(1i*2*pi*freq_mat(nnn)*samp_vec(tt)));
                end
            end
        end
    end
end

Z_DC = K2*(mean(Y_Rcv_Matrix.^2, 1)) + K4*(mean(Y_Rcv_Matrix.^4, 1));
PA_consumption = sum(mean(HPA_Cons_Matrix, 1));
P_BB = 200e-3; P_RFC = 75e-3 + 23e-3 + 5e-3;
P_Consumption = P_in + PA_consumption + P_BB + RFC_num*P_RFC;
P_RCV_Users = (Z_DC.^2)/R_L;

max_cons = eff_max*(As^2);
max_P = P_in + max_cons*RFC_num + P_BB + RFC_num*P_RFC;
HPA_eff = mean(mean(HPA_eff_Matrix, 2));


