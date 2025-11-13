function [w0, Q_dma] = phase_init_func(H_dma, channel_vec)

i_s = 5*10^(-6); eta_0 = 1.05; V_o = 25.86 * 10^(-3); R_ant = 50; R_L = 10e3;  
K2 = R_ant*Rect_K(2, i_s, V_o, eta_0); K4 = (R_ant^2)*Rect_K(4, i_s, V_o, eta_0);
DC_Threshold_main = 2e-5;

RFC_num = size(H_dma, 1);
passive_num = size(H_dma, 2);
sub_number = size(channel_vec, 2);
user_num = size(channel_vec, 1);
N = RFC_num*passive_num;
% Q_dma = zeros(RFC_num, passive_num);
% phi_vec = 0.25.*ones(RFC_num, passive_num);
% for i = 1:RFC_num
%     for l = 1:passive_num
%         Q_dma(i, l) = Q_DMA(2*pi*phi_vec(i,l));
%     end
% end

%% Allocate RF chains
temp_norm = zeros(user_num, sub_number);
for ii = 1:user_num
    for nn = 1:sub_number
        temp_norm(ii, nn) = norm(reshape(channel_vec(ii,nn, :, :), [N 1]));
    end
end
[highest_chan, b1] = max(temp_norm, [], 2);
[~, b2] = min(highest_chan, [], 1);
srch_vec = linspace(0, 2*pi, 1000);
RF_vec = [];
avail_RFs = RFC_num;
RF_cnt = 0;
for m = 1:user_num
    RF_vec(m, 1) = RF_cnt + 1;
    RF_cnt = RF_cnt + 1;
    avail_RFs = avail_RFs - 1;
end
rf_num = [];
for m = 1:user_num
    chan_ratio = 1 - highest_chan(m)/sum(highest_chan);
    rf_num(m, 2) = floor(chan_ratio*avail_RFs);
    rf_num(m, 1) = m;
end
aloc_rf = [];
rf_num = sortrows(rf_num, 2, 'descend');
for m = 1:user_num
    temp_rf = rf_num(m, 2);
    aloc_rf(m, 1) = rf_num(m, 1);
    aloc_rf(m, 2) = ceil((temp_rf/sum(rf_num(:, 2)))*avail_RFs);
end
if user_num ~= 1
    for m = 1:user_num
        user_idx =  aloc_rf(m, 1);
        for mm = 1:aloc_rf(m, 2)
            RF_vec(user_idx, mm + 1) = RF_cnt + 1;
            RF_cnt = RF_cnt + 1;
            avail_RFs = avail_RFs - 1;
            if avail_RFs == 0
                break
            end
        end
        if avail_RFs == 0
            break
        end
    end
else
    RF_vec = [];
    RF_vec(1, :) = 1:RFC_num;
end

%% Start meeting requirement

w0_phase = zeros(RFC_num, sub_number);
w0_amp = zeros(RFC_num, sub_number);
w_beam = zeros(RFC_num, sub_number);

Q_dma = zeros(RFC_num, passive_num);
for m = 1:user_num
    high_idx = b1(m);
    for i = 1:RFC_num
        if ismember(i, RF_vec(m, :))
            for l = 1:passive_num
                temp_phi = angle(channel_vec(m, high_idx, i,l));
                if -pi <= temp_phi <= 0
                    Q_dma(i, l) = Q_DMA(-1i*log(2*exp(-1i*temp_phi) - 1i));
                else
                    Q_dma(i, l) = Q_DMA(pi/2);
                end
            end
        end
    end
end

for mm = 1:user_num
    myamp = 0.001;
    while true
        f0 = 0;
        for iii = 1:RFC_num
            if ismember(iii, RF_vec(mm, :))
%                 for nf = 1:sub_number
%                     tempval = 0;
%                     for nff = 1:sub_number
%                         tempval = tempval + norm(reshape(channel_vec(mm, nff, :, :), [1 N]))^2;
%                     end
%                     w_beam(iii, nf) = exp(-1i*angle(H_dma(iii, 1)*Q_dma(iii, 1)*...
%                         channel_vec(mm, nf, iii, 1)))*...
%                         norm(reshape(channel_vec(mm, nf, :, :), [1 N]))*...
%                         sqrt((2*myamp)/(N*tempval));
%                 end
                for nn = 1:sub_number
                    f = 0;
                    for i = 1:length(srch_vec)
                        ftemp = 0;
                        for lll= 1:passive_num
                            ftemp = ftemp + ...
                                abs(angle(channel_vec(mm, nn, iii, lll)...
                                *Q_dma(iii, lll)*H_dma(iii, lll)*exp(1i*srch_vec(i))));
                        end
                        if ((ftemp)) > f
                            f = ((ftemp));
                            phistar = srch_vec(i);
                        end
                    end
                    w0_phase(iii, nn) = phistar;
                    w0_amp(iii, nn) = myamp;
                end
            end
        end
        w0 = w0_amp.*exp(1i.*w0_phase);
%         w0 = w_beam;
        w_beam_new = [];
        for kk = 1:RFC_num
            for lll = 1:passive_num
                w_beam_new(kk, lll, :) = w0(kk, :);
            end
        end
        w_beam_new = reshape(w_beam_new, [N, sub_number]);
        H_vec_q = zeros(user_num, sub_number, N);
        for mmm = 1:user_num
            for nnn = 1:sub_number
                atemp = reshape(channel_vec(mmm, nnn, :, :), [N, 1]);
                H_vec_q(mmm, nnn, :) = atemp.*reshape(H_dma, [N, 1]).*w_beam_new(:, nnn);
            end
        end
        h = reshape(H_vec_q(mm, :, :), [sub_number N]);
        [coeff, f0value] = taylor_coef_q(reshape(Q_dma, [N 1]), h, K2, K4, sub_number, N);
        zdc0 = real(f0value);

        P_spec = (zdc0^2)/R_L
        if P_spec > DC_Threshold_main
            break
        else
            myamp = myamp*1.1;
        end
    end
end
disp(RF_vec)
end




