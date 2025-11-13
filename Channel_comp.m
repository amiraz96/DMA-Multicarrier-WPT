function [channel_vec] = Channel_comp(freq_mat, DMA_element_loc, user_loc, boresight_gain, sub_number)

c0 = physconst('LightSpeed');
RFC_num = size(DMA_element_loc, 1);
passive_num = size(DMA_element_loc, 2);
user_num = size(user_loc, 1);


%% Channel Coefficients

channel_vec = zeros(user_num, sub_number, RFC_num, passive_num);
for m = 1:user_num
    for n = 1:sub_number
        for i = 1:RFC_num
            for l = 1:passive_num
                d3D  = sqrt(sum((reshape(DMA_element_loc(i, l, :), [1, 3]) - user_loc(m, :)) .^ 2));
                dz = norm(DMA_element_loc(i, l, 3) - user_loc(m, 3));
                channel_rad = 2*(boresight_gain+1)*...
                    ((dz/d3D))^boresight_gain;
                lambda_n = c0/freq_mat(n);
                A_channel = sqrt(channel_rad) * (lambda_n/(4*pi*d3D));
                phase_channel = d3D/lambda_n;
                channel_vec(m, n, i, l) = A_channel*exp(-1i*2*pi*phase_channel);
            end
        end
    end
end

end