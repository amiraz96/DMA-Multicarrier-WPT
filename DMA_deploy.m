function [DMA_element_loc, H,...
    RFC_num, passive_num] = DMA_deploy(f, antenna_length_vec, x_max, y_max, z_max)

%% Initial Parameters
c0 = physconst('LightSpeed');
lambda = c0/f;
RFC_num = floor(antenna_length_vec(1) / (0.5*lambda));  %% Number of RF Chains
passive_num = floor(antenna_length_vec(2) / (0.2*lambda));
x_min = 0; y_min = 0;


%% Deploy the DMA and Calculate element Loss !

DMA_WG_dist = lambda/2; DMA_element_dist = lambda/5;
H = zeros(passive_num*RFC_num, passive_num*RFC_num);
for i = 1:RFC_num
    for l = 1:passive_num
        H((i-1)*passive_num + l, (i-1)*passive_num + l) = H_DMA(f(1), ...
            antenna_length_vec(2), (l-1)*DMA_element_dist);
    end
end

DMA_loc = [(x_min + x_max)/2 - floor(passive_num/2)*DMA_element_dist ...
    (y_min + y_max)/2 - floor(RFC_num/2)*DMA_WG_dist z_max];

DMA_element_loc = zeros(RFC_num*passive_num, 3);

for i = 1:RFC_num
    for j = 1:passive_num
        DMA_element_loc((i-1)*passive_num + j, :) = [DMA_loc(1) + (j - 1)*DMA_element_dist ...
            DMA_loc(2) + (i - 1)*DMA_WG_dist, DMA_loc(3)];
    end
end

DMA_element_loc = reshape(DMA_element_loc, [RFC_num, passive_num, 3]);

end