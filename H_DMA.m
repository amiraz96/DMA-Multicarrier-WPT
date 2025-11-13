function h = H_DMA(f_0, antenna_length, element_loc)

eps_r = 3.54;   %% Datasheet
tans = 0.0021; %% Datasheet
NP_dB = 8.685889638; %% Neper to dB
d_ms = 0.5258e-3; %% Datasheet Tickness
cop = 5.96e7; %% Copper Conductivity
W = antenna_length;
c0 = 3e8; %% Light Speed
lambda = c0/f_0;
wd_ms = W/d_ms;
k0 = (2*pi)/(lambda);
eps_e_0 = (eps_r/2 + 0.5) + (eps_r/2 - 0.5)*(1/(sqrt(1 + 12*wd_ms)));
Z_0 = (120*pi)/(sqrt(eps_e_0) * (wd_ms + 1.393 + 0.667*log(wd_ms + 1.444)));
f_p = Z_0/(8*pi*d_ms*100);
g = 0.6 + 0.009*Z_0;
f_0_GHz = f_0/(10^(9));
Gf = g*((f_0_GHz/f_p)^2);
eps_e_f = eps_r - (eps_r - eps_e_0)/(1 + Gf);
alphaa_d = (k0*eps_r*(eps_e_f - 1)*tans)/(2*(eps_r - 1)*sqrt(eps_e_f));
R_s = sqrt((2*pi*f_0*4*pi*(10^-7))/(2*cop));
alphaa_c = R_s/(Z_0*W);
alphaa = 10^(NP_dB*(alphaa_c + alphaa_d)/(10));
betaa = k0*sqrt(eps_e_f);

h = exp(-element_loc*(alphaa + 1i*betaa));

end