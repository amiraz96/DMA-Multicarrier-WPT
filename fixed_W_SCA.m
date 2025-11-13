function [qvec, fval] = fixed_W_SCA(q0, H_vec_w, deltaval)

i_s = 5*10^(-6); eta_0 = 1.05; V_o = 25.86 * 10^(-3); R_ant = 50; R_L = 10e3;  
K2 = R_ant*Rect_K(2, i_s, V_o, eta_0); K4 = (R_ant^2)*Rect_K(4, i_s, V_o, eta_0);

rho = 1;
G = 1;
user_num = size(H_vec_w, 1);
sub_number = size(H_vec_w, 2);
N = size(H_vec_w, 3);

cvx_begin 
cvx_solver('mosek')
    variable zdc(user_num, 1) nonnegative
    variable qvec(N, 1) complex
%     maximize(min(zdc) - (rho/3)*(pow_pos(sum_square_abs(qvec - q0), 3)))
    maximize(min(zdc))
    subject to
        for minn = 1:user_num
            h = reshape(H_vec_w(minn, :, :), [sub_number N]);
            [coeff, f0value] = taylor_coef_q(q0, h, K2, K4, sub_number, N);
            consval = coeff*(qvec - q0);
            zdc(minn) <= real(f0value + consval);
        end
        for k = 1:N
            (real(qvec(k)))^2 + (imag(qvec(k)) - 0.5)^2 <= 0.5^2;
        end
        norm(qvec - q0) <= deltaval

cvx_end

fval  = min(zdc);


end
