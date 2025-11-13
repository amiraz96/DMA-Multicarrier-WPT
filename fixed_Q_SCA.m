function [w_beam, zdc, yobj, f_opt, cvx_status] = fixed_Q_SCA(w0, yobj0, zdc0, H_vec_w, Q_dma, H_dma, DC_Threshold, deltaval)

i_s = 5*10^(-6); eta_0 = 1.05; V_o = 25.86 * 10^(-3); R_ant = 50; R_L = 10e3;  
K2 = R_ant*Rect_K(2, i_s, V_o, eta_0); K4 = (R_ant^2)*Rect_K(4, i_s, V_o, eta_0);

G = 1;
user_num = size(H_vec_w, 1);
sub_number = size(H_vec_w, 2);
N = size(H_vec_w, 3);
RFC_num = size(Q_dma, 1);
passive_num = size(Q_dma, 2);
rho = 1;
yobj = [];
cvx_begin 
cvx_solver('mosek')
    variable zdc(user_num, 1) nonnegative
    variable w_beam(RFC_num, sub_number) complex
%     variable yobj(RFC_num, 1) nonnegative
    obj = 0;
    obj2 = 0;
    for nn = 1:sub_number
        obj2 = obj2 + w_beam(:, nn)'*w_beam(:, nn);
    end
    for kk = 1:RFC_num
        obj_und = 0;
        obj_und0 = 0;
        for lll = 1:passive_num
            w_beam_new(kk, lll, :) = w_beam(kk, :);
            w0_new(kk, lll, :) = w0(kk, :);
            for nn = 1:sub_number
                obj_und = obj_und + 0.5*(Q_dma(kk, lll)*H_dma(kk, lll)*w_beam(kk, nn))*(Q_dma(kk, lll)*H_dma(kk, lll)*w_beam(kk, nn))';
                obj_und0 = obj_und0 + 0.5*(Q_dma(kk, lll)*H_dma(kk, lll)*w0(kk, nn))*(Q_dma(kk, lll)*H_dma(kk, lll)*w0(kk, nn))';
            end
        end
        obj = obj + sqrt(obj_und0) + (1/(2*sqrt(obj_und0)))*(obj_und - obj_und0);
%         yobj(kk) >= obj_und;
%         obj = obj + sqrt(yobj0(kk)) + (1/(2*sqrt(yobj0(kk))))*(yobj(kk) - yobj0(kk));
    end
    w_beam_new = reshape(w_beam_new, [N sub_number]);
    w0_new = reshape(w0_new, [N sub_number]);
%     minimize(obj+ obj2 + (rho/3)*pow_pos(reshape(w_beam - w0, [RFC_num*sub_number 1])'*reshape(w_beam - w0, [RFC_num*sub_number 1]), 3) + ...
%         (rho/3)*pow_pos((zdc - zdc0)'*(zdc - zdc0), 3))   
    minimize(obj+ obj2)  
    subject to
        for minn = 1:user_num
            h = reshape(H_vec_w(minn, :, :), [sub_number N]);
            [coeff, f0value] = taylor_coef_w(w0_new, h, K2, K4, sub_number, N);
            consval = 0;
            for nin = 1:sub_number
                consval = consval + (coeff(nin, :)*(w_beam_new(:, nin) - w0_new(:, nin)));
            end
            zdc(minn) <= real(f0value + consval);
        end
        for mu = 1:user_num
            R_L*DC_Threshold <= zdc0(mu)^2 + 2*zdc0(mu)*(zdc(mu) - zdc0(mu));
        end
        norm(reshape(w_beam - w0, [RFC_num*sub_number 1])) <= deltaval
cvx_end

f_opt = obj+ obj2;
    
end
