function [final_coef, f0] = taylor_coef_q(q0, h, K2, K4, sub_number, N)

%h should be a row vector
coeff = zeros(1, N);
first_term = 0;
for nnn = 1:sub_number
    first_term = first_term + (h(nnn, :)*q0)*h(nnn, :);
    for nnn1 = 1:sub_number
        for nnn2 = 1:sub_number
            for nnn3 = 1:sub_number
                if nnn + nnn1 == nnn2 + nnn3
                    coeff = coeff +  (h(nnn1, :)*q0)*...
                            (h(nnn2, :)*q0)'*(h(nnn3, :)*q0)'.*(h(nnn, :)) + (h(nnn, :)*q0)*...
                            (h(nnn2, :)*q0)'*(h(nnn3, :)*q0)'.*(h(nnn1, :)) + (h(nnn, :)*q0)*...
                            (h(nnn1, :)*q0)*(h(nnn3, :)*q0)'.*conj(h(nnn2, :)) + (h(nnn, :)*q0)*...
                            (h(nnn1, :)*q0)*(h(nnn2, :)*q0)'.*conj(h(nnn3, :));
                end
            end
        end
    end
end
final_coef = K2*first_term + (3/8)*K4.*coeff;

expw = 0;
expw2 = 0;

for nnn = 1:sub_number
    expw = expw + 0.5*K2*norm(h(nnn, :)*q0)^2;
    for nnn1 = 1:sub_number
        for nnn2 = 1:sub_number
            for nnn3 = 1:sub_number
                if nnn + nnn1 == nnn2 + nnn3
                    expw2 = expw2 +  ((3/8)*K4*(h(nnn, :)*q0)*(h(nnn1, :)*q0)*...
                            (h(nnn2, :)*q0)'*(h(nnn3, :)*q0)');
                end
            end
        end
    end
end

f0 = expw + expw2;
rfpow = expw/K2;

end




