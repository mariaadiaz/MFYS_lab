function next_value = gradient_d(f_Fn, Fn)
% When does it converge? 
% 3 points to compute the next derivative, check

    alpha = 1000; % Specifies the fixed step size
    %----- Gradient descent
    Jn = (f_Fn(end) - f_Fn(end-1))/(Fn(end) - Fn(end-1)); 
    %----- Calculate the next step value using the defined alpha
    next_value = Fn(end) - (alpha * Jn); 
    % If the algorithm does not converge then an extra variable would have to
    % be tune as well
end
