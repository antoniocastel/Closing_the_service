


%calculates Difference between the hazard rate @ t_delta and (mu-delta)/q

function diff = Pol_7_sol(delta, k, thta, mu, q, beta)

    t_delta = log(mu/(mu - delta))/beta;

    numerator = exp(-t_delta/thta) * (t_delta/thta)^(k - 1);
    denominator = thta * gamma(k) * gammainc(t_delta/thta, k, 'upper');
    h_t_delta = numerator / denominator;
    target = max(mu - (h_t_delta * q),0) ;
    
    diff = (delta- target)^2;


end



