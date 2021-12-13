function l = loss_function(x, u, dt, lookup, scale_factors)
    
    %   loss_function estimate the value of the loss given state x 
    %   and others parameters
    %
    %   INPUT:
    %       x: state vector to optimize [Rs, Rp, C, SoC(tau), Qmax]
    %       u: measurement vector
    %       dt: sampling interval
    %       lookup: lookup table
    %       scale_factors: scale factors vector for Rs, Rp and C
    %
    %   OUTPUT:
    %       l: loss value

    SoC2Vocv = @(x)(get_Vocv(x, lookup));
    
    soc = x(4) + (u(:,3).*(dt/x(5)));
    vocv = SoC2Vocv(soc);
    
    V_est = estimate_V(u(2:end,2), vocv(2:end), u(1,2), u(1,1), ...
            vocv(1), dt, x(1), x(2), x(3), ...
            scale_factors(1), scale_factors(2), scale_factors(3));
            
    l = sum(abs(V_est'-u(2:end,1)));
    
end