function [Rs, Rp, C, SoC_tau, Qmax] = she_estimate_battery_model(i_batch, v_batch, ...
                    dt, algorithm, number_restarts, lookup, scale_factors, verbose)
    
    %   estimate_all_params estimate values of equvalent model parameters
    %
    %   INPUT:
    %       i_batch: array of Current measurements
    %       v_batch: array of Voltage measurements
    %       dt: sampling interval
    %       algorithm: algorithm used by fmincon from the set of the matlab
    %       implemented ones
    %       number_restarts: number of restarts to avoid local minima
    %       lookup: lookup table
    %       scale_factors: array [Rs_scale Rp_scale C_scale] scale factors
    %       verbose: verbose integer (1->True, 0->False)
    %
    %   OUTPUT:
    %       Rs, Rp, C, SoC(tau), Qmax
    
    len_i = length(i_batch);
    len_v = length(v_batch);
    
    assert( len_i == len_v );
    
    init_time = 1;
    time = len_i;
    in(1:time-init_time+1, 1) = v_batch(init_time:time);
    in(1:time-init_time+1, 2) = i_batch(init_time:time);
    in(:,3) = - cumsum(i_batch);
    
    % closure to intergrate mandatory information
    to_minimize = @(x)loss_function(x, in, dt, lookup, scale_factors);
    
    % options of the optimizer fmincon
    if verbose == 1
        display_details = 'iter-detailed';
    else
        display_details = 'none';
    end
    minOptions = optimoptions('fmincon', 'Display', display_details, ...
                'Algorithm', algorithm, 'OptimalityTolerance', 1e-10, ...
                'StepTolerance', 1e-10, 'MaxFunctionEvaluations', 1e3, ...
                'MaxIterations', 1e3);
    
    % bounds of the fmincon, the real values are rescaled in order to
    % fit all the possible values of the Rs, Rp, and C eq model
    lb = [1e-2, 1e-2, 1e-2, 1e-2, 1e4];
    ub = [1e2, 1e2, 1e2, 1, 1e6];
    
    % random init in the range of the bounds 
    % (the values are rescaled in the next phase)
    randmatrix = rand(number_restarts, length(ub));
    for jj = 1:length(ub)
       randmatrix(:,jj) = randmatrix(:,jj) .* (ub(jj) - lb(jj)) + lb(jj);  
    end
    
    % performs the multiple restart to avoid local minima
    for jj = 1:number_restarts
        try
        [x(jj,:), fval(jj)] = fmincon(to_minimize, randmatrix(jj,:)', ...
                            [], [], [], [], lb, ub, [], minOptions);
        catch
            fval(jj) = 1e10;
        end
    end
    
    % select the best estimation among the restarts and correct the values
    % according to scale factors, if required
    [~, idx] = min(fval);
    Rs = x(idx, 1) * scale_factors(1);
    Rp = x(idx, 2) * scale_factors(2);
    C = x(idx, 3) * scale_factors(3);
    SoC_tau = x(idx, 4);
    Qmax = x(idx, 5);  
    
end
