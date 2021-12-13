function v_batch_returned = estimate_V(i_batch, vocv_batch, last_i, ...
                last_v, last_vocv, dt, Rs, Rp, C, Rs_scale, Rp_scale, C_scale)

    %   estimate_V estimates V given parameters and Current 
    %   and Voltage measurements
    %
    %   INPUT:
    %       i_batch: I(t)
    %       vocv_batch: Vocv(t)
    %       last_i: Current measurement before the data batch under
    %       analisys
    %       last_v: Voltage measurement before the data batch under
    %       analisys
    %       last_vocv: last value of Vocv estimated
    %       dt: sampling interval
    %       Rs: Rs value
    %       Rp: Rp value
    %       C: C value
    %       Rs_scale: Rs scale value
    %       Rp_scale: Rp scale value
    %       C_Scale: C scale value
    %
    %   OUTPUT:
    %       vocv_batch_returned: Vocv(t)

    assert(length(vocv_batch) == length(i_batch))
    C = C*C_scale;
    Rp = Rp*Rp_scale;
    Rs = Rs*Rs_scale;
    v_batch(1) = last_v;
    vocv_batch = [last_vocv; vocv_batch];
    i_batch = [last_i; i_batch];
    
    
    for ii = 2:length(i_batch)
        v_batch(ii) = ( ...
            ( 1/dt ) * v_batch(ii-1) + ...
            ( (1/dt) + (1/(C*Rp)) ) * vocv_batch(ii) - ...
            ( 1/dt ) * vocv_batch(ii-1) - ...
            ( (Rs/dt) + (1/C) + (Rs/(Rp*C)) ) * i_batch(ii) + ...
            ( Rs/dt ) * i_batch(ii-1)   ) ...
            / ((1./dt)+(1/(C.*Rp)));
    end

    v_batch_returned = v_batch(2:length(i_batch));
    
end