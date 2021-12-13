function vocv_batch_returned = estimate_Vocv( i_batch, v_batch, last_i, ...
                last_v, last_vocv, dt, Rs, Rp, C, Rs_scale, Rp_scale, C_scale)
    
    %   estimate_Vocv estimates Vocv given previously detected
    %   parameters and Current and Voltage measurements
    %
    %   INPUT:
    %       i_batch: I(t)
    %       v_batch: V(t)
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

    assert(length(v_batch) == length(i_batch))
    C = C*C_scale;
    Rp = Rp*Rp_scale;
    Rs = Rs*Rs_scale;
    vocv_batch(1) = last_vocv;
    v_batch = [last_v; v_batch];
    i_batch = [last_i; i_batch];
    
    
    for ii = 2:length(i_batch)
        vocv_batch(ii) = ( ...
            ((1/dt)+(1/(C*Rp))).*v_batch(ii) ...
            - (v_batch(ii-1)./dt) ...
            + (vocv_batch(ii-1)./dt) ...
            + (((Rs./dt)+(1./C)+(Rs./(C.*Rp))).*i_batch(ii)) ...
            - ((Rs./dt).*i_batch(ii-1))    ) ...
            / ((1./dt)+(1/(C.*Rp)));
    end

    vocv_batch_returned = vocv_batch(2:length(i_batch));
    
end

