function [SoC_array] = cc_SoC(i_batch,init_SoC,dt,q)
    
    % CC_VOCV estimate the SoC using Coulomb Counting method
    %   INPUT:
    %   i_batch: array of current measures
    %   init_SoC: the initial value of the SoC
    %   dt: sampling time
    %   q: capacity of the battery at moment t   
    
    aux = init_SoC;
    
    for ii = 1:length(i_batch)
        
        aux = aux - (1/q)*i_batch(ii)*dt;
        SoC_array(ii) = aux;
    
    end
    
end

