function [Vocv] = get_Vocv(SoC,lookup_table)

    Vocv = interp1(lookup_table.SoC,lookup_table.Vocv,SoC, 'linear', 'extrap');
    
end