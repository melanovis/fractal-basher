function dconv_total = run_multifractal(fractal_population, parameters, c_plane, fractal_iters)

blowup_cutoff = 1e30;

dconv_total = single( zeros(size(c_plane)) );

frac_index_cumulative = 0;

for ind_frac = 1:fractal_population

    % a = 2 + 1.5*j;
    % c = 1.609 + 1.2*j;
    % t1 = 1 + 0*j;
    % t2 = 0 + 0*j;
    a = parameters(ind_frac+frac_index_cumulative,1) + parameters(ind_frac+frac_index_cumulative,2)*j;
    c = parameters(ind_frac+frac_index_cumulative+1,1) + parameters(ind_frac+frac_index_cumulative+1,2)*j;
    t1 = parameters(ind_frac+frac_index_cumulative+2,1) + parameters(ind_frac+frac_index_cumulative+2,2)*j;
    t2 = parameters(ind_frac+frac_index_cumulative+3,1) + parameters(ind_frac+frac_index_cumulative+3,2)*j;
    
    frac_index_cumulative = frac_index_cumulative+3;

    start_plane = c_plane.*exp(t1) + t2; %start locations
    z_n = start_plane;
    conv_iters = single( zeros(size(c_plane)) );
    d_p_stop = conv_iters;
    for n=1:fractal_iters
        z_np = z_n.^a + c;
        d_p = abs(z_np-z_n);
        d_p_mask = d_p < blowup_cutoff;
    
        if all( isfinite( reshape(d_p,1,[]) ) )
            d_p_stop = d_p;
        end
    
        conv_iters(~d_p_mask) = conv_iters(~d_p_mask) + 1;
        z_n = z_np;
    end
    d_p_stop = log10(d_p_stop);
    d_conv = -d_p_stop + max(max(abs(d_p_stop)));
    d_conv(conv_iters==0) = 0;

    dconv_total = dconv_total + d_conv;
end



end