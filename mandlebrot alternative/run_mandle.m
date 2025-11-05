function d_conv = run_mandle(res, dims_x, a, c, t1, t2)

iters = 20;
blowup_cutoff = 1e12;

aspect_ratio = res(2)/res(1);
dims_y = dims_x.*aspect_ratio;

x_range = single( linspace(dims_x(1), dims_x(2), res(1)) );
y_range = single( linspace(dims_y(1), dims_y(2), res(2)) );

[x_plane, y_plane] = meshgrid(y_range, x_range);
c_plane = x_plane + y_plane.*j; 


start_plane = c_plane.*exp(t1) + t2; %start locations

z_n = start_plane;
conv_iters = single( zeros(size(c_plane)) );
d_p_stop = conv_iters;
for n=1:iters
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

end