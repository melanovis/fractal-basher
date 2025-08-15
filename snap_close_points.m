function snapped = snap_close_points(z, x)
    bin_size = x;
    bins = containers.Map();
    snapped_list = [];       
    for n = 1:length(z)
        real_part = real(z(n));
        imag_part = imag(z(n));
        
        ind_x = floor(real_part / bin_size);
        ind_y = floor(imag_part / bin_size);
        key = sprintf('%d_%d', ind_x, ind_y);
        
        if isKey(bins, key)
            ind_k = bins(key);
            snapped_list(ind_k) = (snapped_list(ind_k) + z(n)) / 2;
        else
            snapped_list(end+1) = z(n);
            bins(key) = numel(snapped_list);
        end
    end
    snapped = snapped_list(:);
end
