function [converged_root, d_conv] = NR_iter_frame(start_location, roots, convergence_tolerance, max_iters)

    view_transform = roots(1:3); %first three inputs are transformations 
    
    if abs(imag(view_transform(2))) > 2*pi %don't want too much rotation
        view_transform(2) = real(view_transform(2)) + sign(imag(view_transform(2)))*j*2*pi;
    end

    %managing shear
    start_location_tmp = reshape(start_location,1,[]);
    start_location_decomp(1,:) = real(start_location_tmp);
    start_location_decomp(2,:) = imag(start_location_tmp);

    start_location_decomp = start_location_decomp.' * [1, real(view_transform(3)); 0, 1];
    start_location_decomp = start_location_decomp * [1, 0; imag(view_transform(3)), 1];
    
    start_location_tmp = start_location_decomp(:,1) + start_location_decomp(:,2).*j;
    start_location_tmp = reshape(start_location_tmp, size(start_location));
    start_location = start_location_tmp;
    start_location = start_location.*exp(view_transform(2)) + view_transform(1);

    roots = roots(4:end);
    roots = snap_close_points(roots, 0.1);
    
    f_poly = poly(roots);
    fp_poly = polyder(f_poly);

    conv_iterations = zeros(size(start_location));
    d_x = ones(size(start_location));
    convergence_mask = zeros(size(start_location));
    plane_n = start_location;

    iters = 0;
    while iters < max_iters
        plane_np = plane_n - ( polyval(f_poly,plane_n) ) ./ ( polyval(fp_poly,plane_n) );
        plane_np =  sinh(plane_np);
        d_x = abs(plane_np - plane_n);
        convergence_mask = d_x < convergence_tolerance;
        conv_iterations(convergence_mask) = conv_iterations(convergence_mask)+1;
        plane_n = plane_np;
        iters = iters + 1;
    end

    d_conv = abs(start_location - plane_np).*conv_iterations;
    %d_conv = conv_iterations;
    [~, converged_root] = min(abs(plane_np - reshape(roots, 1, 1, [])), [], 3);
end
