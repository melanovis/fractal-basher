function [converged_root, d_conv] = NR_iter_frame(start_location, roots, convergence_tolerance, max_iters)

    view_transform = roots(1:3);

    start_location = (1 - (view_transform(3)/2)*j).*start_location + (view_transform(3)/2)*j.*conj(start_location); %shear 
    start_location = start_location.*exp(view_transform(2)) + view_transform(1); %rotation and translation
    
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
