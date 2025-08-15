function [converged_root, return_matrix, roots_return, fitness] = NR_iter_PSO(start_location, input_vector, convergence_tolerance, max_iters, target_image)

    input_vector = reshape(input_vector,4,[]);
    input_scaling_mag = 4;
    roots =  interp1([0,0.09,0.1,1],[0,0, -10, 10],input_vector(1,:)) .* 10.^interp1([0,1],[-input_scaling_mag,input_scaling_mag],input_vector(2,:)) + interp1([0,0.09,0.1,1],[0,0,-10,10],input_vector(3,:)) .* j.*10.^interp1([0,1],[-input_scaling_mag,input_scaling_mag],input_vector(4,:));
    
    view_transform = roots(1:3); %first three inputs are transformations 
    start_location = (1 - (view_transform(3)/2)*j).*start_location + (view_transform(3)/2)*j.*conj(start_location); %shear
    start_location = start_location.*exp(view_transform(2)) + view_transform(1); %rotation and translation

    roots = roots(4:end);
    roots_return = [view_transform, roots];
    
    roots = snap_close_points(roots, 0.1);

    f_poly = poly(roots);
    fp_poly = polyder(f_poly);

    conv_iterations = zeros(size(start_location));
    dx = ones(size(start_location));
    convergence_mask = zeros(size(start_location));
    plane_n = start_location;

    iters = 0;
    while iters < max_iters
        plane_np = plane_n - ( polyval(f_poly,plane_n) ) ./ ( polyval(fp_poly,plane_n) );
        plane_np =  sinh(plane_np);
        dx = abs(plane_np - plane_n);
        convergence_mask = dx < convergence_tolerance;
        conv_iterations(convergence_mask) = conv_iterations(convergence_mask)+1;
        plane_n = plane_np;
        iters = iters + 1;
    end

    d_conv = abs(start_location - plane_np).*conv_iterations;
    %d_conv = conv_iterations;
    [~, converged_root] = min(abs(plane_np - reshape(roots, 1, 1, [])), [], 3);

    fractal_masked_out = converged_root;
    fractal_masked_out(target_image==0) = 0;
    largest_area_index = mode(fractal_masked_out(fractal_masked_out~=0),"all");
    fractal_fully_masked = single(converged_root);
    fractal_fully_masked(converged_root~=largest_area_index) = 0;
    fractal_fully_masked(converged_root==largest_area_index) = 1;
    
    fitness = double(1 / immse(target_image, fractal_fully_masked));

    return_matrix(:,:,1) = converged_root;
    return_matrix(:,:,2) = d_conv;
    return_matrix(:,:,3) = fractal_fully_masked;
end
