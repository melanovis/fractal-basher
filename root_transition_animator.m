format compact
clear
clc
clf reset
%close all


record_video = true;

pre_render_view = false;
render_mask = false;

animation_frames = 100;

filename = "root_catalog.mat";
load(filename)

if record_video == true
    if ~render_mask
        v = VideoWriter("root_transform_animated", 'MPEG-4');
    else
        v = VideoWriter("root_transform_transparency_mask", 'MPEG-4');
    end
    v.FrameRate = 30;
    open(v);
end

convergence_tolerance = 1e-7;
max_iters = 2e2;

dims = ceil([1920,1080]./20); %control res
aspect_ratio = dims(1)/dims(2);
view_domain_x = [-1,1];
view_domain_y = view_domain_x./aspect_ratio;
canvas_x = linspace(view_domain_x(1),view_domain_x(2),dims(1));
canvas_y = linspace(view_domain_y(1),view_domain_y(2),dims(2));

domain_x = [-1,1];
domain_y = domain_x./aspect_ratio;

[x_plane, y_plane] = meshgrid(canvas_x, canvas_y);
complex_plane = x_plane + y_plane.*j; 

root_quantity = 50;

bounds = 1;
[root_grid_default_x,root_grid_default_y] = meshgrid(linspace(-bounds,bounds,floor(sqrt(root_quantity))), linspace(-bounds,bounds,floor(sqrt(root_quantity))));
root_grid_default = root_grid_default_x + root_grid_default_y.*j;
root_grid_default = reshape(root_grid_default,1,[]);

while length(root_grid_default) < root_quantity
    root_grid_default = [root_grid_default,0];
end

root_sequence = [
0, -3, 0, root_grid_default %first three are transformations
-j*pi/2, pi/2, 0, root_grid_default
root_map(find(root_names=="door_1"),:)
-j*pi/2, pi/2, 0, root_grid_default
];

% convergence map indexes which will be painted black when moving towards that root configuration, 0 corresponds to no indexes
mask_outs = [ 
0
0
0
0
];

if height(mask_outs) ~= height(root_sequence)
    height(mask_outs)
    height(root_sequence)
    error("mask outs and root sequence have inequal rows")
end

mask_out_indexes = round(linspace(1,height(mask_outs),animation_frames));

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

for ind_frame = 1:animation_frames

    for n=1:width(root_sequence)

        %may need to find a better way to find easing in future, should be pre-computed
        sequence = root_sequence(:,n);
        time_range = 0:length(sequence)-1;
        time_interp = linspace(min(time_range),max(time_range),animation_frames);
        sequence_nearest = interp1(time_range,sequence,time_interp,"next");
        tf_offset_real = sequence_nearest(1);

        swing_delay = 1;
        sequence_nearest = [repelem(sequence(1),swing_delay),sequence_nearest];
        sequence_nearest(end-(swing_delay-1):end) = [];

        for m=1:length(time_range)
            [~,ind_nearest(m)] = min(abs(time_interp-time_range(m)));
        end
        
        k = 6e2;
        numerator = [1, k];
        denominator = [2, 100, k];
        sys = tf(numerator,denominator);
        response = lsim(sys, sequence_nearest-tf_offset_real, time_interp); 
        response = response+tf_offset_real;
        for m=1:length(time_range)
            response(ind_nearest(m)) = sequence(m);
        end

        roots(n) = response(ind_frame);
    end

    [converged_root, d_conv] = NR_iter_frame(complex_plane, roots, convergence_tolerance, max_iters);
    converged_root = converged_root-1;

    upcoming_maskout = mask_outs(mask_out_indexes(ind_frame));

    rootmap_with_maskout = zeros(size(converged_root)).*nan;
    rootmap_full_maskout = zeros(size(converged_root));

    if upcoming_maskout ~= 0
        for m=1:height(converged_root)
            rootmap_with_maskout = [converged_root(m,:) == upcoming_maskout].';
            d_conv(m,rootmap_with_maskout) = nan;

            rootmap_full_maskout(m,:) = rootmap_with_maskout;
        end
    end

    colormap(flip(cmap))
    if ~render_mask
        if pre_render_view
            subplot(2,1,1)
        end
        scatter(nan,nan,"w+")
        hold on
        grid on
        if pre_render_view
            axis tight equal
        else
            axis vis3d equal
            axes('Units', 'normalized', 'Position', [0 0 1 1]) 
            d_conv = flip(d_conv);
            rootmap_full_maskout = flip(rootmap_full_maskout);
            imagesc(canvas_x, canvas_y, rootmap_full_maskout, alphadata = rootmap_full_maskout ~= 0)
        end

        imagesc(canvas_x, canvas_y, d_conv.^0.75, alphadata = ~isnan(d_conv))

        set(gca,'Color','k')
        clim([1,max_iters]);
        set(gca,'TickLength',[0 0])
        xlim(domain_x)
        ylim(domain_y)
        if ismember(ind_frame, ind_nearest) && ind_frame ~= 1
            legend(" ", Interpreter="latex", FontSize=25, location= "northwest")
            legend boxoff
        end
        hold off
    
        if pre_render_view
            subplot(2,1,2)
            scatter(nan,nan)
            hold on
            grid on
            axis tight equal
            set(gca,'Color','k')
            surf(canvas_x, canvas_y, converged_root, EdgeColor="none")
            set(gca,"ColorScale","log")
            xlim(domain_x)
            ylim(domain_y)
            hold off
        end
    
        if ismember(ind_frame, ind_nearest) && ind_frame ~= 1 && pre_render_view
            sound(sin(2*pi*400*(0:1/14400:0.1)), 14400);
            disp("paused")
            pause
        end
    else
        colormap("gray")
        rootmap_full_maskout = flip(rootmap_full_maskout);
        scatter(nan,nan,"w+")
        hold on
        grid on
        axis vis3d equal
        axes('Units', 'normalized', 'Position', [0 0 1 1]) 
        imagesc(canvas_x, canvas_y, rootmap_full_maskout, alphadata = rootmap_full_maskout ~= 0)
        set(gca,'Color','k')
        xlim(domain_x)
        ylim(domain_y)
        set(gca,'TickLength',[0 0])
        hold off
    end

    drawnow()

    fprintf("rendered frame %i.\n",ind_frame)

    if record_video == true
        frame = getframe(gcf);
        frame.cdata = imresize(frame.cdata, [1895,3840]);
        writeVideo(v,frame)
    end
    
end

if record_video == true
    close(v);
end
sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);



function [converged_root, d_conv] = NR_iter(start_location, roots, convergence_tolerance, max_iters)

    view_transform = roots(1:3);

    start_location = start_location.*exp(view_transform(2)) + view_transform(1); %rotation and translation

    if abs(view_transform(3)) > pi/2
        view_transform(3) =  (pi/4).*(view_transform(3)./norm(view_transform(3)));
    end

    %start_location = (1 - (view_transform(3)/2)*j).*start_location + (view_transform(3)/2)*j.*conj(start_location); %shear is currently disabled
    
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



