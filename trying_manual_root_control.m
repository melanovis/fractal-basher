format compact
clear
clc
%clf reset
close all

%pick up and move roots around to get a feel for how sinh-modified newtons fractal changes.

convergence_tolerance = 1e-9;
max_iters = 20;

dims = ceil([1920,1080]./20);
aspect_ratio = dims(1)/dims(2);
view_domain_x = [-5,5];
view_domain_y = view_domain_x./aspect_ratio;
canvas_x = linspace(view_domain_x(1),view_domain_x(2),dims(1));
canvas_y = linspace(view_domain_y(1),view_domain_y(2),dims(2));

domain_x = [-10,10];
domain_y = domain_x./aspect_ratio;

h_figure = figure;

x_start = repelem(0,20);
y_start = repelem(0,20);

hold on
grid on
axis equal
h_plot = scatter(x_start, y_start, "r","filled");
xlim([view_domain_x])
ylim([view_domain_y])
ax = gca;
ax.FontSize = 20;

set(h_plot, 'ButtonDownFcn', {@Mouse_Callback, 'down'});
drawnow()
set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

setappdata(h_figure, 'coords_updated', false);

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
colormap(flip(cmap))

ind_1=1;
while isvalid(h_figure)
    pause(0.001)
    if getappdata(h_figure, 'coords_updated')
        coords_new = getappdata(h_figure, 'updated_coords');
        roots_new = (coords_new(:,1) + coords_new(:,2).*j).';
        setappdata(h_figure, 'coords_updated', false);

        [root_result, iterations_result] = build_fractal(roots_new, canvas_x, canvas_y, domain_x, domain_y, convergence_tolerance, max_iters);

        if ind_1>1
            delete(h_surf);
        end
        h_surf = surf(canvas_x, canvas_y, root_result-max(max(root_result))-2,EdgeColor="none");
        %set(gca,"colorscale","log")
        drawnow()

        ind_1 = ind_1+1;

    end
end



function Mouse_Callback(hObj,~,action)
    persistent curobj xdata ydata ind
    pos = get(gca,'CurrentPoint');
    switch action
      case 'down'
          curobj = hObj;
          xdata = get(hObj,'xdata');
          ydata = get(hObj,'ydata');
          [~,ind] = min(sum((xdata-pos(1)).^2+(ydata-pos(3)).^2,1));
          set(gcf,...
              'WindowButtonMotionFcn',  {@Mouse_Callback,'move'},...
              'WindowButtonUpFcn',      {@Mouse_Callback,'up'});
      case 'move'
          xdata(ind) = pos(1);
          set(curobj,'xdata',xdata)
          ydata(ind) = pos(3);
          set(curobj,'ydata',ydata)

        coords = [xdata(:), ydata(:)];
        setappdata(gcf, 'updated_coords', coords);
        setappdata(gcf, 'coords_updated', true);
      case 'up'
          set(gcf,...
              'WindowButtonMotionFcn',  '',...
              'WindowButtonUpFcn',      '');
    end
end


function [root_result, iterations_result] = build_fractal(roots, canvas_x, canvas_y, domain_x, domain_y, convergence_tolerance, max_iters)

    roots = snap_close_points(roots, 0.1);

    f_poly = poly(roots);
    fp_poly = polyder(f_poly);
    
    for n=1:length(canvas_x)
        for m=1:length(canvas_y)
            start_location = canvas_x(n) + canvas_y(m)*j;
            [converged_root, conv_iterations] = NR_iter(start_location, roots, f_poly, fp_poly, convergence_tolerance,max_iters);
            iterations_result(m,n) = conv_iterations;
            root_result(m,n) = converged_root;
        end
    end
    
    iterations_result = iterations_result./max(max(iterations_result));
end

function [converged_root, d_conv] = NR_iter(start_location, roots, f_poly, fp_poly, convergence_tolerance, max_iters)
    conv_iterations = 0;
    dx = inf;
    xn = start_location;
    while dx > convergence_tolerance && conv_iterations < max_iters
        plane_np = xn - polyval(f_poly,xn)/polyval(fp_poly,xn);
        plane_np =  sinh(plane_np);
        dx = abs(plane_np - xn);
        xn = plane_np;
        conv_iterations = conv_iterations+1;
    end
    %which root is closest?
    [~, converged_root] = min(abs(plane_np - roots));
    d_conv = abs(xn - start_location)/conv_iterations;
end
