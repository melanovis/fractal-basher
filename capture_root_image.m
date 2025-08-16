format compact
clear
clc
clf reset
%close all


filename = "root_catalog.mat";
load(filename)

roots = root_map(find(root_names=="newton"),:);
image_raw = imread("newton.jpg");

img_ontop = false;

convergence_tolerance = 1e-7;
max_iters = 2e2;

dims = ceil([1920,1080]./8); %control res
aspect_ratio = dims(1)/dims(2);
view_domain_x = [-1,1]; 
view_domain_y = view_domain_x./aspect_ratio;
canvas_x = linspace(view_domain_x(1),view_domain_x(2),dims(1));
canvas_y = linspace(view_domain_y(1),view_domain_y(2),dims(2));

image_raw = squeeze(image_raw(:,:,1));
image_raw = image_raw./max(max(image_raw));
image_raw = imresize(image_raw,[flip(dims)]);
image_raw = abs(image_raw);
target_image = round(double(image_raw));

% canvas_x = 3.5*(canvas_x*0.67+0.015); %some extra touch-ups
% canvas_y = 3.5*(canvas_y + 0.1);

domain_x = [-1,1];
domain_y = domain_x./aspect_ratio;

[x_plane, y_plane] = meshgrid(canvas_x, canvas_y);
complex_plane = x_plane + y_plane.*j; 

root_quantity = 40;

bounds = 1;
[root_grid_default_x,root_grid_default_y] = meshgrid(linspace(-bounds,bounds,floor(sqrt(root_quantity))), linspace(-bounds,bounds,floor(sqrt(root_quantity))));
root_grid_default = root_grid_default_x + root_grid_default_y.*j;
root_grid_default = reshape(root_grid_default,1,[]);
while length(root_grid_default) < root_quantity
    root_grid_default = [root_grid_default,0];
end

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

[converged_root, d_conv] = NR_iter_frame(complex_plane, roots, convergence_tolerance, max_iters);
converged_root = converged_root-1;

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

maskout_index = 12;

rootmap_full_maskout = zeros(size(converged_root));
for n=1:height(converged_root)
    rootmap_with_maskout = [converged_root(n,:) == maskout_index].';
    d_conv(n,rootmap_with_maskout) = nan;
    rootmap_full_maskout(n,:) = rootmap_with_maskout;
end


colormap(flip(cmap))
%colormap(cmap)

hold on
grid on
axis vis3d equal
axes('Units', 'normalized', 'Position', [0 0 1 1]) 
xlim(domain_x)
ylim(domain_y)
view([0,90])
d_conv = flip(d_conv);

if img_ontop
    diffmask = target_image.*max(max(d_conv)) - d_conv;
    diffmask = abs(diffmask./max(max(diffmask)));
    imagesc(canvas_x, canvas_y,diffmask)
else
    imagesc(canvas_x, canvas_y, d_conv.^0.75, alphadata = ~isnan(d_conv))
    clim([1,max_iters]);
end

set(gca,'Color','k')
set(gca,'TickLength',[0 0])
%set(gca,"ColorScale","log")

sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);






