format compact
clear
clc
%close all
clf reset

load("best_paramlist.mat")

write_to_img = true;

iters = 100;
blowup_cutoff = 1e30;

res = [1080,1920]./1;
dims_x = [-2,2];

aspect_ratio = res(2)/res(1);
dims_y = dims_x.*aspect_ratio;

x_range = single( linspace(dims_x(1), dims_x(2), res(1)) );
y_range = single( linspace(dims_y(1), dims_y(2), res(2)) );

[x_plane, y_plane] = meshgrid(y_range, x_range);
c_plane = x_plane + y_plane.*j; 

parameters=[
global_best_params
];

dconv_total = single( zeros(size(c_plane)) );

frac_index_cumulative = 0;
for ind_frac = 1:fractal_population

    a = parameters(ind_frac+frac_index_cumulative,1) + parameters(ind_frac+frac_index_cumulative,2)*j;
    c = parameters(ind_frac+frac_index_cumulative+1,1) + parameters(ind_frac+frac_index_cumulative+1,2)*j;
    t1 = parameters(ind_frac+frac_index_cumulative+2,1) + parameters(ind_frac+frac_index_cumulative+2,2)*j;
    t2 = parameters(ind_frac+frac_index_cumulative+3,1) + parameters(ind_frac+frac_index_cumulative+3,2)*j;
    
    frac_index_cumulative = frac_index_cumulative+3;

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

    dconv_total = dconv_total + d_conv;
end
dconv_total = flipud(dconv_total);
finite_dconv = dconv_total(~isinf(dconv_total));
dconv_total(isinf(dconv_total)) = max(max(finite_dconv));

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

plot_field = dconv_total;
dot_mask = get_dotmask(plot_field);
plot_field(dot_mask)=0;

frame_colourised = colour_frame(cmap,plot_field);

if write_to_img
    imwrite(cast(flipud(frame_colourised),'uint8'),'fractal_approx.png');
end

hold on
imagesc(frame_colourised)
axis tight equal
colormap(cmap)

function dot_mask = get_dotmask(plot_field)
    logical_plotfield = logical(plot_field);
    scape_vert_1 = circshift(logical_plotfield,[1,0]);
    scape_hor_1 = circshift(logical_plotfield,[0,1]);
    scape_vert_2 = circshift(logical_plotfield,[-1,0]);
    scape_hor_2 = circshift(logical_plotfield,[0,-1]);
    a = xor(logical_plotfield,scape_vert_1);
    b = xor(logical_plotfield,scape_hor_1);
    c = xor(logical_plotfield,scape_vert_2);
    d = xor(logical_plotfield,scape_hor_2);
    edge_mask = logical(a+b+c+d);
    dot_mask = imerode(edge_mask,strel('disk', 1));
end

function return_frame = colour_frame(cmap,dconv_total)
    plot_field = dconv_total./max(max(dconv_total));
    %logscale_map = flip( round(height(cmap) - logspace(log10(1),log10(height(cmap)), height(cmap) )) + 1 );
    plot_field = ceil( plot_field.*height(cmap) );
    plot_field(plot_field==0) = 1;
    for n=1:3
        %cmap_channel_spec = cmap(logscale_map(plot_field),n) .* 255;
        cmap_channel_spec = cmap(plot_field,n) .* 255;
        field_colour_spec = reshape(cmap_channel_spec,[size(plot_field)]);
        return_frame(:,:,n) = uint8(field_colour_spec);
    end
end
