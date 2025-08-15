format compact
clear
clc
clf reset
%close all

image_raw = imread("newton.jpg");
filename = "root_catalog.mat";
root_config_name = "newton";

if ~exist(filename, 'file')
    root_names = [];
    root_map = [];
    save(filename,"root_names","root_map");
end
load(filename)
root_names = unique([root_names; root_config_name]);
save(filename,"root_names","root_map");

convergence_tolerance = 1e-7;
max_iters = 20;

dims = ceil([1920,1080]./30);
aspect_ratio = dims(1)/dims(2);
view_domain_x = [-1,1];
view_domain_y = view_domain_x./aspect_ratio;
canvas_x = linspace(view_domain_x(1),view_domain_x(2),dims(1));
canvas_y = linspace(view_domain_y(1),view_domain_y(2),dims(2));

domain_x = [-1,1];
domain_y = domain_x./aspect_ratio;

image_raw = flip(single(image_raw));
image_raw = squeeze(image_raw(:,:,1));
image_raw = image_raw./max(max(image_raw));
%image_raw = 1-image_raw;
image_raw = imresize(image_raw,[flip(dims)]);
image_raw = abs(image_raw);
target_image = round(image_raw);
%target_image = imdilate(target_image, strel('disk', 1));

root_quantity = 40+3;

n_DOF = root_quantity*4; %how many variables?
var_size = [1, n_DOF]; %solution matrix size
var_min = 0; %variable bounds
var_max = 1;

%clerc kennedy construction coefficient function
phi_1 = 2.05;
phi_2 = 2.05;
phi = phi_1+phi_2;
kappa = 1;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

population = 8*100;
w = chi;
w_damp = 0.7;
w_original = w_damp;
c1 = chi*phi_1;
c2 = chi*phi_2;
max_velocity = (var_max-var_min)*0.25;
min_velocity = -max_velocity;

%initalise particle template
empty_particle.position = [];
empty_particle.velocity = [];
empty_particle.fitness = [];
empty_particle.best.position = [];
empty_particle.best.fitness = [];

%initalise global best
globalbest.fitness = -inf;

%create population matrix
particle = repmat(empty_particle, population, 1);

[x_plane, y_plane] = meshgrid(canvas_x, canvas_y);
complex_plane = x_plane + y_plane.*j; 

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

delete(gcp('nocreate'));
parpool('local',8);

parfor n=1:population

    fitness_out = 0;
    while fitness_out < 2 %we'd like to be starting off with at least half decent options
        particle(n).position = double(unifrnd(var_min,var_max,var_size)); %uniform distribution is used for initalisation
        [converged_root, return_matrix, roots_return, fitness_out] = NR_iter_PSO(complex_plane, particle(n).position, convergence_tolerance, max_iters, target_image);
        particle(n).fitness = fitness_out;
    end
    fprintf("-")

    particle(n).velocity = zeros(var_size); 
    
    particle(n).best.position = particle(n).position;
    particle(n).best.fitness = particle(n).fitness;
end

for n=1:population
    if particle(n).best.fitness > globalbest.fitness
        globalbest = particle(n).best;
        [converged_root, return_matrix, roots_return, fitness_out] = NR_iter_PSO(complex_plane, particle(n).position, convergence_tolerance, max_iters, target_image);
        global_best_img = return_matrix;
    end
end

fprintf("------\n")

stop_fitness = 1e3;

iter = 1;
best_fit_prev = -inf;
ag_countdown = 5;
global_best_fitness = -inf;
best_roothistory = [];
ind_bestroots = 1;

while global_best_fitness < stop_fitness
    parfor n=1:population
        
        particle(n).velocity = w*particle(n).velocity ...
            + c1*rand(var_size).*(particle(n).best.position - particle(n).position) ...
            + c2*rand(var_size).*(globalbest.position - particle(n).position);
        
        particle(n).velocity = max(particle(n).velocity, min_velocity);
        particle(n).velocity = min(particle(n).velocity, max_velocity);

        particle(n).position = particle(n).position + particle(n).velocity;
        
        particle(n).position = max(particle(n).position, var_min);
        particle(n).position = min(particle(n).position, var_max);

        [converged_root, return_matrix, roots_return, fitness_out] = NR_iter_PSO(complex_plane, particle(n).position, convergence_tolerance, max_iters, target_image);
        particle(n).fitness = fitness_out;

        %update personal best
        if particle(n).fitness > particle(n).best.fitness
            particle(n).best.position = particle(n).position;
            particle(n).best.fitness = particle(n).fitness;
        end
        
        %summarise in matrix
        position_summary(n,:) = particle(n).position;
    end

    update_plot = false;
    update_root_map = false;
    for n=1:population
        if particle(n).best.fitness > globalbest.fitness
            globalbest = particle(n).best;

            if globalbest.fitness > best_fit_prev
                global_best_fitness = globalbest.fitness;
                best_fit_prev = globalbest.fitness;
                update_plot = true;
                update_root_map = true;

                [converged_root, return_matrix, roots_return, fitness_out] = NR_iter_PSO(complex_plane, particle(n).position, convergence_tolerance, max_iters, target_image);
                
                global_best_img = return_matrix;
                best_roothistory(ind_bestroots,1:length(roots_return)) = roots_return;
                ind_bestroots = ind_bestroots + 1;
            end
        end
    end
    if iter==1
        update_plot = true;
        best_roothistory(ind_bestroots,1:length(roots_return)) = roots_return;
        ind_bestroots = ind_bestroots + 1;
    end

    best_fitness(iter) = globalbest.fitness;

    if iter > 10
        if std(best_fitness(end-3:end)) < 5e-3
            ag_countdown = ag_countdown-1;
        else
            ag_countdown = 10;
        end
        if ag_countdown <= 0
            fitness_score = [];
            for n=1:population
                fitness_score(n) = particle(n).fitness;
            end
            [~,ind_sort] = sort(fitness_score);

            if rand()<0.75
                ind_ag = ind_sort(1:round(length(ind_sort) * interp1([0,1],[0.3,0.999],rand())) );
            else
                ind_ag = ind_sort(end - length(ind_sort)*0.05 : end-1);
            end

            for n=1:length(ind_ag)
                particle(ind_ag(n)).position = rand(1,n_DOF)*var_max;
                particle(ind_ag(n)).velocity = rand(1,n_DOF)*max_velocity;
            end
            w = w_original;
            ag_countdown = 10;
            fprintf("\n agitating.\n")
        end
    end

    w = w * w_damp;

    fprintf("completed %i, steps forward: %i, best fitness: %3.6f.\n" ,iter, ind_bestroots, best_fitness(iter))

    if update_plot

        root_map = squeeze(global_best_img(:,:,1));
        
        colormap(cmap)

        subplot(2,2,1)
        scatter(nan,nan)
        hold on
        grid on
        axis tight equal
        surf(canvas_x, canvas_y, root_map,EdgeColor="none")
        %set(gca,"colorscale","log")
        hold off

        conv_map = squeeze(global_best_img(:,:,2));
        subplot(2,2,2)
        scatter(nan,nan)
        hold on
        grid on
        axis tight equal
        surf(canvas_x, canvas_y, conv_map,EdgeColor="none")
        %set(gca,"colorscale","log")
        hold off

        subplot(2,2,3)
        scatter(nan,nan)
        hold on
        grid on
        axis tight equal
        surf(canvas_x, canvas_y, target_image,EdgeColor="none")
        %set(gca,"colorscale","log")
        hold off

        mask_map = squeeze(global_best_img(:,:,3));
        subplot(2,2,4)
        scatter(nan,nan)
        hold on
        grid on
        axis tight equal
        surf(canvas_x, canvas_y, mask_map,EdgeColor="none")
        %set(gca,"colorscale","log")
        hold off
        
        set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

        drawnow()        
    end

    if update_root_map 
        load(filename)
        root_catind = find(root_config_name == root_names);
        global_best_roots = best_roothistory(end,:);
        root_map(root_catind,1:length(global_best_roots)) = global_best_roots;
        save(filename,"root_names","root_map");
    end

    iter = iter+1;
end


