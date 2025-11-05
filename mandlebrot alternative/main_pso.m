format compact
clear
clc
%close all
clf reset

%-----

[image_raw] = imread("image.jpg");
image_raw_comp = uint8( zeros(height(image_raw),width(image_raw)) );
for n=1:3
    image_raw_comp = image_raw_comp + squeeze(image_raw(:,:,n));
end
image_raw_comp = single(image_raw_comp);
image_raw_comp = image_raw_comp./max(max(image_raw_comp));
image_target = single(zeros(size(image_raw_comp)));
image_target(image_raw_comp>0.5)=1;
image_target = flipud(image_target);


cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
colormap(cmap)

res = [1080,1920]./19;
dims_x = [-1,1];

%input_vector = rand(1,16);
%[fitness,output_vector] = mandle_fitness(res, dims_x, input_vector,image_target)

n_DOF = 16;
var_size = [1, n_DOF]; %solution matrix size
var_min = 0;
var_max = 1;

%clerc kennedy construction coefficient function
phi_1 = 2.05;
phi_2 = 2.05;
phi = phi_1+phi_2;
kappa = 1;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

population = 8*50;
w = chi;
w_damp = 0.9;
w_original = w_damp;
c1 = chi*phi_1;
c2 = chi*phi_2;
max_velocity = (var_max-var_min)*0.2;
min_velocity = -max_velocity;

%initalise particle template
empty_particle.position = [];
empty_particle.velocity = [];
empty_particle.fitness = [];
empty_particle.best.position = [];
empty_particle.best.fitness = [];

%initalise global best
global_best.fitness = -inf;

%create population matrix
particle = repmat(empty_particle, population, 1);

if ~gcp().Connected
    delete(gcp('nocreate'));
    parpool('local',8);
end

parfor n=1:population
    particle(n).position = unifrnd(var_min,var_max,var_size); %uniform distribution is used for initalisation

    input_vector = particle(n).position;
    [fitness,output_vector,d_conv] = mandle_fitness(res, dims_x, input_vector,image_target);
    particle(n).fitness = fitness;

    particle(n).velocity = zeros(var_size); 
    
    particle(n).best.position = particle(n).position;
    particle(n).best.fitness = particle(n).fitness;
end

for n=1:population
    if particle(n).best.fitness > global_best.fitness
        global_best = particle(n).best;
    end
end

fprintf("-----------\n")

iter = 1;
best_fit_prev = -inf;
ag_countdown = 5;
global_best_fitness = -inf;
steps_forward = 0;
return_history = [];
seq_history = [];
return_lengths_history = [];
best_return_params = nan;
bit_proximity = nan;

while global_best_fitness < 1e9
    parfor n=1:population
        
        particle(n).velocity = w*particle(n).velocity ...
            + c1*rand(var_size).*(particle(n).best.position - particle(n).position) ...
            + c2*rand(var_size).*(global_best.position - particle(n).position);
        
        particle(n).velocity = max(particle(n).velocity, min_velocity);
        particle(n).velocity = min(particle(n).velocity, max_velocity);

        particle(n).position = particle(n).position + particle(n).velocity;
        
        particle(n).position = max(particle(n).position, var_min);
        particle(n).position = min(particle(n).position, var_max);

        input_vector = particle(n).position;
        [fitness,output_vector,d_conv] = mandle_fitness(res, dims_x, input_vector,image_target);
        particle(n).fitness = fitness;

        %update personal best
        if particle(n).fitness > particle(n).best.fitness
            particle(n).best.position = particle(n).position;
            particle(n).best.fitness = particle(n).fitness;
        end
        
        %summarise in matrix
        position_summary(n,:) = particle(n).position;
    end

    update_plot = false;
    for n=1:population
        if particle(n).best.fitness > global_best.fitness
            global_best = particle(n).best;

            if global_best.fitness > best_fit_prev
                global_best_fitness = global_best.fitness;
                best_fit_prev = global_best.fitness;
                %update global best
                steps_forward = steps_forward+1;

                input_vector = particle(n).position;
                [fitness_best,output_vector_best,d_conv_best] = mandle_fitness([1080,1920]./5, dims_x, input_vector,image_target);
                update_plot = true;
            end
        end
    end

    if update_plot
        scatter(nan,nan)
        subplot(1,2,1)
        hold on
        grid on
        imagesc(d_conv_best)
        axis tight equal

        subplot(1,2,2)
        hold on
        grid on
        imagesc(image_target)
        axis tight equal
        drawnow()
        hold off
    end


    bestfitnesss(iter) = global_best.fitness;

    if iter > 10
        if std(bestfitnesss(end-3:end)) < 5e-3
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
    
            ind_ag = ind_sort(1:round(length(ind_sort)*0.9));

            for n=1:length(ind_ag)
                particle(ind_ag(n)).position = rand(1,n_DOF)*var_max;
                particle(ind_ag(n)).velocity = rand(1,n_DOF)*max_velocity;
            end
            
            for n=1:population
                if ~ismember(n, ind_sort(end-3:end))
                    particle(n).position = clamp(particle(n).position + (rand(1,n_DOF)-0.5).*1e-4, 0,1);
                    particle(n).velocity = clamp(particle(n).velocity + (rand(1,n_DOF)-0.5).*1e-4, min_velocity,max_velocity);
                end
            end

            w = w_original;
            ag_countdown = 10;
            fprintf("\n agitating.\n")
        end
    end

    w = w * w_damp;

    fprintf("completed iter %i, best fitness: %3.6f.\n",iter,bestfitnesss(iter))

    iter = iter+1;
end

function b = clamp(a,l,u)
    a(a<l) = l;
    a(a>u) = u;
    b=a;
end