format compact
clear
clc
%close all
clf reset


filename = "best_paramlist.mat";
if ~exist(filename, 'file')
    iter = 1;
    save(filename,"iter")
end

subparameters = 4; %how many numbers control each fractal

fractal_iters = 100;

fractal_population = 8; %how many fractals are we layering ontop of each other

[image_raw] = imread("image.jpg");
image_raw_comp = uint8( zeros(height(image_raw),width(image_raw)) );
for n=1:3
    image_raw_comp = image_raw_comp + squeeze(image_raw(:,:,n));
end
image_raw_comp = single(image_raw_comp);
image_raw_comp = image_raw_comp./max(max(image_raw_comp));
image_target = single(zeros(size(image_raw_comp)));
image_target(image_raw_comp>0.5)=1;
%image_target = flipud(image_target);

res = [1080,1920]./30;
dims_x = [-2,2];

aspect_ratio = res(2)/res(1);
dims_y = dims_x.*aspect_ratio;

x_range = single( linspace(dims_x(1), dims_x(2), res(1)) );
y_range = single( linspace(dims_y(1), dims_y(2), res(2)) );

[x_plane, y_plane] = meshgrid(y_range, x_range);
c_plane = x_plane + y_plane.*j; 

target_small = imresize(image_target,res-1);

parameters_inital = generate_inital(fractal_population,subparameters,c_plane);
parameter_list = parameters_inital;

% dconv_total = run_multifractal(fractal_population, parameter_list, c_plane, fractal_iters);
% change_single_param(parameters_inital,80)

param_search_tree = struct();
param_fitness_tree = struct();

iter = 1;
iter_bitlength_prev = 0;
branch_pointer = iter_bitlength_prev;
best_fit_current = -inf;
global_best_fit = -inf;
stepforward = false;

ind_currentbest_history = 1;

extrascans_indexes = single( round(linspace(1,numel(parameter_list),length(parameter_list)*2)) );

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

while true

    iter_bitlength = length(dec2bin(iter));

    if iter_bitlength ~= iter_bitlength_prev
        %add a new field to struct
        param_search_tree = setfield(param_search_tree,'iter_'+string(iter_bitlength),parameter_list);
        param_fitness_tree = setfield(param_fitness_tree,'iter_'+string(iter_bitlength),best_fit_current);
        
        branch_pointer = branch_pointer + 1;
        if branch_pointer > 7
            branch_pointer = 7;
            iter_bitlength = 7;
        end

        if ~stepforward && branch_pointer > 3
            %we need to roll back to the prev branch
            iter_bitlength = iter_bitlength - 2;
            iter = 2^(iter_bitlength);
            branch_pointer = branch_pointer - 2;

            if rand() < 1e-3
                parameter_list = global_best_params; %chance we just jump back to global best
                best_fit_current = global_best_fit;
                fprintf("restoring from global best.\n")
            else
                
                iters_available = get_tree_inds(param_search_tree,false);
                [~,ind_iterclosest] = min(abs(iters_available-iter_bitlength));
                iter_bitlength = iters_available(ind_iterclosest);
                
                parameter_list = getfield(param_search_tree,'iter_'+string(iter_bitlength) ); %get best from previous branch
                best_fit_current = getfield(param_fitness_tree,'iter_'+string(iter_bitlength) );
            end

        end

        stepforward = false;
        
        currentbest_history(ind_currentbest_history) = best_fit_current;
        ind_currentbest_history = ind_currentbest_history+1;

        if length(currentbest_history) > 40
            if std(currentbest_history(end-40:end)) < 5e-4 || rand() < 5e-4
                parameter_list = generate_inital(fractal_population,subparameters,c_plane); %start again
                currentbest_history = [];
                ind_currentbest_history = 1;
                iter_bitlength = 1;
                iter = 1;
                branch_pointer = 1;
                fprintf("rolling back to start.\n")
            end
        end

        %new parameters for the testing on this branch
        iter_stepfactor = round(interp1([0,1],[1,100],rand())); %how long will this branch test last? larger is faster
        accept_anyway_prob = interp1([0,1],[1e-3,0.2],rand()); %likelihood that we accept a bad configuration anyway
        if rand() < 0.25
            accept_anyway_prob = 0;
        end
        if rand() < 0.01
            accept_anyway_prob = 0.9;
        end

        fprintf("- branch pointer: %i, current best fitness: %3.3f, global best fitness: %3.3f.\n",branch_pointer, best_fit_current, global_best_fit)
    end

    iter_bitlength_prev = iter_bitlength;

    %fractal testing
    for ind_change = 1:length(extrascans_indexes)

        parameter_list_trial = change_single_param(parameter_list,extrascans_indexes(ind_change));
        
        dconv_total = run_multifractal(fractal_population, parameter_list_trial, c_plane, fractal_iters);
        fitness = get_fitness(dconv_total,target_small);

        dconv_linear = reshape(dconv_total,1,[]);
        if all(isfinite(dconv_linear)) && sum(sum(dconv_total==0)) > 10
            exclude = false;
        else
            exclude = true;
        end

        update_plot = false;
        if fitness > best_fit_current && ~exclude
            %accept new parameter list
            parameter_list = parameter_list_trial;
            best_fit_current = fitness;

            param_search_tree = setfield(param_search_tree,'iter_'+string(iter_bitlength),parameter_list);
            param_fitness_tree = setfield(param_fitness_tree,'iter_'+string(iter_bitlength),best_fit_current);

            stepforward = true;

            if fitness > global_best_fit

                global_best_fit = fitness;
                global_best_params = parameter_list_trial;
                global_best_dconv = dconv_total;
                update_plot = true;
            end

        end
        if rand() < accept_anyway_prob && ~exclude
            %accept new parameter list anyway but don't write it.
            parameter_list = parameter_list_trial;
        end

        if update_plot

            colormap(cmap)
            subplot(2,2,1)
            dconv_total = run_multifractal(fractal_population, global_best_params, c_plane, 30);
            imagesc(global_best_dconv)
            axis tight equal

            subplot(2,2,2)
            imagesc(global_best_dconv==0)
            axis tight equal

            subplot(2,2,3)
            imagesc(image_target)
            axis tight equal

            drawnow()

            save(filename,"iter","param_fitness_tree","param_search_tree","global_best_params","subparameters","fractal_population","image_target") %save because why not
        end
    end

    iter = iter + iter_stepfactor;
end




function a = generate_inital(fractal_population,subparameters,c_plane)
    break_good = false;
    while ~break_good
        a = single( rand(fractal_population*subparameters,2).*2 );
        dconv_total = run_multifractal(fractal_population, a, c_plane, 10);
        if sum(sum(dconv_total==0)) > 20 && all(all(isfinite(dconv_total)))
            break_good = true;
        end
    end
end

function fitness = get_fitness(b,target_small)
    b_set = single( b==0 );
    b_set = imresize(b_set,size(target_small));
    fitness = 1/rmse(b_set,target_small, "all");
end


function a = change_single_param(b,ind)
    [col,row] = ind2sub(flip(size(b)),ind);
    a = b;
    a(row,col) = a(row,col) + sign(rand()-0.5) * rand()*10 * 10^interp1([0,1],[-7,2],rand());
end

function a = get_tree_inds(b,get_end)
    c = fieldnames(b);
    c = string(c);
    a = uint16(str2double(erase(c,"iter_")));
    if get_end
        a = a(end);
    end
end
