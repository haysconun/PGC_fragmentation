function [total_cell_count, t_val, y_val, which_cycle_soma, ...
    sorted_cell_output_NC_9_10_11_gp, sorted_cell_output_NC_9_10_11, ...
    div_counter, extra_div_counter, store_index, which_above, which_osc_ind_vect, ...
    gp_vect, cv_9, cv_10, cv_11, cv_14, real_NC9_counter, real_NC10_counter, ...
    real_bud_counter] = sim_examples(iii, randiiiii, color_vect, assert_check_bool, NC11flag, test_case)
%% OUTPUTS:
% total_cell_count
% t_val
% y_val
% which_cycle_soma
% sorted_cell_output_NC_9_10_11_gp
% sorted_cell_output_NC_9_10_11
% div_counter
% extra_div_counter
% store_index
% which_above
% which_osc_ind_vect
% gp_vect
% cv_9
% cv_10
% cv_11
% cv_14
% real_NC9_counter
% real_NC10_counter
% real_bud_counter

%% INPUTS:
% iii
% randiiiii
% color_vect
% assert_check_bool
% NC11flag

% test_case = 1 -> fig. 6 example 1

% test_case = 2 -> fig. 6 example 2

% test_case = 3 -> fig. 5

% test_case = 4
% everything else

a_val = 1.1;
if (assert_check_bool)
    A = 1 / 2 / pi;
else
    A = 1 / 2 / pi;
    A = A / 2;
end

rng(randiiiii(1,iii));

if (test_case == 1)
    rnd_mag = 0.1886;
    rnd_angle = 0.8124;
elseif (test_case == 2)
    rnd_mag = 0.4541;
    rnd_angle = 1.7177;
elseif (test_case == 3)
    rnd_mag = 0.2532;
    rnd_angle = 1.0949;
else
    u = rand;  v = rand;
    offset = u*a_val*[0,1] + v*a_val*[sind(120), cosd(120)];
    rnd_mag   = norm(offset);
    rnd_angle = atan2(offset(2), offset(1));
end

if (test_case == 1)
    rand_seed = 9385969;
elseif (test_case == 2)
    rand_seed = 4164252;
elseif (test_case == 3)
    rand_seed = 2020507;
else
    rand_seed = randiiiii(2,iii);
end


% sample generate_hexagon_tesselation
[ lattice_pos, pairs, percent_germ_plasm, percent_germ_plasm_half, ...
    ~, ~, above_thresh_pairs, above_thresh_pairs_dist_cutoff, ...
    cell_output_NC_9_10_11, cell_output_NC_9_10_11_gp, mean_points ] = ...
    generate_hexagon_tesselation(A, rnd_mag, rnd_angle, a_val, 5, 14, ...
    0, rand_seed, NC11flag, false, false, false, false);

non_lattice = [];
% sum percent germ plasm associated with above thresh pairs
cumsumcum = 0;
for iiij = 1:length(above_thresh_pairs)
    if (~above_thresh_pairs_dist_cutoff(iiij))
        non_lattice = vertcat(non_lattice, lattice_pos(pairs(iiij,1),:));
        non_lattice = vertcat(non_lattice, lattice_pos(pairs(iiij,2),:));
    else
        cumsumcum = cumsumcum + percent_germ_plasm(pairs(iiij,1));
        cumsumcum = cumsumcum + percent_germ_plasm(pairs(iiij,2));
    end
end

real_bud_counter = 0; % NC 11
real_NC9_counter = 0; % NC 9
% store distribution of germ plasm among NC 9 (clones) generating PGCs
gp9 = [];
% store distribution of germ plasm among buds generating PGCs
gp10 = [];
% store distribution of germ plasm among PGCs at NC 11
gp11 = [];
real_NC10_counter = 0; % buds
% check that all above pct_cutoff are above (parameter)
pct_cutoff = 0.015;

num_oscillators_NC9 = size(cell_output_NC_9_10_11_gp, 1);
for osc_ind = 1:num_oscillators_NC9
    % store germ plasm associated with each NC 11 daughter
    inhib_val3 = cell_output_NC_9_10_11_gp{osc_ind,3};
    % store germ plasm associated with each bud
    inhib_val2 = cell_output_NC_9_10_11_gp{osc_ind,2};

    real_bud_counter = real_bud_counter + length(find(inhib_val3>pct_cutoff));
    gp11 = vertcat(gp11, inhib_val3(find(inhib_val3>pct_cutoff)));
    if (~isempty(find(inhib_val3>pct_cutoff)))
        % counter number of clones generating PGCs
        real_NC9_counter = real_NC9_counter + 1;
        % store gp
        gp9 = vertcat(gp9, cell_output_NC_9_10_11_gp{osc_ind,1});
    end

    if (~isempty(find(inhib_val3(1:2)>pct_cutoff)))
        % counter number of buds generating PGCs
        real_NC10_counter = real_NC10_counter + 1;
        % store gp
        gp10 = vertcat(gp10, inhib_val2(1));
    end

    if (~isempty(find(inhib_val3(3:4)>pct_cutoff)))
        % counter number of buds generating PGCs
        real_NC10_counter = real_NC10_counter + 1;
        % store gp
        gp10 = vertcat(gp10, inhib_val2(2));
    end    
end
assert(isequal(length(gp11), real_bud_counter))
assert(isequal(length(gp10), real_NC10_counter))
assert(isequal(length(gp9), real_NC9_counter))
assert(nnz(percent_germ_plasm_half > pct_cutoff) == real_bud_counter, ...
    'Germline buds outside the 14 nearest clones -- increase dist_cut');

% calculate CV of germ plasm distribution at each stage.
cv_9 = std(gp9) /...
    mean(gp9);
cv_10 = std(gp10) /...
    mean(gp10);
cv_11 = std(gp11) /...
    mean(gp11);

kappa_von_mises = 200;
muval_von_mises = pi/5 + pi/10 + pi / 10;
total_time = 20*pi;
Diffuse_val = 2.0;
stepstep = 0.02;
k_time = 1.45/(4*pi*0.15); % proportionality constant in front of inhibitor build-up
theta_init = 0;
[t_val,y_val,which_cycle_soma,sorted_cell_output_NC_9_10_11_gp,sorted_cell_output_NC_9_10_11] = ...
    solve_embryonic_cell_cycle_oscillator_varying_updated(cell_output_NC_9_10_11, ...
    cell_output_NC_9_10_11_gp, kappa_von_mises, muval_von_mises, ...
    k_time, theta_init, total_time, pct_cutoff, Diffuse_val, non_lattice, stepstep);

plot_at_end = 1==1;
[div_counter, store_index, which_above, which_osc_ind_vect, gp_vect, extra_div_counter] = ...
    plot_cycle_traces(sorted_cell_output_NC_9_10_11_gp,color_vect,...
    t_val,y_val, plot_at_end, pct_cutoff, assert_check_bool);

% 1->2->4

if (assert_check_bool)
    total_cell_count = length(store_index)+length(which_above);
    assert(length(which_above) == real_bud_counter)
    assert(length(store_index) == div_counter)
    assert(extra_div_counter == 0)
else
    total_cell_count = length(which_above) + div_counter + 3 * extra_div_counter;
    assert(length(which_above) == real_bud_counter)
    assert(length(gp_vect) == total_cell_count)
end

% initialize cv_14 so that it has value even if assert_check_bool==false
cv_14 = 0;
if (assert_check_bool)
    vertcat_gp = [];
    for iiij = 1:size(sorted_cell_output_NC_9_10_11_gp,1)
        this_vec = sorted_cell_output_NC_9_10_11_gp{iiij,3};
        vertcat_gp = vertcat(vertcat_gp, this_vec([1 3 2 4]));
    end
    non_divide = setdiff(find(vertcat_gp>pct_cutoff),store_index);
    disp('check')
    disp(std(vertcat_gp(union(store_index, non_divide)))/mean(vertcat_gp(union(store_index, non_divide))) - cv_11)

    % preloaded colors
    yellow_color = '#5471b7';
    yellow_color = sscanf(strrep(yellow_color, '#', ''), '%2x%2x%2x')' / 255;
    green_color = '#009395';
    green_color = sscanf(strrep(green_color, '#', ''), '%2x%2x%2x')' / 255;
    brown_color = '#f17098';
    brown_color = sscanf(strrep(brown_color, '#', ''), '%2x%2x%2x')' / 255;
    blue_color = '#f0be6e'; %
    blue_color = sscanf(strrep(blue_color, '#', ''), '%2x%2x%2x')' / 255;
    gray_color = '#5fb4d8';
    gray_color = sscanf(strrep(gray_color, '#', ''), '%2x%2x%2x')' / 255;
    orange_color = '#74c48a';
    orange_color = sscanf(strrep(orange_color, '#', ''), '%2x%2x%2x')' / 255;

    color_vect =[gray_color;brown_color;blue_color;yellow_color;green_color;orange_color];
    colors_extended = distinguishable_colors(length(cell_output_NC_9_10_11),color_vect);
    color_vect = vertcat(color_vect, colors_extended);

    % can figure out color of dividing ones based on this --
    rng(randiiiii(1,iii));
    if (test_case == 1)
        rnd_mag = 0.1886;
        rnd_angle = 0.8124;
    elseif (test_case == 2)
        rnd_mag = 0.4541;
        rnd_angle = 1.7177;
    elseif (test_case == 3)
        rnd_mag = 0.2532;
        rnd_angle = 1.0949;
    else
        u = rand;  v = rand;
        offset = u*a_val*[0,1] + v*a_val*[sind(120), cosd(120)];
        rnd_mag   = norm(offset);
        rnd_angle = atan2(offset(2), offset(1));
    end

    if (test_case == 1)
        rand_seed = 9385969;
    elseif (test_case == 2)
        rand_seed = 4164252;
    elseif (test_case == 3)
        rand_seed = 2020507;
    else
        rand_seed = randiiiii(2,iii);
    end
    [ ~, ~, ~, ~, ...
        ~, ~, ~, ~, ...
        ~, cell_output_NC_9_10_11_gp_new, ~ ] = ...
        generate_hexagon_tesselation(A, rnd_mag, rnd_angle, a_val, 5, 14, ...
        pct_cutoff, rand_seed, NC11flag, false, false, false, false);

    % store copy for sorting again by first entry
    copy_sorted_cell_output_NC_9_10_11_gp = sorted_cell_output_NC_9_10_11_gp;
    num_oscillators_NC9_again = size(cell_output_NC_9_10_11_gp_new, 1);
    color_vect_check = [];
    for osc_ind = 1:num_oscillators_NC9
        which_one = [];
        for osc_ind_again = 1:num_oscillators_NC9_again
            if (abs(copy_sorted_cell_output_NC_9_10_11_gp{osc_ind,1} - cell_output_NC_9_10_11_gp_new{osc_ind_again,1})<1*10^-5)
                which_one = osc_ind_again;
            end
        end
        if (isempty(which_one))
            color_vect_check = vertcat(color_vect_check, [1,1,1]);
        else
            color_vect_check = vertcat(color_vect_check, color_vect(which_one,:));
        end
    end

    color_vect_gp = [];
    for iiij = 1:size(sorted_cell_output_NC_9_10_11_gp,1)
        color_vect_gp = vertcat(color_vect_gp, color_vect_check(iiij,:));
        color_vect_gp = vertcat(color_vect_gp, color_vect_check(iiij,:));
        color_vect_gp = vertcat(color_vect_gp, color_vect_check(iiij,:));
        color_vect_gp = vertcat(color_vect_gp, color_vect_check(iiij,:));
    end

    color_vect_gp_edited = [];
    vertcat_gp_edited = [];
    vertcat_gp_edited = vertcat(vertcat_gp_edited, vertcat_gp(store_index)/2);
    vertcat_gp_edited = vertcat(vertcat_gp_edited, vertcat_gp(store_index)/2);
    vertcat_gp_edited = vertcat(vertcat_gp_edited, vertcat_gp(non_divide));
    color_vect_gp_edited = vertcat(color_vect_gp_edited, color_vect_gp(store_index,:));
    color_vect_gp_edited = vertcat(color_vect_gp_edited, color_vect_gp(store_index,:));
    color_vect_gp_edited = vertcat(color_vect_gp_edited, color_vect_gp(non_divide,:));

    xxx = 1:length(vertcat_gp_edited);
    [~,IIIIII] = sort(vertcat_gp_edited);
    datadata = vertcat_gp_edited(IIIIII);
    datadata = datadata *100;

    % figure; hold all;
    % for i = 1:length(vertcat_gp_edited)
    %     h1111 = bar(xxx(i),datadata(i), 'BarWidth', 0.6); % Makes bars thinner
    %     h1111.FaceColor = color_vect_gp_edited(IIIIII(i),:);
    %     h1111.EdgeColor = color_vect_gp_edited(IIIIII(i),:);
    % end
    % ylim([0,20])
    % xticklabels([])
    % xlim([0.5,max(xxx)+0.5])
    % set(gca, 'FontSize', 25.0);

    disp('difference')
    disp(abs(sum(vertcat_gp_edited) - sum(gp11)))

    cv_14 = std(vertcat_gp_edited) /...
        mean(vertcat_gp_edited);

end

end