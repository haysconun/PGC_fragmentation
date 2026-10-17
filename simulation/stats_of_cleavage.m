function [randiiiii, div_asymm_NC9_10, div_asymm_NC10_11, cv_plasm, real_bud_counter_cat, NC9_count] = ...
    stats_of_cleavage(A, a_val, num_points_pm, dist_cut, ...
    pct_cutoff, NC11_flag, num_realizations)
% this function just generates different outputs of generate_hexagon_tesselation 
% and computes some statistics about them
% NOTE: this never compute solutions of the growing ODE.

% OUTPUTS:
% randiiiii
% div_asymm_NC9_10
% div_asymm_NC10_11
% cv_plasm
% real_bud_counter_cat
% NC9_count

% INPUTS:
% A
% a_val
% num_points_pm
% dist_cut
% pct_cutoff
% NC11_flag
% num_realizations

% A = 1 / 2 / pi
% a_val = 1.1
% num_points_pm = 5
% dist_cut_val = 14
% pct_cutoff = 0.015
% NC11_flag = 1
% num_realizations = 5

% DO NOT PLOT -- this would quickly become a nuisance
plot_NC9 = false;
plot_NC10 = false;
plot_NC11 = false;
should_i_plot = false;

% max seed val
maxval = 10000000;
randiiiii = randi(maxval, 2, num_realizations);
assert(isequal(size(randiiiii),[2,num_realizations]));

% store count of "clones" at NC 9
NC9_count = [];
% compute division asymmetries from NC 9 into NC 10
div_asymm_NC9_10 = [];
% compute division asymmetries from NC 10 into NC 11
div_asymm_NC10_11 = [];
% store cv of integrated plasm distribution at NC 11
cv_plasm = [];
% count the number of cells at NC 11
real_bud_counter_cat = [];

for iii = 1:num_realizations
    % seed with random number for rnd_mag, rnd_angle
    rng(randiiiii(1,iii));

    u = rand;  v = rand;
    offset = u*a_val*[0,1] + v*a_val*[sind(120), cosd(120)];
    rnd_mag   = norm(offset);
    rnd_angle = atan2(offset(2), offset(1));

    % seed with random number for sampling division orientations.
    rand_seed = randiiiii(2,iii);
    % run generate_hexagon_tesselation.
    [ ~, ~, ~, ~, ...
        ~, real_bud_counter, ~, ~, ...
        cell_output_NC_9_10_11, cell_output_NC_9_10_11_gp, mean_points ] = ...
        generate_hexagon_tesselation(A, rnd_mag, rnd_angle, a_val, num_points_pm, dist_cut, ...
        pct_cutoff, rand_seed, NC11_flag, should_i_plot, plot_NC9, plot_NC10, plot_NC11);

    % store % plasm for cellularizing guys at NC 11
    plasm_vect = [];
    for jjj = 1:size(cell_output_NC_9_10_11_gp, 1)
        % compute normalized division asymmetry from NC 9 into NC 10
        this_asymm = abs(cell_output_NC_9_10_11_gp{jjj,2}(2)-cell_output_NC_9_10_11_gp{jjj,2}(1));
        this_asymm = this_asymm / sum(cell_output_NC_9_10_11_gp{jjj,2});
        div_asymm_NC9_10 = vertcat(div_asymm_NC9_10, this_asymm);

        % store vector of daughters at NC 11
        vect1 = cell_output_NC_9_10_11_gp{jjj,3}(1:2);
        vect2 = cell_output_NC_9_10_11_gp{jjj,3}(3:4);

        % check if any are over cutoff for cellularization
        if (any(vect1 >= pct_cutoff))
            % compute normalized division asymmetry
            this_asymm = abs(vect1(2)-vect1(1));
            this_asymm = this_asymm / sum(vect1);
            div_asymm_NC10_11 = vertcat(div_asymm_NC10_11, this_asymm);
            % store plasm percentage for those over cutoff
            plasm_vect = vertcat(plasm_vect, vect1(find(vect1>=pct_cutoff)));
        end
        % check if any are over cutoff for cellularization
        if (any(vect2 >= pct_cutoff))
            % compute normalized division asymmetry
            this_asymm = abs(vect2(2)-vect2(1));
            this_asymm = this_asymm / sum(vect2);
            div_asymm_NC10_11 = vertcat(div_asymm_NC10_11, this_asymm);
            % store plasm percentage for those over cutoff
            plasm_vect = vertcat(plasm_vect, vect2(find(vect2>=pct_cutoff)));
        end
    end
    assert(real_bud_counter == numel(plasm_vect), ...
        'Cellularizing buds outside dist_cut -- increase dist_cut');

    % store cv of distribution at NC 11
    cv_plasm = vertcat(cv_plasm, std(plasm_vect)/mean(plasm_vect));
    % counter number of cells at NC 11
    real_bud_counter_cat = vertcat(real_bud_counter_cat, real_bud_counter);
    % count number of clones.
    NC9_count = vertcat(NC9_count, size(cell_output_NC_9_10_11_gp,1));

    disp(iii);
end

end