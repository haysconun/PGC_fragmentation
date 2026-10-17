function [ dydt ] = odefun_cell_cycle_oscillotile(t, y, cell_output_NC_9_10_11_gp, ...
    kappa_von_mises, muval_von_mises, k_time, pct_cutoff, Diffuse, ...
    sorted_cell_output_NC_9_10_11, non_lattice_pos)

%% OUTPUT:
% dydt: time derivative of theta's and a's

%% INPUTS:
% t: time
% y: state (theta's and a's)
% cell_output_NC_9_10_11_gp: precomputed partitioning of germ plasm through divisions (generate_hexagon_tesselation)
% kappa_von_mises: parameter of non-linear function f(theta) -- sigma^2 ~ 1/kappa
% muval_von_mises: parameter of non-linear function f(theta) -- cell cycle phase at which inhibition occurs
% k_time: prefactor for growth of inhibitor dependent on mRNA
% pct_cutoff: cutoff for half-bud cellularization
% DOCUMENT: Diffuse, sorted_cell_output_NC_9_10_11
% non_lattice_pos: if isempty, globally coupled

% NOTE: assumes that sorted_cell_output_NC_9_10_11 and
% cell_output_NC_9_10_11_gp are sorted/ordered in the same way

% Fprop: strength of somatic driving
Fprop = 0;
% make sure non_lattice_pos contains oscillator positions
assert(~isempty(non_lattice_pos));

%%
num_oscillators_NC9 = size(cell_output_NC_9_10_11_gp, 1);
sorted_cell_output_NC_9_10_11_gp = cell_output_NC_9_10_11_gp;
for osc_ind = 1:num_oscillators_NC9
    assert(abs(sorted_cell_output_NC_9_10_11_gp{osc_ind,1}-sum(sorted_cell_output_NC_9_10_11_gp{osc_ind,2}))<10^-10)
    this_vect = sorted_cell_output_NC_9_10_11_gp{osc_ind, 2};
    this_vect_next = sorted_cell_output_NC_9_10_11_gp{osc_ind, 3};

    tolc = @(q) max(2e-8, 2e-6*abs(q));
    assert(abs(this_vect(1)-sum(this_vect_next(1:2))) < tolc(this_vect(1)))
    assert(abs(this_vect(2)-sum(this_vect_next(3:4))) < tolc(this_vect(2)))

    assert(length(this_vect_next) == 4)
    assert(length(this_vect) == 2)
    % assert sorted
    assert(isequal(sort(this_vect), this_vect));
    % assert sorted
    assert(isequal(sort(this_vect_next(1:2)), this_vect_next(1:2)));
    % assert sorted
    assert(isequal(sort(this_vect_next(3:4)), this_vect_next(3:4)));
end

first_vect = [];
for osc_ind = 1:num_oscillators_NC9
    first_vect = vertcat(first_vect, ...
        sorted_cell_output_NC_9_10_11_gp{osc_ind,1});
end
% assert sorted
assert(isequal(sort(first_vect), first_vect))

%% CONVENTIONS
% Suppose:
% y(1) == theta^1_1
% y(2) == theta^1_2
% y(3) == theta^1_3
% y(4) == theta^1_4

% ...
% y((osc_ind-1)*4+1) == theta^(osc_ind)_1
% y((osc_ind-1)*4+2) == theta^(osc_ind)_2
% y((osc_ind-1)*4+3) == theta^(osc_ind)_3
% y((osc_ind-1)*4+4) == theta^(osc_ind)_4
% ...

% y((num_osc)*4 + (osc_ind-1)*4 + 1) == a^(osc_ind)_1
% y((num_osc)*4 + (osc_ind-1)*4 + 2) == a^(osc_ind)_2
% y((num_osc)*4 + (osc_ind-1)*4 + 3) == a^(osc_ind)_3
% y((num_osc)*4 + (osc_ind-1)*4 + 4) == a^(osc_ind)_4

%% CONVENTIONS
% store right hand side of ODEs
% dimension is number of clones from NC 9 * 4 (from 2 divisions) * 2 (from theta and a)
dydt = zeros(num_oscillators_NC9*4*2,1);
which_oscillators_on = zeros(num_oscillators_NC9*4,1);
which_oscillators_soma = zeros(num_oscillators_NC9*4,1);
oscillator_pos = zeros(num_oscillators_NC9*4,2);

% store d_theta_dt for each oscillator
% also be careful to specify existing conditions correctly (0, 2*pi, 4*pi)
for osc_ind = 1:num_oscillators_NC9
    % store germ plasm associated with each NC 11 daughter
    inhib_val3 = sorted_cell_output_NC_9_10_11_gp{osc_ind,3};

    % NC 9
    pospos1 = sorted_cell_output_NC_9_10_11{osc_ind,1};
    % NC 10
    pospos2 = sorted_cell_output_NC_9_10_11{osc_ind,2};
    % NC 11
    pospos3 = sorted_cell_output_NC_9_10_11{osc_ind,3};

    % for founding oscillator, calculate d_theta_dt
    osc_sub_ind = 1;
    % store indices for a and theta
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    % output of non-linear inhibition of cell cycle
    % non_lin_function(a_inhib, theta_inhib, kappa_von_mises, muval_von_mises)
    inhib_output = non_lin_function(y(calc_ind_a), y(calc_ind_theta), kappa_von_mises, muval_von_mises);
    % pct_cutoff check -- cut off from soma
    if (inhib_val3(1) > pct_cutoff)
        dydt(calc_ind_theta) = 1 - inhib_output + Fprop * sin(t - y(calc_ind_theta)) .* heaviside((y(calc_ind_theta)-4*pi)*-1);
    else
        dydt(calc_ind_theta) = 1 - inhib_output + Fprop * sin(t - y(calc_ind_theta));
        which_oscillators_soma(calc_ind_theta) = 1;
    end
    which_oscillators_on(calc_ind_theta) = 1;

    if (y(calc_ind_theta) >= 0 && y(calc_ind_theta) <= 2*pi)

        oscillator_pos(calc_ind_theta, :) = pospos1;

    elseif (y(calc_ind_theta) > 2*pi && y(calc_ind_theta) <= 4*pi)

        oscillator_pos(calc_ind_theta, :) = pospos2(1,:);

    elseif (y(calc_ind_theta) > 4*pi)

        oscillator_pos(calc_ind_theta, :) = pospos3(1,:);

    end

    % for first daughter of founding oscillator, calculate d_theta_dt
    osc_sub_ind = osc_sub_ind+1;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    % non_lin_function(a_inhib, theta_inhib, kappa_von_mises, muval_von_mises)
    inhib_output = non_lin_function(y(calc_ind_a), y(calc_ind_theta), kappa_von_mises, muval_von_mises);
    if (y((osc_ind-1)*4 + 1)>2*pi)
        if (inhib_val3(3) > pct_cutoff)
            dydt(calc_ind_theta) = 1 - inhib_output + Fprop * sin(t - y(calc_ind_theta)) .* heaviside((y(calc_ind_theta)-4*pi)*-1);
        else
            dydt(calc_ind_theta) = 1 - inhib_output + Fprop * sin(t - y(calc_ind_theta));
            which_oscillators_soma(calc_ind_theta) = 1;
        end
        which_oscillators_on(calc_ind_theta) = 1;
    end
    % osc 1
    if (y((osc_ind-1)*4 + 1) >= 0 && y((osc_ind-1)*4 + 1) <= 2*pi)

        oscillator_pos(calc_ind_theta, :) = pospos1;

    elseif (y(calc_ind_theta) <= 4*pi)

        oscillator_pos(calc_ind_theta, :) = pospos2(2,:);

    elseif (y(calc_ind_theta) > 4*pi)

        oscillator_pos(calc_ind_theta, :) = pospos3(3,:);

    end

    % for NC 11 daughter of founding oscillator, calculate d_theta_dt
    osc_sub_ind = osc_sub_ind+1;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    % non_lin_function(a_inhib, theta_inhib, kappa_von_mises, muval_von_mises)
    inhib_output = non_lin_function(y(calc_ind_a), y(calc_ind_theta), kappa_von_mises, muval_von_mises);
    if (y((osc_ind-1)*4 + 1)>4*pi)
        if (inhib_val3(2) > pct_cutoff)
            dydt(calc_ind_theta) = 1 - inhib_output + Fprop * sin(t - y(calc_ind_theta)) .* heaviside((y(calc_ind_theta)-4*pi)*-1);
        else
            dydt(calc_ind_theta) = 1 - inhib_output + Fprop * sin(t - y(calc_ind_theta));
            which_oscillators_soma(calc_ind_theta) = 1;
        end
        which_oscillators_on(calc_ind_theta) = 1;
    end

    if (y((osc_ind-1)*4 + 1) >= 0 && y((osc_ind-1)*4 + 1) <= 2*pi)

        oscillator_pos(calc_ind_theta, :) = pospos1;

    elseif (y((osc_ind-1)*4 + 1) > 2*pi && y((osc_ind-1)*4 + 1) <= 4*pi)

        oscillator_pos(calc_ind_theta, :) = pospos2(1,:);

    else

        oscillator_pos(calc_ind_theta, :) = pospos3(2,:);

    end

    % for NC 11 daughter of founding oscillator, calculate d_theta_dt
    osc_sub_ind = osc_sub_ind+1;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    % non_lin_function(a_inhib, theta_inhib, kappa_von_mises, muval_von_mises)
    inhib_output = non_lin_function(y(calc_ind_a), y(calc_ind_theta), kappa_von_mises, muval_von_mises);
    if (y((osc_ind-1)*4 + 2)>4*pi)
        if (inhib_val3(4) > pct_cutoff)
            dydt(calc_ind_theta) = 1 - inhib_output + Fprop * sin(t - y(calc_ind_theta)) .* heaviside((y(calc_ind_theta)-4*pi)*-1);
        else
            dydt(calc_ind_theta) = 1 - inhib_output + Fprop * sin(t - y(calc_ind_theta));
            which_oscillators_soma(calc_ind_theta) = 1;
        end
        which_oscillators_on(calc_ind_theta) = 1;
    end

    if (y((osc_ind-1)*4 + 1) >= 0 && y((osc_ind-1)*4 + 1) <= 2*pi)

        oscillator_pos(calc_ind_theta, :) = pospos1;

    elseif (y((osc_ind-1)*4 + 2) > 2*pi && y((osc_ind-1)*4 + 2) <= 4*pi)

        oscillator_pos(calc_ind_theta, :) = pospos2(2,:);

    else

        oscillator_pos(calc_ind_theta, :) = pospos3(4,:);

    end

end

% which_oscillators_on is of size [num_oscillators_NC9*4,1]
% modify coupling terms to allow for local (and/or global)
check_oscillators = find(which_oscillators_on);
store_positions = oscillator_pos(check_oscillators, :);

[C] = unique(store_positions,'rows');
if (size(C, 1) ~= size(store_positions, 1))
    disp(store_positions);
end
% also cut external bonds
assert(size(C, 1) == size(store_positions, 1))
store_positions = vertcat(store_positions, non_lattice_pos);
store_indices = vertcat(check_oscillators, inf*ones(size(non_lattice_pos,1),1));
assert(size(store_indices, 1) == size(store_positions, 1))
DT = delaunayTriangulation(store_positions(:,1),store_positions(:,2));

% % 2. Extract triangle vertices
% tr = DT.ConnectivityList;
% V1 = DT.Points(tr(:,1), :);
% V2 = DT.Points(tr(:,2), :);
% V3 = DT.Points(tr(:,3), :);
assert(isequal(DT.Points, store_positions))

edges_DT = edges(DT);
edges_length_DT = zeros(size(edges_DT,1),1);
for jjjj = 1:length(edges_length_DT)
    edges_length_DT(jjjj) = norm(DT.Points(edges_DT(jjjj,1), :)-...
        DT.Points(edges_DT(jjjj,2), :));
end

% store d_a_dt for each oscillator
for osc_ind = 1:num_oscillators_NC9
    % NC 9
    inhib_val1 = sorted_cell_output_NC_9_10_11_gp{osc_ind,1};
    % NC 10
    inhib_val2 = sorted_cell_output_NC_9_10_11_gp{osc_ind,2};
    % NC 11
    inhib_val3 = sorted_cell_output_NC_9_10_11_gp{osc_ind,3};

    % for founding oscillator, calculate a
    osc_sub_ind = 1;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    if (y(calc_ind_theta) < 0)
        assert(1==2)
    end
    if (y(calc_ind_theta) >= 0 && y(calc_ind_theta) <= 2*pi)

        dydt(calc_ind_a) = k_time * inhib_val1;

    elseif (y(calc_ind_theta) > 2*pi && y(calc_ind_theta) <= 4*pi)

        dydt(calc_ind_a) = k_time * 2 * inhib_val2(1);

    elseif (y(calc_ind_theta) > 4*pi)

        dydt(calc_ind_a) = k_time * 4 * inhib_val3(1);

    end

    dydt = add_diffusive_coupling(dydt, y, calc_ind_theta, edges_DT, edges_length_DT, ...
        store_indices, which_oscillators_on, which_oscillators_soma, num_oscillators_NC9, Diffuse);

    % second oscillator
    osc_sub_ind = osc_sub_ind + 1;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    dydt = add_diffusive_coupling(dydt, y, calc_ind_theta, edges_DT, edges_length_DT, ...
        store_indices, which_oscillators_on, which_oscillators_soma, num_oscillators_NC9, Diffuse);

    % osc 1
    if (y((osc_ind-1)*4 + 1) >= 0 && y((osc_ind-1)*4 + 1) <= 2*pi)

        dydt(calc_ind_a) = dydt(calc_ind_a-1);

    % osc 2
    elseif (y(calc_ind_theta) <= 4*pi)

        dydt(calc_ind_a) = dydt(calc_ind_a) + k_time * 2 * inhib_val2(2);

    elseif (y(calc_ind_theta) > 4*pi)

        dydt(calc_ind_a) = dydt(calc_ind_a) + k_time * 4 * inhib_val3(3);

    else

        assert(1==2)

    end


    % third oscillator
    osc_sub_ind = osc_sub_ind + 1;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    dydt = add_diffusive_coupling(dydt, y, calc_ind_theta, edges_DT, edges_length_DT, ...
        store_indices, which_oscillators_on, which_oscillators_soma, num_oscillators_NC9, Diffuse);

    % osc 1
    if (y((osc_ind-1)*4 + 1) >= 0 && y((osc_ind-1)*4 + 1) <= 2*pi)

        dydt(calc_ind_a) = dydt(calc_ind_a-2);

        % osc 1
    elseif (y((osc_ind-1)*4 + 1) > 2*pi && y((osc_ind-1)*4 + 1) <= 4*pi)

        dydt(calc_ind_a) = dydt(calc_ind_a-2);

    else

        dydt(calc_ind_a) = dydt(calc_ind_a) + k_time * 4 * inhib_val3(2);

    end


    % fourth oscillator
    osc_sub_ind = osc_sub_ind + 1;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    dydt = add_diffusive_coupling(dydt, y, calc_ind_theta, edges_DT, edges_length_DT, ...
        store_indices, which_oscillators_on, which_oscillators_soma, num_oscillators_NC9, Diffuse);

    % osc 1
    if (y((osc_ind-1)*4 + 1) >= 0 && y((osc_ind-1)*4 + 1) <= 2*pi)

        dydt(calc_ind_a) = dydt(calc_ind_a-3);

    % osc 2
    %elseif (y((osc_ind-1)*4 + 2) > 2*pi && y((osc_ind-1)*4 + 2) <= 4*pi)
    elseif (y((osc_ind-1)*4 + 1) >= 2*pi && y((osc_ind-1)*4 + 2) <= 4*pi)

        dydt(calc_ind_a) = dydt(calc_ind_a-2);

    else

        dydt(calc_ind_a) = dydt(calc_ind_a) + k_time * 4 * inhib_val3(4);

    end

end

%% NO-BACKWARD-CROSSING GUARD
% Once theta passes 2*pi*n it may not return below it: forbid negative
% velocity in a band of width delta_guard just above each multiple.
delta_guard = 0.25;                 % must be < muval_von_mises (see notes)
n_theta     = num_oscillators_NC9*4;
theta       = y(1:n_theta);
n_cycle     = floor(theta/(2*pi) + 1e-9);      % tolerance: 6*pi/(2*pi) may round to 2.999...
phase       = theta - 2*pi*n_cycle;            % in [~0, 2*pi)
guard       = (phase < delta_guard) & (dydt(1:n_theta) < 0);
dydt(guard) = 0;


end

function inhib_output = non_lin_function(a_inhib, theta_inhib, kappa_von_mises, muval_von_mises)

inhib_output = a_inhib;

inhib_output = inhib_output .* ...
    exp(kappa_von_mises * cos(theta_inhib - muval_von_mises)) / ...
    exp(kappa_von_mises);

end

function dydt = add_diffusive_coupling(dydt, y, calc_ind_theta, edges_DT, edges_length_DT, ...
    store_indices, which_oscillators_on, which_oscillators_soma, num_oscillators_NC9, Diffuse)
% Add diffusive coupling of inhibitor a to oscillator calc_ind_theta from its
% Delaunay neighbors. Coupling is active only if both this oscillator and the
% neighbor are on, and each is either somatic or has not yet passed 4*pi.
% Edges touching non-lattice points (store_indices == inf) are ignored.

if (which_oscillators_on(calc_ind_theta))
    if (which_oscillators_soma(calc_ind_theta) || (y(calc_ind_theta) < 4*pi))
        for jjjj = 1:size(edges_DT,1)
            this_1 = edges_DT(jjjj,1);
            this_2 = edges_DT(jjjj,2);

            if (~isinf(store_indices(this_1)) && ~isinf(store_indices(this_2)))
                if (store_indices(this_1) == calc_ind_theta)
                    flag_other = (which_oscillators_soma(store_indices(this_2)) || (y(store_indices(this_2)) < 4*pi));

                    if (which_oscillators_on(store_indices(this_2)) && flag_other)
                        %store_angle = coupling_function(DT, angles, this_1, this_2, store_indices);
                        calc_ind_a_this_1 = num_oscillators_NC9*4 + store_indices(this_1);
                        calc_ind_a_this_2 = num_oscillators_NC9*4 + store_indices(this_2);
                        % prefactor_val = 1 / 2 / voro_area_vect(this_1);
                        % prefactor_val = prefactor_val * sum(store_angle);
                        prefactor_val = 1/(edges_length_DT(jjjj))^2;
                        dydt(calc_ind_a_this_1) = dydt(calc_ind_a_this_1) + ...
                            Diffuse * prefactor_val * (y(calc_ind_a_this_2)-y(calc_ind_a_this_1));
                    end
                end
                if (store_indices(this_2) == calc_ind_theta)
                    flag_other = (which_oscillators_soma(store_indices(this_1)) || (y(store_indices(this_1)) < 4*pi));

                    if (which_oscillators_on(store_indices(this_1)) && flag_other)
                        %store_angle = coupling_function(DT, angles, this_1, this_2, store_indices);
                        calc_ind_a_this_1 = num_oscillators_NC9*4 + store_indices(this_1);
                        calc_ind_a_this_2 = num_oscillators_NC9*4 + store_indices(this_2);
                        % prefactor_val = 1 / 2 / voro_area_vect(this_2);
                        % prefactor_val = prefactor_val * sum(store_angle);
                        prefactor_val = 1/(edges_length_DT(jjjj))^2;
                        dydt(calc_ind_a_this_2) = dydt(calc_ind_a_this_2) + ...
                            Diffuse * prefactor_val * (y(calc_ind_a_this_1)-y(calc_ind_a_this_2));
                    end
                end
            end
        end
    end
end

end

