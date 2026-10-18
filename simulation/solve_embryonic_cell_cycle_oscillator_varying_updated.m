function [t_val,y_val,which_cycle_soma,sorted_cell_output_NC_9_10_11_gp,sorted_cell_output_NC_9_10_11] = ...
    solve_embryonic_cell_cycle_oscillator_varying_updated(cell_output_NC_9_10_11, ...
    cell_output_NC_9_10_11_gp, kappa_von_mises, muval_von_mises, ...
    k_time, theta_init, total_time, pct_cutoff, Diffuse, non_lattice_pos, tp_step)
% NOTE: check Phase reversal

%% outputs:
% t_val
% y_val
% which_cycle_soma
% sorted_cell_output_NC_9_10_11_gp
% sorted_cell_output_NC_9_10_11

%% inputs:
% cell_output_NC_9_10_11
% cell_output_NC_9_10_11_gp
% kappa_von_mises
% muval_von_mises
% k_time
% theta_init
% total_time
% pct_cutoff
% Diffuse
% non_lattice_pos
% tp_step


%% sort vectors

% count number of founders
num_oscillators_NC9 = size(cell_output_NC_9_10_11_gp, 1);

% sort by amount of germ plasm (at each depth in lineage)
sorted_cell_output_NC_9_10_11_gp = cell_output_NC_9_10_11_gp;
sorted_cell_output_NC_9_10_11 = cell_output_NC_9_10_11;

% for all oscillators
for osc_ind = 1:num_oscillators_NC9
    % store vectors
    this_vect = sorted_cell_output_NC_9_10_11_gp{osc_ind, 2};
    this_vect_pos = sorted_cell_output_NC_9_10_11{osc_ind, 2};

    % store values at NC 11
    this_vect_next = sorted_cell_output_NC_9_10_11_gp{osc_ind, 3};
    this_vect_next_pos = sorted_cell_output_NC_9_10_11{osc_ind, 3};

    % sort second entry
    [~,sorted_NC10_I] = sort(this_vect);
    
    assert(isequal(size(sorted_NC10_I), size([1,2]')))
    % sort third entry
    if (~isequal(sorted_NC10_I, [1,2]'))
        % sort secondary element
        sorted_cell_output_NC_9_10_11_gp{osc_ind, 2} = this_vect(sorted_NC10_I);
        sorted_cell_output_NC_9_10_11{osc_ind, 2} = this_vect_pos(sorted_NC10_I, :);

        % sort tertiary element
        [ssss, IIII] = sort(this_vect_next(3:4));
        sorted_cell_output_NC_9_10_11_gp{osc_ind, 3}(1:2) = ssss;
        this_vect_next_pos_34 = this_vect_next_pos(3:4,:);
        sorted_cell_output_NC_9_10_11{osc_ind, 3}(1:2,:) = this_vect_next_pos_34(IIII,:);
        [ssss, IIII] = sort(this_vect_next(1:2));
        sorted_cell_output_NC_9_10_11_gp{osc_ind, 3}(3:4) = ssss;
        this_vect_next_pos_12 = this_vect_next_pos(1:2,:);
        sorted_cell_output_NC_9_10_11{osc_ind, 3}(3:4,:) = this_vect_next_pos_12(IIII,:);
    else
        % sort tertiary element
        [ssss, IIII] = sort(this_vect_next(1:2));
        sorted_cell_output_NC_9_10_11_gp{osc_ind, 3}(1:2) = ssss;
        this_vect_next_pos_12 = this_vect_next_pos(1:2,:);
        sorted_cell_output_NC_9_10_11{osc_ind, 3}(1:2,:) = this_vect_next_pos_12(IIII,:);
        [ssss, IIII] = sort(this_vect_next(3:4));
        sorted_cell_output_NC_9_10_11_gp{osc_ind, 3}(3:4) = ssss;
        this_vect_next_pos_34 = this_vect_next_pos(3:4,:);
        sorted_cell_output_NC_9_10_11{osc_ind, 3}(3:4,:) = this_vect_next_pos_34(IIII,:);
    end
end
% store copy for sorting again by first entry
copy_sorted_cell_output_NC_9_10_11_gp = sorted_cell_output_NC_9_10_11_gp;
copy_sorted_cell_output_NC_9_10_11 = sorted_cell_output_NC_9_10_11;
first_vect = [];
for osc_ind = 1:num_oscillators_NC9
    first_vect = vertcat(first_vect, ...
        sorted_cell_output_NC_9_10_11_gp{osc_ind,1});
end
% sort by first element
[~,sort_first] = sort(first_vect);
% for all oscillators
for osc_ind = 1:num_oscillators_NC9
    sorted_cell_output_NC_9_10_11_gp{osc_ind,1} = copy_sorted_cell_output_NC_9_10_11_gp{sort_first(osc_ind),1};
    sorted_cell_output_NC_9_10_11_gp{osc_ind,2} = copy_sorted_cell_output_NC_9_10_11_gp{sort_first(osc_ind),2};
    sorted_cell_output_NC_9_10_11_gp{osc_ind,3} = copy_sorted_cell_output_NC_9_10_11_gp{sort_first(osc_ind),3};

    sorted_cell_output_NC_9_10_11{osc_ind,1} = copy_sorted_cell_output_NC_9_10_11{sort_first(osc_ind),1};
    sorted_cell_output_NC_9_10_11{osc_ind,2} = copy_sorted_cell_output_NC_9_10_11{sort_first(osc_ind),2};
    sorted_cell_output_NC_9_10_11{osc_ind,3} = copy_sorted_cell_output_NC_9_10_11{sort_first(osc_ind),3};

    assert(abs(sorted_cell_output_NC_9_10_11_gp{osc_ind,1}-sum(sorted_cell_output_NC_9_10_11_gp{osc_ind,2}))<10^-10)
    this_vect = sorted_cell_output_NC_9_10_11_gp{osc_ind, 2};
    this_vect_next = sorted_cell_output_NC_9_10_11_gp{osc_ind, 3};
    tolc = @(q) max(2e-8, 2e-6*abs(q));
    assert(abs(this_vect(1)-sum(this_vect_next(1:2))) < tolc(this_vect(1)))
    assert(abs(this_vect(2)-sum(this_vect_next(3:4))) < tolc(this_vect(2)))
    this_v = sorted_cell_output_NC_9_10_11{osc_ind, 1};
    this_vect = sorted_cell_output_NC_9_10_11{osc_ind, 2};
    this_vect_next = sorted_cell_output_NC_9_10_11{osc_ind, 3};
    assert(abs(this_v(1) - mean(this_vect(:,1)))<10^-10)
    assert(abs(this_v(2) - mean(this_vect(:,2)))<10^-10)
    assert(abs(this_vect(1,1) - mean(this_vect_next(1:2,1)))<10^-10)
    assert(abs(this_vect(1,2) - mean(this_vect_next(1:2,2)))<10^-10)
    assert(abs(this_vect(2,1) - mean(this_vect_next(3:4,1)))<10^-10)
    assert(abs(this_vect(2,2) - mean(this_vect_next(3:4,2)))<10^-10)
end

assert(theta_init<2*pi && theta_init>=0)
% solve ode fun oscillator
num_oscillators_NC9 = size(cell_output_NC_9_10_11_gp, 1);
% check how system state is passed
system_state = zeros(1, num_oscillators_NC9*4*2);
which_cycle_soma = zeros(num_oscillators_NC9*4,1);

for osc_ind = 1:num_oscillators_NC9
    inhib_val3 = sorted_cell_output_NC_9_10_11_gp{osc_ind,3};

    osc_sub_ind = 1;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    system_state(calc_ind_a) = 0;
    system_state(calc_ind_theta) = theta_init;
    if (inhib_val3(1) <= pct_cutoff)
        which_cycle_soma(calc_ind_theta) = 1;
    end

    osc_sub_ind = 2;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    system_state(calc_ind_a) = 0;
    system_state(calc_ind_theta) = 2*pi;
    if (inhib_val3(3) <= pct_cutoff)
        which_cycle_soma(calc_ind_theta) = 1;
    end

    osc_sub_ind = 3;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    system_state(calc_ind_a) = 0;
    system_state(calc_ind_theta) = 4*pi;
    if (inhib_val3(2) <= pct_cutoff)
        which_cycle_soma(calc_ind_theta) = 1;
    end

    osc_sub_ind = 4;
    calc_ind_a = num_oscillators_NC9*4 + (osc_ind-1)*4 + osc_sub_ind;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    system_state(calc_ind_a) = 0;
    system_state(calc_ind_theta) = 4*pi;
    if (inhib_val3(4) <= pct_cutoff)
        which_cycle_soma(calc_ind_theta) = 1;
    end
end

options_chosen = ...
    odeset('MaxStep',tp_step,'RelTol', 1e-6, 'Stats','on');
[t_val,y_val] = ode45(@(t,y)odefun_cell_cycle_oscillotile(t, y, ...
    sorted_cell_output_NC_9_10_11_gp, kappa_von_mises, muval_von_mises, ...
    k_time, pct_cutoff, Diffuse, sorted_cell_output_NC_9_10_11, non_lattice_pos), ...
    theta_init:(tp_step):(total_time+theta_init), system_state, options_chosen );

n_theta = size(sorted_cell_output_NC_9_10_11_gp, 1) * 4;   % theta columns come first in y_val
tol = 1e-3;                                                % allow for solver round-off
all_ok  = true;

for k = 1:n_theta
    theta = y_val(:, k);
    thr   = floor(max(theta)/2/pi)*2*pi;

    % first time point at which theta reaches the threshold
    first_cross = find(theta >= thr, 1, 'first');

    % any later point that falls back below it?
    bad = find(theta(first_cross:end) < thr - tol, 1, 'first');
    if ~isempty(bad)
        bad = first_cross - 1 + bad;
        fprintf('Oscillator %d drops below %d*2pi at t = %.4f (theta = %.6f)\n', ...
            k, round(thr/(2*pi)), t_val(bad), theta(bad));
        all_ok = false;
    end
end

if all_ok
    disp('No oscillator crosses back below its highest 2*pi multiple.');
end

end