function [div_counter, store_index, which_above, which_osc_ind_vect, gp_vect, extra_div_counter] = ...
    plot_cycle_traces(cell_output_NC_9_10_11_gp,color_vect,...
    t_val,y_val, plot_at_end, pct_cutoff, assert_check_bool)

%% OUTPUTS:
% div_counter: count number of PGCs that divide once after NC 11
% store_index: store index for dividing PGCs
% which_above: stores which oscillators above threshold at NC 11
% which_osc_ind_vect: stores indices for clones that generate at least one PGC
% gp_vect: distribution of germ plasm per PGC at NC 14
% extra_div_counter: count number of PGCs that divide twice after NC 11

%% INPUTS:
% cell_output_NC_9_10_11_gp: cell array for storing germ plasm amount over lineages
% color_vect: store vector of colors -- for visualization purposes
% t_val: time points in simulation output
% y_val: phases and inhibitor concentrations (see odefun_cell_cycle_oscillotile)
% plot_at_end: should I plot at end? If running many realizations at once,
%   set to false.
% pct_cutoff: if the daughter of a bud has more than this percentage of
%   germ plasm at NC11, cellularize.
% assert_check_bool: if wild-type germ plasm amount, then check that at max
%   one division per PGC; otherwise, count additional divisions

color_vect = vertcat(color_vect, zeros(size(color_vect)));
color_vect = vertcat(color_vect, zeros(size(color_vect)));
color_vect = vertcat(color_vect, zeros(size(color_vect)));
color_vect = vertcat(color_vect, zeros(size(color_vect)));


which_above = [];
store_index = [];
which_osc_ind_vect = [];
div_counter = 0;
extra_div_counter = 0;
gp_vect = [];

num_oscillators_NC9 = size(cell_output_NC_9_10_11_gp, 1);

which_cycle_soma = zeros(num_oscillators_NC9*4,1);

for osc_ind = 1:num_oscillators_NC9
    inhib_val3 = cell_output_NC_9_10_11_gp{osc_ind,3};

    osc_sub_ind = 1;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    if (inhib_val3(1) <= pct_cutoff)
        which_cycle_soma(calc_ind_theta) = 1;
    end

    osc_sub_ind = 2;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    if (inhib_val3(3) <= pct_cutoff)
        which_cycle_soma(calc_ind_theta) = 1;
    end

    osc_sub_ind = 3;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    if (inhib_val3(2) <= pct_cutoff)
        which_cycle_soma(calc_ind_theta) = 1;
    end

    osc_sub_ind = 4;
    calc_ind_theta = (osc_ind-1)*4 + osc_sub_ind;
    if (inhib_val3(4) <= pct_cutoff)
        which_cycle_soma(calc_ind_theta) = 1;
    end
end

for osc_ind = 1:size(cell_output_NC_9_10_11_gp, 1)

    % NC 11
    inhib_val3 = cell_output_NC_9_10_11_gp{osc_ind,3};

    if (~which_cycle_soma((osc_ind-1)*4 + 1) || ~which_cycle_soma((osc_ind-1)*4 + 2) || ...
            ~which_cycle_soma((osc_ind-1)*4 + 3) || ~which_cycle_soma((osc_ind-1)*4 + 4))
        if (plot_at_end)
            figure; hold all;
            title(osc_ind)
        end

        if (~which_cycle_soma((osc_ind-1)*4 + 1))
            if (plot_at_end)
                plot(t_val, y_val(:,(osc_ind-1)*4 + 1),'Color',color_vect(osc_ind,:));
            end
            max_val1 = max(y_val(:,(osc_ind-1)*4 + 1));
            which_above = vertcat(which_above, (osc_ind-1)*4 + 1);
            which_osc_ind_vect = vertcat(which_osc_ind_vect, osc_ind);
            if (max_val1>6*pi && max_val1<8*pi)
                div_counter = div_counter +1;
                store_index = vertcat(store_index, (osc_ind-1)*4 + 1);
                gp_vect = vertcat(gp_vect, inhib_val3(1)/2);
                gp_vect = vertcat(gp_vect, inhib_val3(1)/2);
            elseif (max_val1>8*pi && max_val1<10*pi)
                extra_div_counter = extra_div_counter + 1;
                store_index = vertcat(store_index, (osc_ind-1)*4 + 1);
                gp_vect = vertcat(gp_vect, inhib_val3(1)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(1)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(1)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(1)/4);
            else
                gp_vect = vertcat(gp_vect, inhib_val3(1));
            end
            assert(max_val1<10*pi)
            assert(max_val1>4*pi)
            if (assert_check_bool)
                assert(max_val1<8*pi);
            end
        else
            if (plot_at_end)
                find_t_val = find(t_val<=4*pi);
                hold all; plot(t_val(find_t_val), y_val(find_t_val,(osc_ind-1)*4 + 1),'Color',color_vect(osc_ind,:));
            end
        end
        if (~which_cycle_soma((osc_ind-1)*4 + 2))
            if (plot_at_end)
                find_t_val = find(y_val(:,(osc_ind-1)*4 + 2)>2*pi);
                hold all; plot(t_val(find_t_val), y_val(find_t_val,(osc_ind-1)*4 + 2),'Color',color_vect(osc_ind,:));
            end
            max_val2 = max(y_val(:,(osc_ind-1)*4 + 2));
            which_above = vertcat(which_above, (osc_ind-1)*4 + 2);
            which_osc_ind_vect = vertcat(which_osc_ind_vect, osc_ind);
            if (max_val2>6*pi && max_val2<8*pi)
                div_counter = div_counter +1;
                store_index = vertcat(store_index, (osc_ind-1)*4 + 2);
                gp_vect = vertcat(gp_vect, inhib_val3(3)/2);
                gp_vect = vertcat(gp_vect, inhib_val3(3)/2);
            elseif (max_val2>8*pi && max_val2<10*pi)
                extra_div_counter = extra_div_counter + 1;
                store_index = vertcat(store_index, (osc_ind-1)*4 + 2);
                gp_vect = vertcat(gp_vect, inhib_val3(3)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(3)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(3)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(3)/4);
            else
                gp_vect = vertcat(gp_vect, inhib_val3(3));
            end
            assert(max_val2<10*pi)
            assert(max_val2>4*pi)
            if (assert_check_bool)
                assert(max_val2<8*pi);
            end
        else
            if (plot_at_end)
                find_t_val = find(t_val<=4*pi & y_val(:,(osc_ind-1)*4 + 2)>2*pi);
                hold all; plot(t_val(find_t_val), y_val(find_t_val,(osc_ind-1)*4 + 2),'Color',color_vect(osc_ind,:));
            end
        end
        if (~which_cycle_soma((osc_ind-1)*4 + 3))
            if (plot_at_end)
                find_t_val = find(y_val(:,(osc_ind-1)*4 + 3)>4*pi);
                hold all; plot(t_val(find_t_val), y_val(find_t_val,(osc_ind-1)*4 + 3),'Color',color_vect(osc_ind,:));
            end
            max_val3 = max(y_val(:,(osc_ind-1)*4 + 3));
            which_above = vertcat(which_above, (osc_ind-1)*4 + 3);
            which_osc_ind_vect = vertcat(which_osc_ind_vect, osc_ind);
            if (max_val3>6*pi && max_val3<8*pi)
                div_counter = div_counter +1;
                store_index = vertcat(store_index, (osc_ind-1)*4 + 3);
                gp_vect = vertcat(gp_vect, inhib_val3(2)/2);
                gp_vect = vertcat(gp_vect, inhib_val3(2)/2);
            elseif (max_val3>8*pi && max_val3<10*pi)
                extra_div_counter = extra_div_counter + 1;
                store_index = vertcat(store_index, (osc_ind-1)*4 + 3);
                gp_vect = vertcat(gp_vect, inhib_val3(2)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(2)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(2)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(2)/4);
            else
                gp_vect = vertcat(gp_vect, inhib_val3(2));
            end
            assert(max_val3<10*pi)
            assert(max_val3>4*pi)
            if (assert_check_bool)
                assert(max_val3<8*pi);
            end
        else
            if (plot_at_end)
                find_t_val = find(t_val<=4*pi & y_val(:,(osc_ind-1)*4 + 3)>4*pi);
                hold all; plot(t_val(find_t_val), y_val(find_t_val,(osc_ind-1)*4 + 3),'Color',color_vect(osc_ind,:));
            end
        end
        if (~which_cycle_soma((osc_ind-1)*4 + 4))
            if (plot_at_end)
                find_t_val = find(y_val(:,(osc_ind-1)*4 + 4)>4*pi);
                hold all; plot(t_val(find_t_val), y_val(find_t_val,(osc_ind-1)*4 + 4),'Color',color_vect(osc_ind,:));
            end
            max_val4 = max(y_val(:,(osc_ind-1)*4 + 4));
            which_above = vertcat(which_above, (osc_ind-1)*4 + 4);
            which_osc_ind_vect = vertcat(which_osc_ind_vect, osc_ind);
            if (max_val4>6*pi && max_val4<8*pi)
                div_counter = div_counter +1;
                store_index = vertcat(store_index, (osc_ind-1)*4 + 4);
                gp_vect = vertcat(gp_vect, inhib_val3(4)/2);
                gp_vect = vertcat(gp_vect, inhib_val3(4)/2);
            elseif (max_val4>8*pi && max_val4<10*pi)
                extra_div_counter = extra_div_counter + 1;
                store_index = vertcat(store_index, (osc_ind-1)*4 + 4);
                gp_vect = vertcat(gp_vect, inhib_val3(4)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(4)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(4)/4);
                gp_vect = vertcat(gp_vect, inhib_val3(4)/4);
            else
                gp_vect = vertcat(gp_vect, inhib_val3(4));
            end
            assert(max_val4<10*pi)
            assert(max_val4>4*pi)
            if (assert_check_bool)
                assert(max_val4<8*pi);
            end
        else
            if (plot_at_end)
                find_t_val = find(t_val<=4*pi & y_val(:,(osc_ind-1)*4 + 4)>4*pi);
                hold all; plot(t_val(find_t_val), y_val(find_t_val,(osc_ind-1)*4 + 4),'Color',color_vect(osc_ind,:));
            end
        end
        if (plot_at_end)
            hold all; plot(t_val, 6*pi*ones(size(t_val)),'--k')
            hold all; plot(t_val, 8*pi*ones(size(t_val)),'--k')
            ylim([0,10*pi]);
            yticks(0:(2*pi):(10*pi))
            yticklabels({'0','2 \pi','4 \pi', '6 \pi','8 \pi','10 \pi'})
            xlim([0,max(t_val)])
            xticks(0:(2*pi):max(t_val))
            xticklabels({'0','2 \pi','4 \pi', '6 \pi','8 \pi','10 \pi','12 \pi',...
                '14 \pi','16 \pi', '18 \pi','20 \pi'})
            ylabel('\Theta_i')
            xlabel('time')
        end
    end
end

which_osc_ind_vect = unique(which_osc_ind_vect);

end