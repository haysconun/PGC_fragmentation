function [ lattice_pos, pairs, percent_germ_plasm, percent_germ_plasm_half, ...
    summed_density, real_bud_counter, above_thresh_pairs, above_thresh_pairs_dist_cutoff, ...
    cell_output_NC_9_10_11, cell_output_NC_9_10_11_gp, mean_points ] = ...
    generate_hexagon_tesselation(A, rnd_mag, rnd_angle, a_val, num_points_pm, dist_cut, ...
    pct_cutoff, rand_seed, NC11_flag, should_i_plot, plot_NC9, plot_NC10, plot_NC11)
%% OUTPUTS: 

% lattice_pos: regular array of nuclear positions at NC 10
% pairs: paired NC 10 sites -- assumption: NC 9 sites between pairs NC 10 sites
% percent_germ_plasm: percent germ plasm associated with NC 10 sites (same ordering as lattice_pos)
% percent_germ_plasm_half: percent germ plasm associated with NC 11 sites (same ordering as lattice_pos, but with 2)
% summed_density: percent germ plasm associated with NC 9 nucleus (by summing two NC 10 daughters)
% real_bud_counter: count number of NC 11 daughters that cellularize
% above_thresh_pairs: count which pairs of NC 10 (aka NC 9 site) has at least one daughter at NC 11 that is a cell
% above_thresh_pairs_dist_cutoff: which NC 9 nuclei are sufficiently close to the center (based on dist_cut)
% cell_output_NC_9_10_11: positions of nuclei at 9, 10, 11 for each clone
% cell_output_NC_9_10_11_gp: germ plasm % associated with nuclei at 9, 10, 11 for each clone
% mean_points: mean position of paired NC 10 sites

%% INPUTS: 

% sigma = 1 (see below)
% A = 1 / 2 / pi / sigma^2 % Gaussian amplitude
% rnd_mag = 1/4 of sigma
% rnd_angle = pi/2
% a_val = 1.1
% num_points_pm = 4
% dist_cut = 14
% pct_cutoff = 0.015
% rand_seed = 8696
% NC11_flag = 2

% A: amplitude of gaussian
% rnd_mag: rand displacement magnitude
% rnd_angle: rand displacement angle
% a_val: inter-nuclear spacing at NC 10 (relative to germ plasm gradient)
% num_points_pm: build lattice -num_points_pm:num_points_pm along lattice vectors
% dist_cut: distance cutoff in units of sigma for lattice points (UPDATE UPDATE UPDATE)
% pct_cutoff: pct_cutoff for what is a cell
% rand_seed: seed for random number generation

% % 3 OPTIONS FOR NC 11 divisions:
% NC11_flag == 1 -- divide each bud based on the vector connecting the two related buds
% NC11_flag == 2 -- divide each bud with random orientation but with same random orientation of sister bud
% NC11_flag == 3 -- divide each bud with independent random orientation
% NC11_flag == 4 -- equal partitioning of germ plasm starting at division of buds

% should_i_plot: should I generate plot

% plot_NC9, plot_NC10, plot_NC11: should I plot different germ plasm distributions?

% NOTE (for below):
% dist_cut is actually a number of nuclei to keep
% integrate_gaussian_polygon does not generalize well to really large polygons (well outside of parameter range of interest)

%%

% seed random number generator
rng(rand_seed);

% generate centers for hexagons
sigma = 1; % always fix to be 1 -- change a_val relative to this.
mu = [0,0]; % center of Gaussian -- lattice displaced relative to this
lattice_pos = []; % store nuclear positions in lattice
percent_germ_plasm = []; % store percent germ plasm associated with each lattice point
R = a_val / sqrt(3); % store circumradius of each hexagon (voronoi cell)
% absolute and relative tolerances
ABSTOL = 1e-10;
RELTOL = 1e-8;
% don't rotate lattice of hexagons!
phi = 0;

% build lattice about center
for mmm = -num_points_pm:num_points_pm
    for nnn = -num_points_pm:num_points_pm
        c = mmm*a_val*[0,1]+nnn*a_val*[sind(120), cosd(120)];
        c(1) = c(1) + rnd_mag * cos(rnd_angle);
        c(2) = c(2) + rnd_mag * sin(rnd_angle);
        % store position in lattice.
        lattice_pos = vertcat(lattice_pos, c);
        % store hexagon vertices
        [vx, vy] = hexagon_vertices(R, c, phi);
        verts     = [vx(:), vy(:)];
        I_quad = integrate_gaussian_polygon(A, mu, sigma, verts, ABSTOL, RELTOL);
        % [I_mc, ~] = MC_integrate(verts, sigma, mu, A);
        % assert(abs(I_mc - I_quad)<=tol_val);

        % store percent germ plasm.
        percent_germ_plasm = vertcat(percent_germ_plasm, I_quad);
    end
end

% plot gradient with no mask
plot_r  = max(6*sigma, R + 0.5);
xv      = linspace(- plot_r, plot_r, 400);
yv      = linspace(- plot_r, plot_r, 400);
[xx, yy] = meshgrid(xv, yv);
zz      = A * exp(-((xx).^2 + (yy).^2) / (2*sigma^2));
zz_masked = zz;
% compute voronoi cells
[vMGFPoints, vMGFcells] = voronoin(lattice_pos);

% should I plot
if (should_i_plot)
    % load figure
    figure('Color','w','Position',[100 100 1100 450]);

    % --- Left: Gaussian heatmap inside hexagon ---
    subplot(1,2,1);
    imagesc(xv, yv, zz_masked);
    axis xy; axis equal tight;
    n = 256;
    cmap = interp1([0; 1], [0 0 0; 1 0 0], linspace(0,1,n));
    colormap(cmap);
    colorbar;
    hold on;

    for vcell = 1:length(lattice_pos)
        patch(vMGFPoints(vMGFcells{vcell},1),...
            vMGFPoints(vMGFcells{vcell},2), 'w','FaceAlpha',0.2);
    end
end

% Do random matching of NC 10 sites -- this determines NC 9 positions
pairs = random_hex_matching(lattice_pos, a_val);
% calculate NC 9 positions
mean_points = zeros(size(pairs,1), 2);
% store random div angle
random_div_angle = zeros(size(lattice_pos,1), 1);
% store summed density between paired NC 10 sites.
summed_density = zeros(size(pairs,1), 1);
% loop through pairs.
for iii = 1:size(pairs,1)
    % plot bond between paired lattice sites (with cyan dot on line)
    if (should_i_plot)
        hold all; plot(lattice_pos(pairs(iii,:),1), lattice_pos(pairs(iii,:),2), 'k','LineWidth',3)
        hold all; scatter(mean(lattice_pos(pairs(iii,:),1)), mean(lattice_pos(pairs(iii,:),2)),'cyan','filled')
    end

    % calculate mean point for each bond
    mean_points(iii,1) = mean(lattice_pos(pairs(iii,:),1));
    mean_points(iii,2) = mean(lattice_pos(pairs(iii,:),2));

    % based on NC11_flag, determine bud division
    % flag == 1 -- divide each bud based on the vector connecting the two related buds
    % flag == 2 -- divide each bud with random orientation but with same random orientation of sister bud
    % flag == 3 -- divide each bud with independent random orientation
    % flag == 4 -- equal partitioning of germ plasm starting at division of buds
    if (NC11_flag == 1)
        random_div_angle(pairs(iii,1)) = atan2(lattice_pos(pairs(iii,1),2) - lattice_pos(pairs(iii,2),2), lattice_pos(pairs(iii,1),1) - lattice_pos(pairs(iii,2),1));
        random_div_angle(pairs(iii,2)) = atan2(lattice_pos(pairs(iii,1),2) - lattice_pos(pairs(iii,2),2), lattice_pos(pairs(iii,1),1) - lattice_pos(pairs(iii,2),1));
    elseif (NC11_flag == 2)
        rand_val = rand * 2 * pi;
        random_div_angle(pairs(iii,1)) = rand_val;
        random_div_angle(pairs(iii,2)) = rand_val;
    elseif (NC11_flag == 3)
        random_div_angle(pairs(iii,1)) = rand * 2 * pi;
        random_div_angle(pairs(iii,2)) = rand * 2 * pi;
    elseif (NC11_flag == 4)
        rand_val = rand * 2 * pi;
        random_div_angle(pairs(iii,1)) = rand_val;
        random_div_angle(pairs(iii,2)) = rand_val;
    else
        error('NC11_flag not specified correctly')
    end

    % sum plasm over two sister buds at NC 10
    summed_density(iii) = percent_germ_plasm(pairs(iii,1)) + percent_germ_plasm(pairs(iii,2));
end

% count number of bud daughters that cellularize.
real_bud_counter = 0;
% store germ plasm in each half of dividing bud
percent_germ_plasm_half = zeros(length(lattice_pos), 2);
% which pairs of NC 10 (aka NC 9 site) has at least one daughter at NC 11 that is a cell
above_thresh_pairs = zeros(size(pairs,1), 1);
above_thresh_pairs_dist_cutoff = zeros(size(pairs,1), 1);

% loop through lattice sites.
for i = 1:length(lattice_pos)

    % calculate vertices of hexagon
    [vx, vy] = hexagon_vertices(R, lattice_pos(i,:), 0);
    % calculate polyshape
    poly1 = polyshape(vx, vy);
    % calculate line segment cutting hexagon
    lineseg = [lattice_pos(i,1)+2*a_val*cos(random_div_angle(i)) lattice_pos(i,2)+2*a_val*sin(random_div_angle(i)); lattice_pos(i,1)-2*a_val*cos(random_div_angle(i)) lattice_pos(i,2)-2*a_val*sin(random_div_angle(i))];
    % calculate intersect.
    [in,~] = intersect(poly1,lineseg);
    % plot line
    if (should_i_plot)
        hold on
        plot(in(:,1),in(:,2),'b');
    end

    %% --- Parameters ---
    % calculate vertices of hexagon
    [vx, vy] = hexagon_vertices(R, lattice_pos(i,:), 0);
    verts     = [vx(:), vy(:)];

    % cut polygon in half
    [half1, half2] = clip_polygon_by_line(verts, lattice_pos(i,:), random_div_angle(i));

    % integrate gaussian over polygon
    I1     = integrate_gaussian_polygon(A, mu, sigma, half1, ABSTOL, RELTOL);
    % integrate gaussian over polygon
    I2     = integrate_gaussian_polygon(A, mu, sigma, half2, ABSTOL, RELTOL);
    % [I_mc, ~] = MC_integrate(half1, sigma, mu, A);
    % assert(abs(I_mc - I1)<=tol_val);
    % calculate polyshape (half)
    if (NC11_flag == 4)
        mean_I1_I2 = mean([I1, I2]);
        I1 = mean_I1_I2;
        I2 = mean_I1_I2;
        assert(abs(2*mean_I1_I2 - percent_germ_plasm(i)) < ...
            max(200*ABSTOL, 200*RELTOL*abs(percent_germ_plasm(i))));
        %assert(abs(2*mean_I1_I2 - percent_germ_plasm(i)) < 20*ABSTOL);
    end
    pshap1 = polyshape(half1(:,1), half1(:,2));
    % if amount of germ plasm is above cutoff
    if (I1 > pct_cutoff)
        % plot polyshape
        if (should_i_plot)
            hold all; plot(pshap1,'FaceAlpha',0.2,'EdgeColor','w');
        end
        % store that this is above thresh
        above_thresh_pairs(find(pairs(:,1) == i | pairs(:,2) == i)) = 1;
        % increment counter
        real_bud_counter = real_bud_counter + 1;
    end
    % [I_mc, ~] = MC_integrate(half2, sigma, mu, A);
    % assert(abs(I_mc - I2)<=tol_val);
    % calculate polyshape (half)
    pshap2 = polyshape(half2(:,1), half2(:,2));
    if (I2 > pct_cutoff)
        % plot polyshape
        if (should_i_plot)
            hold all; plot(pshap2,'FaceAlpha',0.2,'EdgeColor','w');
        end
        % store that this is above thresh
        above_thresh_pairs(find(pairs(:,1) == i | pairs(:,2) == i)) = 1;
        % increment counter
        real_bud_counter = real_bud_counter + 1;
    end
    % store germ plasm in each half
    percent_germ_plasm_half(i,1) = I1;
    percent_germ_plasm_half(i,2) = I2;


end

% distance from center of gaussian... to NC 9 nucleus
mean_dist_check = zeros(size(mean_points, 1), 1);
for iii = 1:size(pairs,1)
    mean_dist_check(iii) = norm(mean_points(iii,:));
end

[~, sort_III] = sort(mean_dist_check);

% dist_cut is actually a number to keep!!
% store which to keep based on distance of NC 9 nucleus from center.
for iii = 1:min(dist_cut, numel(sort_III))
    above_thresh_pairs_dist_cutoff(sort_III(iii)) = 1;
end

% store the positions of nuclei at 9, 10, 11 for each clone.
cell_output_NC_9_10_11 = cell(length(find(above_thresh_pairs & above_thresh_pairs_dist_cutoff)), 3);
% store the germ plasm associated with nuclei at 9, 10, 11 for each clone.
cell_output_NC_9_10_11_gp = cell(length(find(above_thresh_pairs & above_thresh_pairs_dist_cutoff)), 3);

% counter for which NC9 clones have above_thresh_pairs(i) == true
counter_val_val = 0;
% loop over NC 9 clones.
for i = 1:size(pairs,1)

    % store centroid
    store_centroid_val = [];
    % store germ plasm
    store_gp = [];

    % if at least one NC 11 offspring of this NC 9 clone cellularizes
    if (above_thresh_pairs(i) && above_thresh_pairs_dist_cutoff(i))
        % increment counter
        counter_val_val = counter_val_val + 1;
        % loop over each in pair (NC 10) 
        for j = 1:size(pairs,2)

            % store position at NC 9
            cell_output_NC_9_10_11{counter_val_val,1} = mean_points(i,:);
            % store summed gp at NC 9
            cell_output_NC_9_10_11_gp{counter_val_val,1} = summed_density(i);
            % store positions at NC 10
            cell_output_NC_9_10_11{counter_val_val,2} = lattice_pos(pairs(i,:),:);
            % store percent germ plasm at NC 10
            cell_output_NC_9_10_11_gp{counter_val_val,2} = percent_germ_plasm(pairs(i,:));
            % calculate hexagon vertices
            [vx, vy] = hexagon_vertices(R, lattice_pos(pairs(i,j),:), 0);
            % store verts
            verts = [vx(:), vy(:)];
            % cut hexagon in half
            [half1, half2] = clip_polygon_by_line(verts, ...
                lattice_pos(pairs(i,j),:), random_div_angle(pairs(i,j)));

            I1 = percent_germ_plasm_half(pairs(i,j), 1);
            I2 = percent_germ_plasm_half(pairs(i,j), 2);
            pshap1 = polyshape(half1(:,1), half1(:,2));
            [xvv,yvv] = centroid(pshap1);
            store_centroid_val = vertcat(store_centroid_val, [xvv,yvv]);
            store_gp = vertcat(store_gp, I1);
            pshap2 = polyshape(half2(:,1), half2(:,2));
            [xvv,yvv] = centroid(pshap2);
            store_centroid_val = vertcat(store_centroid_val, [xvv,yvv]);
            store_gp = vertcat(store_gp, I2);
            % store positions and percent gp at NC 11.
            cell_output_NC_9_10_11{counter_val_val, 3} = store_centroid_val;
            cell_output_NC_9_10_11_gp{counter_val_val,3} = store_gp;          
            % check that each half hexagon has the same area.
            assert(abs(area(pshap1) - area(pshap2))<=10^-10)
        end
    end
end

if (should_i_plot)
    subplot(1,2,2);

    % plot gp
    imagesc(xv, yv, zz_masked);
    axis xy; axis equal tight; n = 256;
    cmap = interp1([0; 1], [0 0 0; 1 0 0], linspace(0,1,n));
    colormap(cmap);
    colorbar;
    hold on;
end


% check for conservation of density
[vMGFPoints, vMGFcells] = voronoin(mean_points);

if (should_i_plot)
    % loop through cell output for positions -- plot lineages on gradient
    for vcell = 1:size(cell_output_NC_9_10_11, 1)
        this_point = cell_output_NC_9_10_11{vcell,1};
        other_points = cell_output_NC_9_10_11{vcell,2};
        other_other_points = cell_output_NC_9_10_11{vcell,3};
        hold all; scatter(this_point(1), this_point(2), 'c');
        hold all; quiver(this_point(1), this_point(2), other_points(1,1) - this_point(1), other_points(1,2) - this_point(2), 0,'w');
        hold all; quiver(this_point(1), this_point(2), other_points(2,1) - this_point(1), other_points(2,2) - this_point(2), 0,'w');

        hold all; quiver(other_points(1,1), other_points(1,2), other_other_points(1,1) - other_points(1,1), other_other_points(1,2) - other_points(1,2), 0,'w');
        hold all; quiver(other_points(1,1), other_points(1,2), other_other_points(2,1) - other_points(1,1), other_other_points(2,2) - other_points(1,2), 0,'w');

        hold all; quiver(other_points(2,1), other_points(2,2), other_other_points(3,1) - other_points(2,1), other_other_points(3,2) - other_points(2,2), 0,'w');
        hold all; quiver(other_points(2,1), other_points(2,2), other_other_points(4,1) - other_points(2,1), other_other_points(4,2) - other_points(2,2), 0,'w');
    end
end

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

color_vect = [gray_color;brown_color;blue_color;yellow_color;green_color;orange_color];
colors_extended = distinguishable_colors(size(cell_output_NC_9_10_11,1),color_vect);
color_vect = vertcat(color_vect, colors_extended);
if (plot_NC9)
    figure;
    imagesc(xv, yv, zz_masked);
    axis xy; axis equal tight; n = 256;
    cmap = interp1([0; 1], [0 0 0; 1 0 0], linspace(0,1,n));
    colormap(cmap);
    hold on;
    xticks([]) 
    yticks([]) 

    for vcell = 1:length(mean_points)

        if (above_thresh_pairs(vcell) && above_thresh_pairs_dist_cutoff(vcell))
            
            hold all; patch(vMGFPoints(vMGFcells{vcell},1),vMGFPoints(vMGFcells{vcell},2), 'w','FaceAlpha',0.2); %0.3

        end

    end
    gp_calc = [];
    for vcell = 1:size(cell_output_NC_9_10_11,1)
        this_point = cell_output_NC_9_10_11{vcell,1};
        gp_calc = vertcat(gp_calc, cell_output_NC_9_10_11_gp{vcell,1});
        hold all; scatter(this_point(1), this_point(2),200,color_vect(vcell,:), 'filled');
    end

    storestore_pos = [];
    storestore_pos = vertcat(storestore_pos, mean_points(find(above_thresh_pairs & above_thresh_pairs_dist_cutoff),:));

    gp_calc = gp_calc / sum(gp_calc);

    xxx = 1:length(gp_calc);
    [~,IIIIII] = sort(gp_calc);
    datadata = gp_calc(IIIIII);
    datadata = datadata *100;

    figure; hold all;
    for i = 1:length(gp_calc)
        h1111 = bar(xxx(i),datadata(i));
        h1111.FaceColor = color_vect(IIIIII(i),:);
        h1111.EdgeColor = color_vect(IIIIII(i),:);
    end
    ylim([0,40])
    xticklabels([])
    xlim([0.5,max(xxx)+0.5])
    set(gca, 'FontSize', 25.0);


end

if (plot_NC10)
    figure;
    imagesc(xv, yv, zz_masked);
    axis xy; axis equal tight; n = 256;
    cmap = interp1([0; 1], [0 0 0; 1 0 0], linspace(0,1,n));
    colormap(cmap);
    hold on;
    xticks([])
    yticks([])

    % compute voronoi cells
    [vMGFPoints, vMGFcells] = voronoin(lattice_pos);

    for vcell = 1:length(lattice_pos)
        if (percent_germ_plasm_half(vcell,1)>pct_cutoff||percent_germ_plasm_half(vcell,2)>pct_cutoff)
            patch(vMGFPoints(vMGFcells{vcell},1),vMGFPoints(vMGFcells{vcell},2), 'w','FaceAlpha',0.2,'EdgeColor','w');
        end
    end

    % store positions within both specified bounds.
    storestore_pos = [];
    storestore_pos = vertcat(storestore_pos, lattice_pos(pairs(find(above_thresh_pairs & above_thresh_pairs_dist_cutoff),1),:));
    storestore_pos = vertcat(storestore_pos, lattice_pos(pairs(find(above_thresh_pairs & above_thresh_pairs_dist_cutoff),2),:));

    gp_calc = [];
    color_vect_vect = [];
    for vcell = 1:size(cell_output_NC_9_10_11,1)
        this_point = cell_output_NC_9_10_11{vcell,2};
        hold all; scatter(this_point(:,1), this_point(:,2),100,color_vect(vcell,:), 'filled');
        %hold all; quiver(this_point(1,1), this_point(1,2), this_point(2,1) - this_point(1,1), this_point(2,2) - this_point(1,2), 0, 'Color', color_vect(vcell,:));
  
        this_other_other_vect = cell_output_NC_9_10_11_gp{vcell,3};
        if (this_other_other_vect(1)>pct_cutoff || this_other_other_vect(2)>pct_cutoff)
            gp_calc = vertcat(gp_calc, sum(this_other_other_vect(1:2)));
            color_vect_vect = vertcat(color_vect_vect, color_vect(vcell,:));
        end
        if (this_other_other_vect(3)>pct_cutoff || this_other_other_vect(4)>pct_cutoff)
            gp_calc = vertcat(gp_calc, sum(this_other_other_vect(3:4)));
            color_vect_vect = vertcat(color_vect_vect, color_vect(vcell,:));
        end
    end

    disp(length(gp_calc))
    gp_calc = gp_calc / sum(gp_calc);
    assert(length(gp_calc) == length(color_vect_vect))

    xxx = 1:length(gp_calc);
    [~,IIIIII] = sort(gp_calc);
    datadata = gp_calc(IIIIII);
    datadata = datadata *100;

    figure; hold all;
    for i = 1:length(gp_calc)
        h1111 = bar(xxx(i),datadata(i));
        h1111.FaceColor = color_vect_vect(IIIIII(i),:);
        h1111.EdgeColor = color_vect_vect(IIIIII(i),:);
    end
    ylim([0,20])
    xticklabels([])
    xlim([0.5,max(xxx)+0.5])
    set(gca, 'FontSize', 25.0);
end

if (plot_NC11)

    figure;
    imagesc(xv, yv, zz_masked);
    axis xy; axis equal tight; n = 256;
    cmap = interp1([0; 1], [0 0 0; 1 0 0], linspace(0,1,n));
    colormap(cmap);
    hold on;
    xticks([])
    yticks([])

    for i = 1:size(lattice_pos, 1)
        
        [vx, vy] = hexagon_vertices(R, lattice_pos(i,:), 0);
        verts     = [vx(:), vy(:)];

        % cut polygon in half
        [half1, half2] = clip_polygon_by_line(verts, lattice_pos(i,:), random_div_angle(i));

        % calculate polyshape (half)
        pshap1 = polyshape(half1(:,1), half1(:,2));
        pshap2 = polyshape(half2(:,1), half2(:,2));

        if (percent_germ_plasm_half(i,1)>pct_cutoff)
            plot(pshap1,'FaceAlpha',0.2,'EdgeColor','w','FaceColor','w');
        end
        if (percent_germ_plasm_half(i,2)>pct_cutoff)
            plot(pshap2,'FaceAlpha',0.2,'EdgeColor','w','FaceColor','w');
        end

    end


    gp_calc = [];
    color_vect_vect = [];
    for vcell = 1:size(cell_output_NC_9_10_11, 1)
        this_point = cell_output_NC_9_10_11{vcell,3};
        hold all; scatter(this_point(:,1), this_point(:,2),100,color_vect(vcell,:), 'filled');
        hold all; quiver(this_point(1,1), this_point(1,2), this_point(2,1) - this_point(1,1), this_point(2,2) - this_point(1,2), 0, 'Color', color_vect(vcell,:));
        hold all; quiver(this_point(3,1), this_point(3,2), this_point(4,1) - this_point(3,1), this_point(4,2) - this_point(3,2), 0, 'Color', color_vect(vcell,:));
        this_other_other_vect = cell_output_NC_9_10_11_gp{vcell,3};
        if (this_other_other_vect(1)>pct_cutoff)
            gp_calc = vertcat(gp_calc, this_other_other_vect(1));
            color_vect_vect = vertcat(color_vect_vect, color_vect(vcell,:));
        end
        if (this_other_other_vect(2)>pct_cutoff)
            gp_calc = vertcat(gp_calc, this_other_other_vect(2));
            color_vect_vect = vertcat(color_vect_vect, color_vect(vcell,:));
        end
        if (this_other_other_vect(3)>pct_cutoff)
            gp_calc = vertcat(gp_calc, this_other_other_vect(3));
            color_vect_vect = vertcat(color_vect_vect, color_vect(vcell,:));
        end
        if (this_other_other_vect(4)>pct_cutoff)
            gp_calc = vertcat(gp_calc, this_other_other_vect(4));
            color_vect_vect = vertcat(color_vect_vect, color_vect(vcell,:));
        end
        if (NC11_flag == 4)
            assert(this_other_other_vect(1) == this_other_other_vect(2));
            assert(this_other_other_vect(3) == this_other_other_vect(4));
        end
    end

    disp(length(gp_calc))
    gp_calc = gp_calc / sum(gp_calc);
    assert(length(gp_calc) == length(color_vect_vect))

    xxx = 1:length(gp_calc);
    [~,IIIIII] = sort(gp_calc);
    datadata = gp_calc(IIIIII);
    datadata = datadata *100;

    figure; hold all;
    for i = 1:length(gp_calc)
        h1111 = bar(xxx(i),datadata(i));
        h1111.FaceColor = color_vect_vect(IIIIII(i),:);
        h1111.EdgeColor = color_vect_vect(IIIIII(i),:);
    end
    ylim([0,20])
    xticklabels([])
    xlim([0.5,max(xxx)+0.5])
    set(gca, 'FontSize', 25.0);

end

tol = 0.01*a_val;
D = sqrt((lattice_pos(:,1) - lattice_pos(:,1)').^2 + (lattice_pos(:,2) - lattice_pos(:,2)').^2);
is_boundary = sum(D > tol & D < a_val + tol, 2) < 6;
for vcell = 1:size(pairs,1)
    if above_thresh_pairs(vcell) && above_thresh_pairs_dist_cutoff(vcell)
        assert(~any(is_boundary(pairs(vcell,:))), 'Kept clone touches lattice boundary');
    end
end


end

%% =======================================================================
%  LOCAL FUNCTIONS
%% =======================================================================

% Assumes lattice_pos is already built as above.
% Returns pairs: an (N/2) x 2 array of indices into lattice_pos.
function pairs = random_hex_matching(lattice_pos, a_val)

N = size(lattice_pos, 1);

% --- Build neighbor list ---
% Two hexagons are neighbors if their distance is a_val (within tolerance)
tol = a_val * 0.01;
neighbors = cell(N, 1);
for i = 1:N
    dists = vecnorm(lattice_pos - lattice_pos(i,:), 2, 2);
    neighbors{i} = find(dists > tol & dists < a_val + tol);
end

% --- Handle odd N: drop the node with fewest neighbors (corner) ---
if mod(N, 2) == 1
    n_neighbors = cellfun(@numel, neighbors);
    [~, drop_idx] = min(n_neighbors);
    fprintf('Odd number of hexagons -- dropping index %d (has %d neighbors)\n', ...
            drop_idx, n_neighbors(drop_idx));
    % Remove from neighbor lists
    for i = 1:N
        neighbors{i} = neighbors{i}(neighbors{i} ~= drop_idx);
    end
    active = setdiff(1:N, drop_idx);
else
    active = 1:N;
end

% --- Greedy random matching with restart ---
max_attempts = 10000000;
pairs = [];

for attempt = 1:max_attempts
    order     = active(randperm(numel(active)));  % random node order
    paired    = false(N, 1);
    pairs_try = zeros(numel(active)/2, 2);
    pair_count = 0;
    success   = true;

    for i = 1:numel(order)
        u = order(i);
        if paired(u), continue; end

        % Find unpaired neighbors in random order
        nbrs_u    = neighbors{u};
        nbrs_free = nbrs_u(~paired(nbrs_u));

        if isempty(nbrs_free)
            success = false;
            break;
        end

        % Pick a random unpaired neighbor
        v = nbrs_free(randi(numel(nbrs_free)));
        paired(u) = true;
        paired(v) = true;
        pair_count = pair_count + 1;
        pairs_try(pair_count, :) = [u, v];
    end

    if success
        pairs = pairs_try;
        fprintf('Matched successfully in %d attempt(s)\n', attempt);
        return;
    end
end

error('Failed to find a perfect matching in %d attempts -- try increasing max_attempts', ...
      max_attempts);
end

function [half1, half2] = clip_polygon_by_line(verts, c, theta)
    half1 = clip_to_halfplane(verts, c, theta,  1);
    half2 = clip_to_halfplane(verts, c, theta, -1);
end

function poly_out = clip_to_halfplane(verts, c, theta, side)
% Sutherland-Hodgman clip of polygon to one side of the directed line.
    n_vec = side * [-sin(theta), cos(theta)];
    sgn   = @(p) (p(1)-c(1))*n_vec(1) + (p(2)-c(2))*n_vec(2);
    nv    = size(verts,1);
    %output = zeros(nv+2, 2);
    output = zeros(2*nv, 2);
    nout   = 0;
    for k = 1:nv
        A_pt = verts(k,:);
        B_pt = verts(mod(k,nv)+1,:);
        sA   = sgn(A_pt);
        sB   = sgn(B_pt);
        if sA >= 0
            nout = nout+1; output(nout,:) = A_pt;
            if sB < 0
                nout = nout+1;
                output(nout,:) = seg_line_intersect(A_pt, B_pt, c, theta);
            end
        else
            if sB >= 0
                nout = nout+1;
                output(nout,:) = seg_line_intersect(A_pt, B_pt, c, theta);
            end
        end
    end
    poly_out = output(1:nout,:);
end

function pt = seg_line_intersect(A, B, c, theta)
    d     = [cos(theta), sin(theta)];
    AB    = B - A;
    denom = d(1)*(-AB(2)) - d(2)*(-AB(1));
    if abs(denom) < 1e-14, pt = (A+B)/2; return; end
    rhs = A - c;
    s   = (d(1)*rhs(2) - d(2)*rhs(1)) / denom;
    pt  = A + s*AB;
end

% integrate_gaussian_polygon(A, mu, sigma, verts, ABSTOL, RELTOL)
function I = integrate_gaussian_polygon(A, mu, sigma, verts, abstol, reltol)
% Integrate Gaussian over polygon via signed triangle fan from centroid.
    P0 = mean(verts, 1);
    n  = size(verts, 1);
    I  = 0;
    for k = 1:n
        k2 = mod(k,n)+1;
        I  = I + signed_triangle_integral(A, mu, sigma, P0, ...
                     verts(k,:), verts(k2,:), abstol, reltol);
    end
end

function I = signed_triangle_integral(A, mu, sigma, P0, P1, P2, abstol, reltol)
% Signed Gaussian integral over triangle (P0,P1,P2) via adaptive 2D quadrature.
    e1 = P1-P0;  e2 = P2-P0;
    J  = e1(1)*e2(2) - e1(2)*e2(1);   % signed cross product
    if abs(J) < 1e-14, I = 0; return; end

    d  = P0 - mu;
    s2 = 2*sigma^2;

    % Early exit: triangle is negligibly far from Gaussian peak.
    rmax = max(norm(e1), norm(e2));   % triangle lies in ball of radius rmax about P0
    dmin = max(0, norm(d) - rmax);
    if A*abs(J)/2 * exp(-dmin^2/s2) < abstol*1e-3, I = 0; return; end

    % Gaussian evaluated directly in triangle coordinates (values in [0,1])
    f = @(s,t) exp(-((d(1) + s*e1(1) + t*e2(1)).^2 + ...
                     (d(2) + s*e1(2) + t*e2(2)).^2) / s2);

    I = A * J * integral2(f, 0, 1, 0, @(s) 1 - s, ...
                          'AbsTol', abstol, 'RelTol', reltol);
end

function [vx, vy] = hexagon_vertices(R, c, phi)
    angles = phi + (0:5)*pi/3;
    vx = c(1) + R*cos(angles);
    vy = c(2) + R*sin(angles);
end
