%% Use it For positive curvature--------------------------------------------------------------------------
% AFP GENERATION - UNIFIED BINARY MASK WITH DUAL DIMENSIONS
% Handles Positive/Negative Curvature + Separate Canvas Trimming
% Includes Automatic Flipped Image Generation for Opposite Curvature
%--------------------------------------------------------------------------
clear; clc;
warning('off', 'MATLAB:polyshape:repairedBySimplify');

% -------------------------------------------------------------------------
% 1. USER INPUTS
% -------------------------------------------------------------------------
% --- PHYSICAL PLATE DIMENSIONS (For Curvature Math & Physics) ---
plate_width = 330;     % Physical plate width (mm)
plate_height = 580;    % Physical plate height (mm)

% --- IMAGE EXTRACTION / CANVAS DIMENSIONS (Frame size in saved image) ---
canvas_width = 254;    % Desired extraction image width (mm)
canvas_height = 406;   % Desired extraction image height (mm)

% --- STEERING & ANGLE PARAMETERS ---
T0_deg = 41;           % Center angle (degrees) from Y-axis
T1_deg = 53.714;         % Edge angle (degrees) from Y-axis

% --- AFP MACHINE PARAMETERS ---
Nt = 8;                % Number of tows in a course
tw = 3.175;            % Tow width (mm) - 1/8 inch
w = Nt * tw;           % Total course width (mm)
tow_edges = linspace(-w/2, w/2, Nt + 1);

% _________________________________________________________________________
% Extended Dimensions to ensure math curves reach past all boundaries
new_plate_width = 1.5 * max(plate_width, canvas_width);
new_plate_height = 1.5 * max(plate_height, canvas_height);

% --- DIRECTORY & FILENAME SETUP ---
folder_name = 'AFP_Defect_image';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end
% Filenames for both Positive and Negative variants
save_filename_pos = fullfile(folder_name, sprintf('BinaryMask_pos_T0_%g_T1_%g.png', T0_deg, T1_deg));
save_filename_neg = fullfile(folder_name, sprintf('BinaryMask_neg_T0_%g_T1_%g.png', T0_deg, T1_deg));

%--------------------------------------------------------------------------
% 2. DYNAMIC CALCULATION OF KAPPA AND P (BASED ON PHYSICAL PLATE)
%--------------------------------------------------------------------------
T0 = deg2rad(T0_deg);
T1 = deg2rad(T1_deg);
if T0_deg > T1_deg
    kappa = (cos(T1) - cos(T0)) / (plate_width / 2);
elseif T1_deg > T0_deg
    kappa = (cos(T0) - cos(T1)) / (plate_width / 2);
else
    kappa = 0; 
end

if kappa == 0, p = inf; else, p = 1 / kappa; end
x = linspace(-new_plate_width/2, new_plate_width/2, 2500);

%--------------------------------------------------------------------------
% 3. CALCULATE BASE PATH AND CONSTANT VERTICAL SHIFT (w_e)
%--------------------------------------------------------------------------
if isinf(p)
    cos_theta = ones(size(x)) * cos(T0);
    theta = acos(cos_theta);
    y_base = x * cot(T0); 
    w_e = w / cos(T0); 
else
    if T0_deg > T1_deg
        cos_theta = cos(T0) + kappa .* abs(x);
    elseif T1_deg > T0_deg
        cos_theta = cos(T0) - kappa .* abs(x);
    end
    
    cos_limit = cos(deg2rad(0.01)); 
    cos_theta = max(min(cos_theta, cos_limit), -cos_limit); 
    
    theta = acos(cos_theta); 
    dy_dx = cot(theta);
    y_base = cumtrapz(x, dy_dx);
    y_base = y_base - interp1(x, y_base, 0); 
    
    rho = abs(p);
    inner_term = max(0, (rho + w/2)^2 - (rho * cos(T0))^2); %here +w/2
    w_e = 2 * (sqrt(inner_term) - rho * sin(T0)); %sign changed here
end

%--------------------------------------------------------------------------
% 4. BUILD BASE MASTER CURVES
%--------------------------------------------------------------------------
X_RED = x + (w/2) .* cos(theta);
Y_RED = y_base - (w/2) .* sin(theta);

use_Y = abs(Y_RED(end) - Y_RED(1)) > abs(X_RED(end) - X_RED(1));
if use_Y
    [MASTER_VALS, sort_idx] = sort(Y_RED);
    MASTER_SEC = X_RED(sort_idx);
else
    [MASTER_VALS, sort_idx] = sort(X_RED);
    MASTER_SEC = Y_RED(sort_idx);
end
[UNI_MASTER_VALS, id_uni] = unique(MASTER_VALS);
UNI_MASTER_SEC = MASTER_SEC(id_uni);

k_right = -1; 
y_c_right = y_base + k_right * w_e;
TOW_X_CELL = cell(1, Nt + 1);
TOW_Y_CELL = cell(1, Nt + 1);
for i = 1:(Nt + 1)
    TOW_X_CELL{i} = x + tow_edges(i) .* cos(theta);
    TOW_Y_CELL{i} = y_c_right - tow_edges(i) .* sin(theta);
end

%--------------------------------------------------------------------------
% 5. FIND INTERSECTIONS & TANGENTS (Only populates during Overlaps)
%--------------------------------------------------------------------------
INTERSECTIONS = []; 
for i = 1:(Nt + 1)
    x_edge = TOW_X_CELL{i};
    y_edge = TOW_Y_CELL{i};
    
    if length(UNI_MASTER_VALS) > 1 
        if use_Y
            SEC_bound = interp1(UNI_MASTER_VALS, UNI_MASTER_SEC, y_edge, 'linear', NaN);
            d_diff = x_edge - SEC_bound; 
        else
            SEC_bound = interp1(UNI_MASTER_VALS, UNI_MASTER_SEC, x_edge, 'linear', NaN);
            d_diff = y_edge - SEC_bound; 
        end
        
        crossings = find(d_diff(1:end-1) .* d_diff(2:end) <= 0 & ~isnan(d_diff(1:end-1)) & ~isnan(d_diff(2:end)));
        for c_idx = 1:length(crossings)
            idx = crossings(c_idx);
            t = d_diff(idx) / (d_diff(idx) - d_diff(idx+1));
            x_int = x_edge(idx) + t * (x_edge(idx+1) - x_edge(idx));
            y_int = y_edge(idx) + t * (y_edge(idx+1) - y_edge(idx));
            theta_int = theta(idx) + t * (theta(idx+1) - theta(idx));
            INTERSECTIONS = [INTERSECTIONS; x_int, y_int, theta_int, i];
        end
        
        tangent_tol = 0.5; 
        for k = 2:(length(d_diff)-1)
            if sign(d_diff(k-1)) == sign(d_diff(k+1)) && sign(d_diff(k)) == sign(d_diff(k-1)) && ~isnan(d_diff(k))
                if abs(d_diff(k)) < abs(d_diff(k-1)) && abs(d_diff(k)) < abs(d_diff(k+1)) && abs(d_diff(k)) < tangent_tol
                    x_tan = x_edge(k); y_tan = y_edge(k); theta_tan = theta(k);
                    INTERSECTIONS = [INTERSECTIONS; x_tan, y_tan, theta_tan, i];
                end
            end
        end
    end
end

%--------------------------------------------------------------------------
% 6. DYNAMIC SORTING & EXTRACTING TOW-DROP POLYGONS
%--------------------------------------------------------------------------
INT_PTS_NEG = cell(1, Nt + 1); 
INT_PTS_POS = cell(1, Nt + 1); 
if ~isempty(INTERSECTIONS)
    for i = 1:(Nt + 1)
        tow_mask = (INTERSECTIONS(:, 4) == i);
        tow_data = INTERSECTIONS(tow_mask, :);
        if ~isempty(tow_data)
            x_vals = tow_data(:, 1);
            [~, sort_idx] = sort(x_vals);
            tow_data_sorted = tow_data(sort_idx, :);
            
            r_neg = tow_data_sorted(1, :);
            INT_PTS_NEG{i} = [r_neg(1), r_neg(2), r_neg(1) + tw*cos(r_neg(3)), r_neg(2) - tw*sin(r_neg(3))];
            
            r_pos = tow_data_sorted(end, :);
            INT_PTS_POS{i} = [r_pos(1), r_pos(2), r_pos(1) + tw*cos(r_pos(3)), r_pos(2) - tw*sin(r_pos(3))];
        end
    end
end

BASE_GAPS = {}; 
for gap_side = 1:2
    if gap_side == 1, CURRENT_INT_PTS = INT_PTS_NEG; else, CURRENT_INT_PTS = INT_PTS_POS; end
    
    for i = 1:Nt
        if ~isempty(CURRENT_INT_PTS{i}) && ~isempty(CURRENT_INT_PTS{i+1})
            pt1 = CURRENT_INT_PTS{i}; 
            x_int1 = pt1(1); y_int1 = pt1(2); x_norm_end = pt1(3); y_norm_end = pt1(4); 
            
            pt2 = CURRENT_INT_PTS{i+1}; 
            x_int2 = pt2(1); y_int2 = pt2(2); 
            
            X_next = TOW_X_CELL{i+1}; Y_next = TOW_Y_CELL{i+1};
            [~, id_t1] = min((X_next - x_norm_end).^2 + (Y_next - y_norm_end).^2);
            [~, id_t2] = min((X_next - x_int2).^2 + (Y_next - y_int2).^2);
            if id_t1 <= id_t2, range_t = id_t1:id_t2; else, range_t = id_t1:-1:id_t2; end
            
            arc_tow_x = X_next(range_t); arc_tow_y = Y_next(range_t);
            arc_tow_x(1) = x_norm_end; arc_tow_y(1) = y_norm_end; 
            arc_tow_x(end) = x_int2;   arc_tow_y(end) = y_int2;
            
            [~, id_b1] = min((X_RED - x_int2).^2 + (Y_RED - y_int2).^2);
            [~, id_b2] = min((X_RED - x_int1).^2 + (Y_RED - y_int1).^2);
            if id_b1 <= id_b2, range_b = id_b1:id_b2; else, range_b = id_b1:-1:id_b2; end
            
            arc_blue_x = X_RED(range_b); arc_blue_y = Y_RED(range_b);
            arc_blue_x(1) = x_int2; arc_blue_y(1) = y_int2; 
            arc_blue_x(end) = x_int1; arc_blue_y(end) = y_int1;
            
            BASE_GAPS{end+1} = struct('poly_x', [arc_tow_x, arc_blue_x], 'poly_y', [arc_tow_y, arc_blue_y]);
        end
    end
end

%--------------------------------------------------------------------------
% 7. UNIFIED BINARY MASK RENDERING (TRIMMED TO CANVAS)
%--------------------------------------------------------------------------
aspect_ratio = canvas_height / canvas_width;
fig_width = 800; 
fig_height = round(fig_width * aspect_ratio);

fig = figure('Color','k', 'Name', 'AFP - Binary Defect Mask', 'Visible', 'off', ...
             'Position', [100, 100, fig_width, fig_height]); 

ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0 0 1 1]);
hold(ax, 'on');
axis(ax, 'equal', 'off');
% Restrict the visual axes strictly to the CANVAS dimensions
xlim(ax, [-canvas_width/2, canvas_width/2]); 
ylim(ax, [-canvas_height/2, canvas_height/2]);

% STEP A: Paint the CANVAS background BLACK (Handles natural separation gaps)
fill(ax, [-canvas_width/2, canvas_width/2, canvas_width/2, -canvas_width/2], ...
         [-canvas_height/2, -canvas_height/2, canvas_height/2, canvas_height/2], ...
         'k', 'EdgeColor', 'none');

% Define the outer boundary polygon of a single full course
X_left_course = x - (w/2) .* cos(theta);
Y_left_course = y_base + (w/2) .* sin(theta);
poly_course_X = [X_RED, fliplr(X_left_course)];

% Ensure we shift enough times to cover whichever is taller (plate or canvas)
max_shifts = ceil((max(plate_height, canvas_height)/2 + abs(min(Y_RED))) / w_e) + 5; 

% STEP B: Draw all courses as solid WHITE 
for k = -max_shifts : max_shifts
    Y_shift = k * w_e;
    poly_course_Y = [Y_RED + Y_shift, fliplr(Y_left_course + Y_shift)];
    
    % Draw the solid white course
    fill(ax, poly_course_X, poly_course_Y, 'w', 'FaceAlpha', 1.0, 'EdgeColor', 'none');
end

% STEP C: Carve out the explicit Tow-Drop defects in BLACK (For overlapping cases)
for k = -max_shifts : max_shifts
    Y_shift = k * w_e;
    
    for p_idx = 1:length(BASE_GAPS)
        g = BASE_GAPS{p_idx};
        poly_y_shifted = g.poly_y + Y_shift;
        
        % Only render if the black triangle sits inside the CANVAS frame
        if max(g.poly_x) >= -canvas_width/2 && min(g.poly_x) <= canvas_width/2 && ...
           max(poly_y_shifted) >= -canvas_height/2 && min(poly_y_shifted) <= canvas_height/2
            
            fill(ax, g.poly_x, poly_y_shifted, 'k', 'FaceAlpha', 1.0, 'EdgeColor', 'none');
        end
    end
end

%--------------------------------------------------------------------------
% 8. SAVE PERFECTLY TRIMMED IMAGES (POSITIVE & NEGATIVE)
%--------------------------------------------------------------------------
% Use print to strictly obey the window dimensions without auto-cropping
set(fig, 'PaperPositionMode', 'auto'); 
set(fig, 'InvertHardcopy', 'off'); 

% 1. Print and save the Positive Curvature image
print(fig, save_filename_pos, '-dpng', '-r300'); 

% 2. Read, Flip, and Save the Negative Curvature image
img_original = imread(save_filename_pos);
img_flipped = fliplr(img_original);
imwrite(img_flipped, save_filename_neg);

fprintf('Process Complete.\n');
fprintf('Trimmed Positive image saved as: %s\n', save_filename_pos);
fprintf('Trimmed Negative image saved as: %s\n', save_filename_neg);

% --- VISUALIZATION OVERLAY ---
set(ax, 'Position', [0.1 0.1 0.8 0.8]); 
axis(ax, 'on'); 
grid(ax, 'on');
set(ax, 'Color', 'k'); % Keep preview background black

% Zoom out slightly to show both boundaries
xlim(ax, [-new_plate_width/2, new_plate_width/2]);
ylim(ax, [-new_plate_height/2, new_plate_height/2]);

% Red Dashed = Actual Physical Cut Plate
rectangle(ax, 'Position', [-plate_width/2, -plate_height/2, plate_width, plate_height], ...
    'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
% Green Dot-Dash = Image Canvas Extraction Boundary
rectangle(ax, 'Position', [-canvas_width/2, -canvas_height/2, canvas_width, canvas_height], ...
    'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '-.');

title(ax, sprintf('Unified Binary Mask: T_0 = %g^\\circ, T_1 = %g^\\circ\nRed = Physical Plate | Green = Extracted Canvas', T0_deg, T1_deg), 'Color', 'w');
set(fig, 'Visible', 'on'); 
hold(ax, 'off');