%% Use + B/w image + Area compute  ======================================================================
% Doesnot trim the black slight region from all corners for negative angles
% VAT DEFECT GENERATION BASED ON SECTION 2.1 (Zamani et al. 2025)
% BLACK DEFECTS, BORDERLESS EXPORT, AND PIXEL-BASED AREA CALCULATION
% ======================================================================
clear all; close all; clc;
tic
% 1. GLOBAL PARAMETERS (in m)
a = 0.254;
Tow_Width = 3.125e-3;
Num_Tows = 16;
w_course = Num_Tows * Tow_Width;
p_l = w_course / 2;
p_r = w_course / 2;

% --- NEW USER PARAMETERS ---
line_thickness_calc = 0.39;  % (mm) Thickness for area math
line_thickness_plot = 1.0;   % (pt) Thickness for visual display

%% If T0<T1, Edge=Gap,Center=Overlap...If T0>T1, Edge=Overlap,Center=Gap
% Define your cases [T0_deg, T1_deg, Tangency_Type]
cases = {
60,  15, 'Center';    % Case 1 (Positive)
-60, -15, 'Center'     % Case 2 (Negative)
};

% --- EXTENDED DOMAIN ---
x_buffer = w_course * 2;
x = linspace(-a/2 - x_buffer, a/2 + x_buffer, 2000);

% Loop through each case to generate a separate window
for i = 1:size(cases,1)
T0_deg = cases{i,1};
T1_deg = cases{i,2};
tangency = cases{i,3};

fig_name = sprintf('Case %d: T0=%d, T1=%d (%s Tangency)', i, T0_deg, T1_deg, tangency);
fig = figure('Name', fig_name, 'Position', [50 + i*50, 50 + i*50, 700, 700], 'Color', 'w');

% Force print to respect exact figure dimensions and keep background colors
set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');

% Use exact axes handles to guarantee clean, borderless exports
ax = axes('Position', [0 0 1 1]); % Force axes to fill the figure
hold(ax, 'on'); 
axis(ax, 'equal'); 
axis(ax, 'off'); % Silences the grid, axes, and background

% Clip the view exactly to the [-0.25, 0.25] square plate dimensions
xlim(ax, [-a/2, a/2]); 
ylim(ax, [-a/2, a/2]);

% --- OFF-WHITE BOUNDING BOX ---
% By using a color just 1 bit off from pure white [0.99 0.99 0.99], 
% we force the bounding box to the exact mathematical corners of the plate.
% Since it is > 128 grayscale, the defect area calculator completely ignores it.
plot(ax, [-a/2, a/2, a/2, -a/2, -a/2], [-a/2, -a/2, a/2, a/2, -a/2], ...
    'Color', [1 1 1], 'LineWidth', 0.1, 'HandleVisibility', 'off');

T0 = deg2rad(T0_deg);
T1 = deg2rad(T1_deg);

% Fiber Orientation Theta(x)
theta = T0 + ((T1 - T0) / (a/2)) * abs(x);
theta = max(-pi/2 + 1e-4, min(pi/2 - 1e-4, theta));

% Centerline Path y_CL(x)
y_CL = zeros(size(x));
idx_pos = x >= 0;
idx_neg = x < 0;

if abs(T1 - T0) > 1e-6
    k = (T1 - T0) / (a/2);
    y_CL(idx_pos) = (-1 / k) * (log(cos(T0 + k*x(idx_pos))) - log(cos(T0)));
    y_CL(idx_neg) = (1 / k) * (log(cos(T0 - k*x(idx_neg))) - log(cos(T0)));
else
    y_CL = tan(T0) * x;
end

% Shifting Logic
if strcmp(tangency, 'Center')
    w_e = w_course / cos(T0); 
else
    w_e = w_course / cos(T1); 
end

w_proj = w_course ./ cos(theta);
Num_Shifts = ceil((1.5 * a) / w_e) + 10;
shifts = -Num_Shifts : Num_Shifts;

% Store line handles so we can change their thickness programmatically
line_handles = [];

% Plotting loop
for k = 1:length(shifts)-1
    shift_current = shifts(k) * w_e;
    shift_next = shifts(k+1) * w_e;
    
    y_center_current = y_CL + shift_current;
    y_center_next = y_CL + shift_next;
    
    x_l_current = x - p_l * sin(theta);
    y_l_current = y_center_current + p_l * cos(theta);
    
    x_r_next = x + p_r * sin(theta);
    y_r_next = y_center_next - p_r * cos(theta);
    
    poly_x = [x_l_current, fliplr(x_r_next)];
    poly_y = [y_l_current, fliplr(y_r_next)];
    
    gap_magnitude = w_e - w_proj;
    
    % Shade defects BLACK
    if abs(mean(gap_magnitude)) > 1e-5
        fill(ax, poly_x, poly_y, 'k', 'FaceAlpha', 1.0, 'EdgeColor', 'none');
    end
    
    % Plot boundary curves BLACK and collect handles
    h1 = plot(ax, x_l_current, y_l_current, 'k-'); 
    h2 = plot(ax, x_r_next, y_r_next, 'k-');       
    line_handles = [line_handles, h1, h2];
end

%% ==================================================================
% AREA CALCULATION (PIXEL-BASED CORRELATION METHOD)
% ==================================================================
% 1. Convert mm to standard MATLAB points (1 pt = 1/72 inch = 0.3528 mm)
lw_calc_pt = line_thickness_calc * (72 / 25.4);

% 2. Set lines to mathematical thickness
set(line_handles, 'LineWidth', lw_calc_pt);

% 3. Export high-res temporary image of the exact plate boundaries using print
temp_img_name = 'temp_calc.png';
print(fig, temp_img_name, '-dpng', '-r600');

% 4. Read image and calculate ratio of black pixels (defects + lines) to total pixels
img = imread(temp_img_name);
if size(img,3) == 3
    img = rgb2gray(img);
end

bw = img < 128; % Black pixels evaluate to True (1)
defect_percentage = (sum(bw(:)) / numel(bw)) * 100;

% Output strictly to the console
fprintf('Case %d [%d, %d] %s Tangency: Total Defect Area = %.2f%%\n', ...
    i, T0_deg, T1_deg, tangency, defect_percentage);

% Delete the temporary calculation image
delete(temp_img_name);

%% ==================================================================
% VISUAL EXPORT
% ==================================================================
% Update lines to the visual display thickness
set(line_handles, 'LineWidth', line_thickness_plot);

% Save the final borderless file using print
final_img_name = sprintf('LVF_%d_%d_%d_%s.png', a*1000, T0_deg, T1_deg, tangency);
print(fig, final_img_name, '-dpng', '-r300');
fprintf('Saved borderless visual image to: %s\n\n', final_img_name);

% Add the title ONLY after saving, so the saved image remains perfectly borderless
title(ax, sprintf('[\\pm\\langle%d, %d\\rangle] %s Tangency | Area: %.2f%%', ...
    T0_deg, T1_deg, tangency, defect_percentage), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');


end
toc