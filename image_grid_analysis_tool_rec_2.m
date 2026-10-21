function image_grid_analysis_tool_rec_2()
% IMAGE GRID ANALYSIS TOOL - MATLAB Version
% Modified to handle custom Row/Column inputs, generating rectangular 
% blocks that align perfectly with rectangular FEA meshes.
    try
        inputs = get_user_input();
        success = process_image_with_grid(inputs);
        
        if success
            display_summary(inputs.output_excel_path, inputs.original_filename);
            
            fprintf('\n🎉 PROCESSING COMPLETED SUCCESSFULLY!\n');
            fprintf('📄 Original file: %s\n', inputs.original_filename);
            fprintf('🖼️  Output Image: %s\n', inputs.output_image_name);
            fprintf('📊 Excel Report: %s (Includes "Summary" sheet)\n', inputs.output_excel_name);
            fprintf('📁 Location: %s\n', inputs.output_folder);
            
            % Ask if user wants to open results
            open_results = input('\n🔍 Open results folder? (y/n): ', 's');
            if lower(open_results) == 'y'
                winopen(inputs.output_folder);
            end
        else
            fprintf('\n❌ Processing failed.\n');
        end
        
    catch ME
        if strcmp(ME.identifier, 'MATLAB:operation_cancelled')
            fprintf('\n⚠️  Operation cancelled by user.\n');
        else
            fprintf('\n❌ Unexpected error: %s\n', ME.message);
        end
    end
end

function inputs = get_user_input()
    fprintf('--- IMAGE GRID ANALYSIS TOOL ---\n');
    
    while true
        image_path = input('\n📁 Enter image path: ', 's');
        image_path = strtrim(image_path);
        if exist(image_path, 'file')
            break;
        else
            fprintf('❌ File not found!\n');
        end
    end
    
    fprintf('\n📏 Select Meshing Density:\n');
    fprintf('  1. Coarse  (Auto Square blocks)\n');
    fprintf('  2. Medium  (Auto Square blocks)\n');
    fprintf('  3. Fine    (Auto Square blocks)\n');
    fprintf('  4. Custom  (Specify exact number of rows and columns)\n');
    
    while true
        choice = input('Enter choice (1-4): ');
        if ismember(choice, [1, 2, 3, 4])
            inputs.mesh_choice = choice;
            if choice == 4
                inputs.custom_rows = input('Enter required number of ROWS (n): ');
                inputs.custom_cols = input('Enter required number of COLUMNS (m): ');
            else
                inputs.custom_rows = 0;
                inputs.custom_cols = 0;
            end
            break;
        else
            fprintf('❌ Invalid choice. Please enter 1, 2, 3, or 4.\n');
        end
    end
    
    output_folder = input('\n📂 Output folder (Enter for image folder): ', 's');
    output_folder = strtrim(output_folder);
    if isempty(output_folder)
        output_folder = fileparts(image_path);
        if isempty(output_folder), output_folder = '.'; end
    end
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    
    [~, name, ext] = fileparts(image_path);
    output_image_name = [name '_grid_output.png'];
    output_excel_name = [name '_grid_analysis.xlsx'];
    
    inputs.image_path = image_path;
    inputs.output_folder = output_folder;
    inputs.output_image_path = fullfile(output_folder, output_image_name);
    inputs.output_excel_path = fullfile(output_folder, output_excel_name);
    inputs.output_image_name = output_image_name;
    inputs.output_excel_name = output_excel_name;
    inputs.original_filename = [name ext];
end

function success = process_image_with_grid(inputs)
    fprintf('\n🔄 PROCESSING IMAGE...\n');
    
    try
        % Read image, colormap, and alpha channel if present
        [image, map, ~] = imread(inputs.image_path);
        
        % Handle indexed images (GIFs, some PNGs)
        if ~isempty(map)
            image = ind2rgb(image, map);
            image = im2uint8(image);
        end
        
        % Handle Bit-depth / Channels
        if ndims(image) == 3
            if size(image, 3) == 4
                fprintf('🖼️ 32-bit RGBA image detected. Dropping alpha channel...\n');
                image = image(:, :, 1:3); % Drop alpha to make it standard RGB
            end
            if size(image, 3) == 3
                fprintf('🎨 Color image detected. Converting to grayscale...\n');
                image = rgb2gray(image);
            end
        end
        
    catch ME
        fprintf('❌ Error loading image: %s\n', ME.message);
        success = false; 
        return;
    end
    
    [orig_height, orig_width] = size(image);
    max_dim = max(orig_height, orig_width);
    
    % Determine grid density based on user mesh choice
    if inputs.mesh_choice == 4
        % CUSTOM ROWS AND COLS
        grid_rows = inputs.custom_rows;
        grid_cols = inputs.custom_cols;
    else
        % AUTO SQUARE BLOCKS
        switch inputs.mesh_choice
            case 1
                block_size = max(1, floor(max_dim / 15)); % Coarse
            case 2
                block_size = max(1, floor(max_dim / 35)); % Medium
            case 3
                block_size = max(1, floor(max_dim / 75)); % Fine
        end
        grid_rows = floor(orig_height / block_size);
        grid_cols = floor(orig_width / block_size);
    end
    
    if grid_rows == 0 || grid_cols == 0
        fprintf('❌ Error: Requested grid is too large for the image resolution.\n');
        success = false;
        return;
    end
    
    % =====================================================================
    % THE FIX: Use linspace to calculate exact pixel boundaries without 
    % dropping remainders. NO CROPPING ALLOWED.
    % =====================================================================
    x_edges = round(linspace(1, orig_width + 1, grid_cols + 1));
    y_edges = round(linspace(1, orig_height + 1, grid_rows + 1));
    
    fprintf('📄 Image: %s\n', inputs.original_filename);
    fprintf('📊 Dimensions: %d x %d pixels (W x H)\n', orig_width, orig_height);
    fprintf('🔲 Grid Density: %d rows x %d cols\n', grid_rows, grid_cols);
    
    % Create output image matching original dimensions perfectly
    output_image = ones(orig_height, orig_width, 'uint8') * 255;
    excel_data = cell(grid_rows * grid_cols, 7);
    total_blocks = grid_rows * grid_cols;
    
    fprintf('\n⏳ Analyzing blocks silently. Please wait...\n');
    
    for row = 1:grid_rows
        for col = 1:grid_cols
            % Calculate exact block boundaries
            start_x = x_edges(col);
            end_x   = x_edges(col+1) - 1;
            start_y = y_edges(row);
            end_y   = y_edges(row+1) - 1;
            
            % Extract block and analyze
            block = image(start_y:end_y, start_x:end_x);
            total_pixels = numel(block);
            
            % Thresholding for black pixels (adjust 128 if needed)
            binary_block = block < 128;
            black_pixels = sum(binary_block(:));
            black_percentage = black_pixels / total_pixels;
            white_percentage = 1 - black_percentage;
            grayscale_value = uint8(white_percentage * 255);
            
            % Fill output block perfectly back into the uncropped matrix
            output_image(start_y:end_y, start_x:end_x) = grayscale_value;
            
            % Store data for Excel
            idx = (row-1) * grid_cols + col;
            excel_data{idx, 1} = row;
            excel_data{idx, 2} = col;
            excel_data{idx, 3} = round(black_percentage, 4);
            excel_data{idx, 4} = round(white_percentage, 4);
            excel_data{idx, 5} = grayscale_value;
            excel_data{idx, 6} = total_pixels;
            excel_data{idx, 7} = black_pixels;
        end
    end
    
    fprintf('--------------------------------------------------\n');
    
    imwrite(output_image, inputs.output_image_path);
    fprintf('✅ Image saved: %s\n', inputs.output_image_name);
    
    % =====================================================================
    % EXCEL GENERATION: SUMMARY & RAW DATA TABS
    % =====================================================================
    column_names = {'Block_Row', 'Block_Column', 'Black_Percentage', ...
                   'White_Percentage', 'Grayscale_Value', 'Total_Pixels', 'Black_Pixels'};
    data_table = cell2table(excel_data, 'VariableNames', column_names);
    
    total_black_pixels = sum(data_table.Black_Pixels);
    total_blocks_analyzed = total_blocks;
    
    % Use mean of total pixels per block for the summary to account for slight linspace variations
    average_pixels_per_block = mean(data_table.Total_Pixels);
    
    average_black_per_block = total_black_pixels / total_blocks_analyzed;
    overall_black_percentage = (average_black_per_block / average_pixels_per_block) * 100;
    
    % Build the Summary Table
    Summary_Metric = {
        'Total Blocks Analyzed'; 
        'Average Pixels per Block';
        'Total Black Pixels (Sum)'; 
        'Average Black Pixels per Block'; 
        'Overall Black Percentage (%)'
    };
    Value = [
        total_blocks_analyzed;
        round(average_pixels_per_block, 2);
        total_black_pixels; 
        average_black_per_block; 
        overall_black_percentage
    ];
    summary_table = table(Summary_Metric, Value);
    
    try
        writetable(summary_table, inputs.output_excel_path, 'Sheet', 'Summary');
        writetable(data_table, inputs.output_excel_path, 'Sheet', 'Block_Data');
    catch
        writetable(data_table, inputs.output_excel_path);
    end
    
    fprintf('✅ Excel saved: %s\n', inputs.output_excel_name);
    
    success = true;
end

function display_summary(excel_path, original_filename)
    % Read the Block_Data sheet specifically for the terminal summary
    try
        data = readtable(excel_path, 'Sheet', 'Block_Data');
    catch
        data = readtable(excel_path);
    end
    
    overall_black_pct = (sum(data.Black_Pixels) / sum(data.Total_Pixels)) * 100;
    
    fprintf('\n📈 ANALYSIS SUMMARY\n');
    fprintf('📄 Original file: %s\n', original_filename);
    fprintf('📊 Total blocks analyzed: %d\n', height(data));
    fprintf('⚫ Overall Black Percentage: %.2f%%\n', overall_black_pct);
    fprintf('🎨 Grayscale range: %d to %d\n', min(data.Grayscale_Value), max(data.Grayscale_Value));
    fprintf('📏 Total pixels analyzed: %d\n', sum(data.Total_Pixels));
    fprintf('🔲 Black pixels total: %d\n', sum(data.Black_Pixels));
end