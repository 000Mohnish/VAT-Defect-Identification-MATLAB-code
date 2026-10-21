# User Manual: Image-Based VAT Defect Analysis

This manual outlines the step-by-step procedure for generating Automated Fiber Placement (AFP) manufacturing defects in Variable Angle Tow (VAT) composite laminates and exporting the resulting area fraction data for finite element analysis.

## 1. Software Requirements
*    MATLAB (Tested on R2024b or newer).


## 2. Workflow Overview

The analysis is split into two distinct phases:
1.  **Phase 1: Defect Generation.** Using the mathematical definition of the fiber path (Constant Curvature or Linearly Varying Fiber), the code generates a highly precise binary mask where the black area represents the defect (gap/overlap) and the white area represents the pristine composite region.
2.  **Phase 2: Grid Analysis.** The binary image is passed through a grid analysis tool, which breaks the plate into rectangular blocks matching your FEA mesh. It computes the area percentage of the defect in each element and exports the results to an Excel spreadsheet.

---

## 3. Phase 1: Defect Generation

Depending on the specific steering kinematics, choose one of the three generation scripts.

### 3.1 Script Selection
*   **`CC_negative_curvature.m`**: Use when operating with a constant curvature path where the center angle is greater than the edge angle ($T_0 > T_1$).
*   **`CC_positive_curvature.m`**: Use when operating with a constant curvature path where the edge angle is greater than the center angle ($T_1 > T_0$).
*   **`LVF_new_area2.m`**: Use when operating with a Linearly Varying Fiber path. This script generates borderless exports and pixel-based area calculations for specific tangency conditions.

### 3.2 Key Variables to Configure Before Running
Open the chosen script in MATLAB and modify the variables under the `USER INPUTS` section:

| Variable | Description |
| :--- | :--- |
| `plate_width` / `plate_height` | Actual physical dimensions of the composite plate in mm. |
| `canvas_width` / `canvas_height` | Dimensions of the final extracted image frame in mm. |
| `T0_deg` | The fiber angle at the center of the plate (degrees). |
| `T1_deg` | The fiber angle at the edges of the plate (degrees). |
| `Nt` | Number of tows in a course (default is 8 for CC, 16 for LVF). |
| `tw` | Tow width in mm (default is 3.175 mm or 1/8 inch). |
| `line_thickness_calc` | Line thickness used for area calculation math (mm). |

*Note for LVF Users: You will need to define your configurations in the `cases` cell array, specifying `[T0_deg, T1_deg, Tangency_Type]` (e.g., `'Center'`).*

### 3.3 Expected Outputs
Running these scripts will automatically create a folder named `AFP_Defect_image_working` (or `AFP_Defect_image`). Inside, you will find:
*   A "Positive" trimmed image (e.g., `Gap_pos_T0_45_T1_26.png`).
*   An automatically generated "Negative" flipped image.
*   The MATLAB console will print the total defect area and percentage relative to the physical plate size.

---

## 4. Phase 2: Grid Analysis

Once the binary images are generated, run the `image_grid_analysis_tool_rec_2.m` function. This tool processes the image interactively.

### 4.1 Running the Grid Analysis
1.  Run the function in the MATLAB command window: `image_grid_analysis_tool_rec_2()`
2.  **Input the Image Path:** When prompted (`📁 Enter image path:`), type or paste the full file path to the binary image generated in Phase 1.
3.  **Select Meshing Density:** You will be prompted to choose the grid density to match your FEA mesh:
    *   `1`: Coarse (Auto Square blocks)
    *   `2`: Medium (Auto Square blocks)
    *   `3`: Fine (Auto Square blocks)
    *   `4`: Custom (You will be prompted to enter the exact number of rows and columns).
4.  **Define Output Folder:** Press Enter to save the output files in the same directory as the original image, or type a custom path.

### 4.2 Expected Outputs
The script will output two files:
1.  **Grayscale Image:** (e.g., `filename_grid_output.png`) A patch-based image where each element's grayscale intensity corresponds to the percentage of pristine material (white = 100% pristine).
2.  **Excel Spreadsheet:** (e.g., `filename_grid_analysis.xlsx`) Contains two sheets:
    *   **`Summary`:** High-level metrics, including the overall black percentage and total blocks analyzed.
    *   **`Block_Data`:** A row-by-row breakdown for every element in the grid, explicitly listing the row, column, black percentage (defect fraction), and grayscale value.

This Excel spreadsheet is now ready to be imported into your finite element analysis module to compute the defect-based material properties using the Mori-Tanaka scheme.
