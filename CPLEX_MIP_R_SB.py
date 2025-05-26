from docplex.mp.model import Model
import matplotlib.pyplot as plt
import pandas as pd
import time
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import math
import timeit
import sys
import signal
import os
import json
import subprocess
import traceback

# Global variables to track best solution found so far
best_height = float('inf')
best_positions = []
best_rotations = []
solve_time = 0
upper_bound = 0

# Signal handler for graceful interruption (e.g., by runlim)
def handle_interrupt(signum, frame):
    print(f"\nReceived interrupt signal {signum}. Saving current best solution.")
    
    # Lấy chiều cao tốt nhất (hoặc là giá trị tìm được, hoặc là upper_bound)
    current_height = best_height if best_height != float('inf') else upper_bound
    print(f"Best height found before interrupt: {current_height}")
    
    # Save result as JSON for the controller to pick up
    result = {
        'Instance': instances[instance_id],
        'Runtime': timeit.default_timer() - start,
        'Optimal_Height': current_height,
        'Status': 'TIMEOUT'
    }
    
    with open(f'results_{instance_id}.json', 'w') as f:
        json.dump(result, f)
    
    sys.exit(0)

# Register signal handlers
signal.signal(signal.SIGTERM, handle_interrupt)  # Termination signal
signal.signal(signal.SIGINT, handle_interrupt)   # Keyboard interrupt (Ctrl+C)

# Create output folder if it doesn't exist
if not os.path.exists('CPLEX_MIP_R_SB'):
    os.makedirs('CPLEX_MIP_R_SB')

def read_file_instance(n_instance):
    """Read problem instance from file."""
    s = ''
    filepath = f"inputs/ins-{n_instance}.txt"
    try:
        with open(filepath, 'r') as file:
            s = file.read()
        return s.splitlines()
    except FileNotFoundError:
        print(f"Error: File {filepath} not found")
        return None

instances = [ "",
  "HT01(c1p1)", "HT02(c1p2)", "HT03(c1p3)", "HT04(c2p1)", "HT05(c2p2)", "HT06(c2p3)", 
  "HT07(c3p1)", "HT08(c3p2)", "HT09(c3p3)", 
  "CGCUT01", "CGCUT02", "CGCUT03", 
  "GCUT01", "GCUT02", "GCUT03", "GCUT04", 
  "NGCUT01", "NGCUT02", "NGCUT03", "NGCUT04", "NGCUT05", "NGCUT06", "NGCUT07", 
  "NGCUT08", "NGCUT09", "NGCUT10", "NGCUT11", "NGCUT12", 
  "BENG01", "BENG02", "BENG03", "BENG04", "BENG05", "BENG06", "BENG07", "BENG08", "BENG09", "BENG10", "HT10(c4p1)", "HT11(c4p2)", "HT12(c4p3)"
]

# Thêm hàm save_checkpoint để lưu tiến trình giải
def save_checkpoint(instance_id, height, status="IN_PROGRESS"):
    checkpoint = {
        'Runtime': timeit.default_timer() - start,
        'Optimal_Height': height if height != float('inf') else upper_bound,
        'Status': status
    }
    
    # Ghi ra file checkpoint
    with open(f'checkpoint_{instance_id}.json', 'w') as f:
        json.dump(checkpoint, f)

def calculate_first_fit_upper_bound(width, rectangles, allow_rotation=True):
    # Sort rectangles by height (non-increasing)
    if allow_rotation:
        # Consider both orientations and choose the one with smaller height
        sorted_rects = sorted([(min(w, h), max(w, h)) for w, h in rectangles], key=lambda r: -r[1])
    else:
        sorted_rects = sorted(rectangles, key=lambda r: -r[1])
    
    levels = []  # (y-position, remaining_width)
    
    for w, h in sorted_rects:
        # Try to place on existing level
        placed = False
        for i in range(len(levels)):
            if levels[i][1] >= w:
                levels[i] = (levels[i][0], levels[i][1] - w)
                placed = True
                break
        
        # Create new level if needed
        if not placed:
            y_pos = 0 if not levels else max(level[0] + sorted_rects[i][1] for i, level in enumerate(levels))
            levels.append((y_pos, width - w))
    
    # Calculate total height
    if not levels:
        return 0
        
    return max(level[0] + sorted_rects[levels.index(level)][1] for level in levels)

def solve_strip_packing(widths, heights, strip_width, time_limit=1800, allow_rotation=True):
    """
    Solve the Strip Packing Problem using CPLEX MIP with time limit and rotation option.
    Returns the optimal height and positions of rectangles.
    """
    global best_height, best_positions, best_rotations
    
    # Create model
    mdl = Model("StripPacking")
    
    # Number of rectangles
    n = len(widths)
    
    # Calculate bounds
    total_area = sum(widths[i] * heights[i] for i in range(n))
    if allow_rotation:
        max_height = max(min(widths[i], heights[i]) for i in range(n))
        max_width = max(min(widths[i], heights[i]) for i in range(n))
        area_lb = total_area / strip_width
        lower_bound = max(max_height, area_lb)
        upper_bound_h = min(sum(max(widths[i], heights[i]) for i in range(n)), 
                           calculate_first_fit_upper_bound(strip_width, list(zip(widths, heights)), allow_rotation))
    else:
        max_height = max(heights)
        max_width = max(widths)
        area_lb = total_area / strip_width
        lower_bound = max(max_height, area_lb)
        upper_bound_h = min(sum(heights), calculate_first_fit_upper_bound(strip_width, list(zip(widths, heights)), allow_rotation))
    
    # Variables
    x = mdl.continuous_var_list(n, lb=0, ub=strip_width, name="x")
    y = mdl.continuous_var_list(n, lb=0, name="y")
    
    # Maximum height (objective to minimize)
    H = mdl.continuous_var(lb=lower_bound, ub=upper_bound_h, name="H")
    
    # Binary variables for non-overlapping constraints
    lr = mdl.binary_var_matrix(range(n), range(n), name="lr")  # i left of j
    rl = mdl.binary_var_matrix(range(n), range(n), name="rl")  # j left of i
    ab = mdl.binary_var_matrix(range(n), range(n), name="ab")  # i above j
    ba = mdl.binary_var_matrix(range(n), range(n), name="ba")  # j above i
    
    # Rotation variables (1 if rotated, 0 if not)
    rotate = []
    if allow_rotation:
        rotate = mdl.binary_var_list(n, name="rotate")
        
        # Width and height after rotation
        w_eff = mdl.continuous_var_list(n, name="w_eff")  # Effective width
        h_eff = mdl.continuous_var_list(n, name="h_eff")  # Effective height
        
        # Define effective width and height based on rotation
        for i in range(n):
            mdl.add_constraint(w_eff[i] == widths[i] * (1 - rotate[i]) + heights[i] * rotate[i], f"w_eff_{i}")
            mdl.add_constraint(h_eff[i] == heights[i] * (1 - rotate[i]) + widths[i] * rotate[i], f"h_eff_{i}")
    else:
        # No rotation, use original dimensions
        w_eff = widths
        h_eff = heights
    
    # Big M constant
    M = strip_width + upper_bound_h
    
    # Domain constraints
    for i in range(n):
        # Stay within strip boundaries
        mdl.add_constraint(x[i] + w_eff[i] <= strip_width, f"width_bound_{i}")
        mdl.add_constraint(y[i] + h_eff[i] <= H, f"height_bound_{i}")
    
    # Non-overlapping constraints
    for i in range(n):
        for j in range(i + 1, n):
            # At least one of these conditions must be true
            mdl.add_constraint(lr[i,j] + rl[i,j] + ab[i,j] + ba[i,j] >= 1, f"no_overlap_{i}_{j}")
            
            # Define relationships
            mdl.add_constraint(x[i] + w_eff[i] <= x[j] + M * (1 - lr[i,j]), f"left_{i}_{j}")
            mdl.add_constraint(x[j] + w_eff[j] <= x[i] + M * (1 - rl[i,j]), f"right_{i}_{j}")
            mdl.add_constraint(y[i] + h_eff[i] <= y[j] + M * (1 - ab[i,j]), f"above_{i}_{j}")
            mdl.add_constraint(y[j] + h_eff[j] <= y[i] + M * (1 - ba[i,j]), f"below_{i}_{j}")
    
    # SYMMETRY BREAKING CONSTRAINTS
    
    # 1. Large Rectangles Constraint (Horizontal) - C2 configuration
    for i in range(n):
        for j in range(i+1, n):
            if allow_rotation:
                # With rotation: check minimum possible width of each rectangle
                min_width_i = min(widths[i], heights[i])
                min_width_j = min(widths[j], heights[j])
                
                # If rectangles can't fit side by side even with rotation
                if min_width_i + min_width_j > strip_width:
                    # Disable horizontal placements
                    mdl.add_constraint(lr[i,j] == 0, f"large_rect_h1_{i}_{j}")
                    mdl.add_constraint(rl[i,j] == 0, f"large_rect_h2_{i}_{j}")
            else:
                # Without rotation: check actual widths
                if widths[i] + widths[j] > strip_width:
                    # Disable horizontal placements
                    mdl.add_constraint(lr[i,j] == 0, f"large_rect_h1_{i}_{j}")
                    mdl.add_constraint(rl[i,j] == 0, f"large_rect_h2_{i}_{j}")
    
    # 2. Large Rectangles Constraint (Vertical) - NEW
    for i in range(n):
        for j in range(i+1, n):
            if allow_rotation:
                # With rotation: check minimum possible height of each rectangle
                min_height_i = min(widths[i], heights[i])
                min_height_j = min(widths[j], heights[j])
                
                # If rectangles can't be stacked vertically even with rotation
                if min_height_i + min_height_j > upper_bound_h:
                    # Disable vertical placements
                    mdl.add_constraint(ab[i,j] == 0, f"large_rect_v1_{i}_{j}")
                    mdl.add_constraint(ba[i,j] == 0, f"large_rect_v2_{i}_{j}")
            else:
                # Without rotation: check actual heights
                if heights[i] + heights[j] > upper_bound_h:
                    # Disable vertical placements
                    mdl.add_constraint(ab[i,j] == 0, f"large_rect_v1_{i}_{j}")
                    mdl.add_constraint(ba[i,j] == 0, f"large_rect_v2_{i}_{j}")
    
    # 3. Same Rectangles Constraint - NEW
    for i in range(n):
        for j in range(i+1, n):
            if widths[i] == widths[j] and heights[i] == heights[j]:
                # For identical rectangles, fix their relative positions to break symmetry
                # Rectangle j can't be to the left of rectangle i
                mdl.add_constraint(rl[i,j] == 0, f"same_rect_lr_{i}_{j}")
                
                # Either i is to the left of j, or j is below i
                mdl.add_constraint(lr[i,j] + ba[i,j] >= 1, f"same_rect_pos_{i}_{j}")
    
    # 4. Domain Reduction (Maximum Rectangle) - NEW
    if n > 0:
        # Find the rectangle with the maximum area
        max_area_idx = 0
        max_area = widths[0] * heights[0]
        for i in range(1, n):
            area = widths[i] * heights[i]
            if area > max_area:
                max_area = area
                max_area_idx = i
        
        # Limit the x-position of the maximum rectangle to the first half of the strip
        if allow_rotation:
            # Only apply if the rectangle is not rotated
            mdl.add_constraint(x[max_area_idx] <= (strip_width - widths[max_area_idx]) / 2 + 
                              strip_width * rotate[max_area_idx], f"domain_reduction_{max_area_idx}")
        else:
            mdl.add_constraint(x[max_area_idx] <= (strip_width - widths[max_area_idx]) / 2, 
                             f"domain_reduction_{max_area_idx}")
    
    # 5. One Pair of Rectangles Constraint (already implemented)
    if n >= 2:
        # Rectangle 1 cannot be to the left of rectangle 0
        mdl.add_constraint(rl[0,1] == 0, "pair_constraint1")
        # Rectangle 1 cannot be below rectangle 0 
        mdl.add_constraint(ba[0,1] == 0, "pair_constraint2")
    
    # Objective: minimize height
    mdl.minimize(H)
    
    # Set time limit if provided
    if time_limit:
        mdl.set_time_limit(time_limit)
    
    # Save checkpoint before solving
    save_checkpoint(instance_id, upper_bound_h)
    
    # Solve
    solution = mdl.solve(log_output=True)
    
    if solution:
        # Process results including rotation information
        if allow_rotation:
            rotations = [solution.get_value(rotate[i]) for i in range(n)]
            actual_widths = [heights[i] if rotations[i] > 0.5 else widths[i] for i in range(n)]
            actual_heights = [widths[i] if rotations[i] > 0.5 else heights[i] for i in range(n)]
        else:
            rotations = [0] * n
            actual_widths = widths
            actual_heights = heights
        
        current_height = solution[H]
        positions = [(solution[x[i]], solution[y[i]]) for i in range(n)]
        
        # Update global best solution
        if current_height < best_height:
            best_height = current_height
            best_positions = positions
            best_rotations = rotations
            
            # Save checkpoint after finding better solution
            save_checkpoint(instance_id, best_height)
        
        return {
            'height': current_height,
            'positions': positions,
            'rotations': rotations,
            'actual_widths': actual_widths,
            'actual_heights': actual_heights,
            'objective_value': mdl.objective_value
        }
    else:
        return None

def display_solution(strip, rectangles, pos_circuits, rotations, instance_name):
    """Visualize the packing solution with improved styling and rotation"""
    strip_width, strip_height = strip
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.title(f"Strip Packing Solution for {instance_name} (Width: {strip_width}, Height: {strip_height:.2f})")

    if len(pos_circuits) > 0:
        # Use colormap to make visualization more pleasing
        colors = cm.viridis(np.linspace(0, 1, len(rectangles)))
        
        for i in range(len(rectangles)):
            w, h = rectangles[i]
            # Apply rotation if needed
            if rotations[i] > 0.5:
                w, h = h, w
            
            # Create rectangle with nice coloring
            rect = plt.Rectangle(pos_circuits[i], w, h, 
                               edgecolor="#333", facecolor=colors[i], alpha=0.7)
            ax.add_patch(rect)
            
            # Add rectangle number at center
            rx, ry = pos_circuits[i]
            cx, cy = rx + w/2, ry + h/2
            
            # Calculate text color based on background brightness
            rgb = mcolors.to_rgb(colors[i])
            brightness = 0.299*rgb[0] + 0.587*rgb[1] + 0.114*rgb[2]
            text_color = 'white' if brightness < 0.6 else 'black'
            
            # Add rotation info
            rot_info = 'R' if rotations[i] > 0.5 else 'NR'
            ax.annotate(f"{i}\n{rot_info}", (cx, cy), color=text_color, 
                        ha='center', va='center', fontsize=10, fontweight='bold')

    # Use tight layout to maximize figure space
    ax.set_xlim(0, strip_width)
    ax.set_ylim(0, strip_height + 1)
    
    # Set grid
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Add axis ticks with appropriate spacing
    max_ticks = 20
    if strip_width <= max_ticks:
        ax.set_xticks(range(strip_width + 1))
    else:
        ax.set_xticks(np.linspace(0, strip_width, max_ticks).astype(int))
    
    if int(strip_height) <= max_ticks:
        ax.set_yticks(range(int(strip_height) + 2))
    else:
        ax.set_yticks(np.linspace(0, int(strip_height) + 1, max_ticks).astype(int))
    
    # Add labels
    ax.set_xlabel('Width', fontsize=12)
    ax.set_ylabel('Height', fontsize=12)
    
    # Save the plot to the output folder
    plt.savefig(f'CPLEX_MIP_R_SB/{instance_name}.png')
    plt.close()

if __name__ == "__main__":
    # Phần controller mode
    if len(sys.argv) == 1:
        # This is the controller mode - running without arguments
        # Create output folder if it doesn't exist
        if not os.path.exists('CPLEX_MIP_R_SB'):
            os.makedirs('CPLEX_MIP_R_SB')
        
        # Đọc file Excel hiện có để kiểm tra instances đã hoàn thành
        excel_file = 'CPLEX_MIP_R_SB.xlsx'
        if os.path.exists(excel_file):
            # Đọc file Excel hiện có nếu nó tồn tại
            existing_df = pd.read_excel(excel_file)
            # Lấy danh sách các instance đã hoàn thành
            completed_instances = existing_df['Instance'].tolist() if 'Instance' in existing_df.columns else []
        else:
            # Tạo DataFrame trống nếu chưa có file
            existing_df = pd.DataFrame()
            completed_instances = []
        
        # Set timeout in seconds
        TIMEOUT = 1200  # 30 minutes timeout
        
        for instance_id in range(1, 42):
            instance_name = instances[instance_id]
            
            # Kiểm tra xem instance này đã được chạy chưa
            if instance_name in completed_instances:
                print(f"\nSkipping instance {instance_id}: {instance_name} (already completed)")
                continue
                
            print(f"\n{'=' * 50}")
            print(f"Running instance {instance_id}: {instance_name}")
            print(f"{'=' * 50}")
            
            # Clean up any previous result file
            if os.path.exists(f'results_{instance_id}.json'):
                os.remove(f'results_{instance_id}.json')
            
            # Run the instance with runlim, but use THIS script with the instance_id
            command = f"./runlim --real-time-limit={TIMEOUT} python3 CPLEX_MIP_R_SB.py {instance_id}"
            
            try:
                # Run the command and wait for it to complete
                process = subprocess.Popen(command, shell=True)
                process.wait()
                
                # Wait a moment to ensure file is written
                time.sleep(1)
                
                # Kiểm tra kết quả
                result = None
                
                # Thử đọc file results trước (kết quả hoàn chỉnh)
                if os.path.exists(f'results_{instance_id}.json'):
                    with open(f'results_{instance_id}.json', 'r') as f:
                        result = json.load(f)
                
                # Nếu không tìm thấy file results, kiểm tra file checkpoint
                elif os.path.exists(f'checkpoint_{instance_id}.json'):
                    with open(f'checkpoint_{instance_id}.json', 'r') as f:
                        result = json.load(f)
                    # Đánh dấu đây là kết quả timeout
                    result['Status'] = 'TIMEOUT'
                    result['Instance'] = instance_name
                    print(f"Instance {instance_name} timed out. Using checkpoint data.")
                
                # Xử lý kết quả (nếu có)
                if result:
                    print(f"Instance {instance_name} - Status: {result['Status']}")
                    print(f"Optimal Height: {result['Optimal_Height']}, Runtime: {result['Runtime']}")
                    
                    # Cập nhật Excel
                    if os.path.exists(excel_file):
                        try:
                            existing_df = pd.read_excel(excel_file)
                            instance_exists = instance_name in existing_df['Instance'].tolist() if 'Instance' in existing_df.columns else False
                            
                            if instance_exists:
                                # Cập nhật instance đã tồn tại
                                instance_idx = existing_df.index[existing_df['Instance'] == instance_name].tolist()[0]
                                for key, value in result.items():
                                    existing_df.at[instance_idx, key] = value
                            else:
                                # Thêm instance mới
                                result_df = pd.DataFrame([result])
                                existing_df = pd.concat([existing_df, result_df], ignore_index=True)
                        except Exception as e:
                            print(f"Lỗi khi đọc file Excel hiện có: {str(e)}")
                            existing_df = pd.DataFrame([result])
                    else:
                        # Tạo DataFrame mới nếu chưa có file Excel
                        existing_df = pd.DataFrame([result])
                    # Lưu DataFrame vào Excel
                    existing_df.to_excel(excel_file, index=False)
                    print(f"Results saved to {excel_file}")
                        
                else:
                    print(f"No results or checkpoint found for instance {instance_name}")
                    
            except Exception as e:
                print(f"Error running instance {instance_name}: {str(e)}")
            
            # Clean up the results file to avoid confusion
            for file in [f'results_{instance_id}.json', f'checkpoint_{instance_id}.json']:
                if os.path.exists(file):
                    os.remove(file)
        
        print(f"\nAll instances completed. Results saved to {excel_file}")
    
    # Phần single instance mode
    else:
        # Single instance mode
        instance_id = int(sys.argv[1])
        instance_name = instances[instance_id]
        
        start = timeit.default_timer()  # start clock
        
        try:
            print(f"\nProcessing instance {instance_name}")
            
            # Reset global best solution for this instance
            best_height = float('inf')
            best_positions = []
            best_rotations = []

            # read file input
            input_data = read_file_instance(instance_id)
            
            if input_data is None:
                print(f"Failed to read instance {instance_id}")
                sys.exit(1)
                
            strip_width = int(input_data[0])
            n_rec = int(input_data[1])
            rectangles = []
            
            for i in range(2, 2+n_rec):
                if i < len(input_data):
                    rect = [int(val) for val in input_data[i].split()]
                    rectangles.append(rect)
                else:
                    raise IndexError(f"Missing rectangle data at line {i}")
            
            widths = [rect[0] for rect in rectangles]
            heights = [rect[1] for rect in rectangles]
            
            # Calculate initial bounds
            total_area = sum(w * h for w, h in rectangles)
            allow_rotation = True  # Enable rotation for this solver
            
            if allow_rotation:
                # With rotation: lower bound uses max of minimum dimensions
                max_dimension = max(min(w, h) for w, h in rectangles)
                area_lb = math.ceil(total_area / strip_width)
                lower_bound = max(max_dimension, area_lb)
                
                # Upper bound: either sum of all heights or first-fit
                upper_bound = min(sum(max(w, h) for w, h in rectangles), 
                                calculate_first_fit_upper_bound(strip_width, rectangles, allow_rotation))
            else:
                # Without rotation: standard bounds
                max_height = max(h for _, h in rectangles)
                area_lb = math.ceil(total_area / strip_width)
                lower_bound = max(max_height, area_lb)
                
                # Upper bound: either sum of all heights or first-fit
                upper_bound = min(sum(heights), 
                                calculate_first_fit_upper_bound(strip_width, rectangles, allow_rotation))

            print(f"Solving 2D Strip Packing with CPLEX MIP (with rotation) for instance {instance_name}")
            print(f"Width: {strip_width}")
            print(f"Number of rectangles: {n_rec}")
            print(f"Lower bound: {lower_bound}")
            print(f"Upper bound: {upper_bound}")
            
            # Solve with MIP
            result = solve_strip_packing(widths, heights, strip_width, time_limit=1800, allow_rotation=True)
            
            stop = timeit.default_timer()
            runtime = stop - start

            # Display and save the solution if we found one
            if result:
                optimal_height = result['height']
                display_solution((strip_width, optimal_height), rectangles, 
                               result['positions'], result['rotations'], instance_name)
            else:
                optimal_height = best_height if best_height != float('inf') else upper_bound
                print(f"No solution found, using best height: {optimal_height}")
                if best_positions and best_rotations:
                    display_solution((strip_width, optimal_height), rectangles, 
                                   best_positions, best_rotations, instance_name)

            # Tạo result object
            result_data = {
                'Instance': instance_name,
                'Runtime': runtime,
                'Optimal_Height': optimal_height,
                'Status': 'COMPLETE' if result else 'ERROR',
                'Allow_Rotation': 'Yes'
            }
            
            # Ghi kết quả vào Excel trực tiếp
            excel_file = 'CPLEX_MIP_R_SB.xlsx'
            if os.path.exists(excel_file):
                try:
                    existing_df = pd.read_excel(excel_file)
                    instance_exists = instance_name in existing_df['Instance'].tolist() if 'Instance' in existing_df.columns else False
                    
                    if instance_exists:
                        # Cập nhật instance đã tồn tại
                        instance_idx = existing_df.index[existing_df['Instance'] == instance_name].tolist()[0]
                        for key, value in result_data.items():
                            existing_df.at[instance_idx, key] = value
                    else:
                        # Thêm instance mới
                        result_df = pd.DataFrame([result_data])
                        existing_df = pd.concat([existing_df, result_df], ignore_index=True)
                except Exception as e:
                    print(f"Lỗi khi đọc file Excel hiện có: {str(e)}")
                    existing_df = pd.DataFrame([result_data])
            else:
                # Tạo DataFrame mới nếu chưa có file Excel
                existing_df = pd.DataFrame([result_data])
            
            # Lưu DataFrame vào Excel
            existing_df.to_excel(excel_file, index=False)
            print(f"Results saved to {excel_file}")
            
            # Save result to a JSON file that the controller will read
            with open(f'results_{instance_id}.json', 'w') as f:
                json.dump(result_data, f)
            
            print(f"Instance {instance_name} completed - Runtime: {runtime:.2f}s, Height: {optimal_height}")

        except Exception as e:
            print(f"Error in instance {instance_name}: {str(e)}")
            traceback.print_exc()  # Print the traceback for the error
            
            # Save error result - use upper_bound if no best_height
            current_height = best_height if best_height != float('inf') else upper_bound
            result_data = {
                'Instance': instance_name,
                'Runtime': timeit.default_timer() - start,
                'Optimal_Height': current_height,
                'Status': 'ERROR',
                'Allow_Rotation': 'Yes'
            }
            
            # Ghi kết quả lỗi vào Excel
            excel_file = 'CPLEX_MIP_R_SB.xlsx'
            if os.path.exists(excel_file):
                try:
                    existing_df = pd.read_excel(excel_file)
                    instance_exists = instance_name in existing_df['Instance'].tolist() if 'Instance' in existing_df.columns else False
                    
                    if instance_exists:
                        # Cập nhật instance đã tồn tại
                        instance_idx = existing_df.index[existing_df['Instance'] == instance_name].tolist()[0]
                        for key, value in result_data.items():
                            existing_df.at[instance_idx, key] = value
                    else:
                        # Thêm instance mới
                        result_df = pd.DataFrame([result_data])
                        existing_df = pd.concat([existing_df, result_df], ignore_index=True)
                except Exception as ex:
                    print(f"Lỗi khi đọc file Excel hiện có: {str(ex)}")
                    existing_df = pd.DataFrame([result_data])
            else:
                # Tạo DataFrame mới nếu chưa có file Excel
                existing_df = pd.DataFrame([result_data])
            
            # Lưu DataFrame vào Excel
            existing_df.to_excel(excel_file, index=False)
            print(f"Error results saved to {excel_file}")
            
            # Save result to a JSON file that the controller will read
            with open(f'results_{instance_id}.json', 'w') as f:
                json.dump(result_data, f)