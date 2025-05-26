from ortools.linear_solver import pywraplp
import matplotlib.pyplot as plt
import numpy as np
import time
import math
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import timeit
import sys
import signal
import os
import json
import subprocess
import traceback
import pandas as pd

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
        'Status': 'TIMEOUT',
        'Allow_Rotation': 'Yes'
    }
    
    with open(f'results_{instance_id}.json', 'w') as f:
        json.dump(result, f)
    
    sys.exit(0)

# Register signal handlers
signal.signal(signal.SIGTERM, handle_interrupt)  # Termination signal
signal.signal(signal.SIGINT, handle_interrupt)   # Keyboard interrupt (Ctrl+C)

# Create output folder if it doesn't exist
if not os.path.exists('OR-TOOLS_MIP_R_SB'):
    os.makedirs('OR-TOOLS_MIP_R_SB')

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
        'Status': status,
        'Allow_Rotation': 'Yes'
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

def solve_strip_packing(strip_width, items, time_limit=1800):
    """
    Solves the strip packing problem with rotation option and time limit.
    
    Args:
        strip_width: Width of the strip
        items: List of tuples (width, height) for each rectangle
        time_limit: Time limit in seconds
    
    Returns:
        Dictionary with solution or None if no solution found
    """
    global best_height, best_positions, best_rotations, upper_bound
    
    # Create the solver
    solver = pywraplp.Solver.CreateSolver('SCIP')
    if not solver:
        return None
    
    # Set time limit
    solver.set_time_limit(time_limit * 1000)  # Time limit in milliseconds
    
    n = len(items)
    
    # Calculate bounds
    total_area = sum(w * h for w, h in items)
    # For rotation, the max height is the minimum dimension of each rectangle
    max_height = max(min(w, h) for w, h in items)
    area_lb = math.ceil(total_area / strip_width)
    lower_bound = max(area_lb, max_height)
    
    # Upper bound with rotation - can use the smallest dimension in height
    ub = min(sum(max(w, h) for w, h in items), 
             calculate_first_fit_upper_bound(strip_width, items, allow_rotation=True))
    upper_bound = ub
    
    # Save checkpoint with initial upper bound
    save_checkpoint(instance_id, ub)
    
    # Variables
    x = [solver.NumVar(0, strip_width, f'x_{i}') for i in range(n)]
    y = [solver.NumVar(0, ub, f'y_{i}') for i in range(n)]
    rotate = [solver.BoolVar(f'rotate_{i}') for i in range(n)]  # 1 if rotated, 0 if not
    
    # Binary variables for relative positions
    left = [[solver.BoolVar(f'left_{i}_{j}') for j in range(n)] 
            for i in range(n)]
    below = [[solver.BoolVar(f'below_{i}_{j}') for j in range(n)] 
             for i in range(n)]
    
    # Height to minimize
    height = solver.NumVar(lower_bound, ub, 'height')

    # Constraints
    for i in range(n):
        w, h = items[i]
        # Width constraint considering rotation
        solver.Add(x[i] + w * (1 - rotate[i]) + h * rotate[i] <= strip_width)
        # Height constraint considering rotation
        solver.Add(y[i] + h * (1 - rotate[i]) + w * rotate[i] <= height)
        
        for j in range(i + 1, n):
            wj, hj = items[j]
            M = strip_width + ub  # Big-M value
            
            # Left constraint considering rotations
            solver.Add(x[i] + 
                      w * (1 - rotate[i]) + h * rotate[i] <= 
                      x[j] + M * (1 - left[i][j]))
            solver.Add(x[j] + 
                      wj * (1 - rotate[j]) + hj * rotate[j] <= 
                      x[i] + M * (1 - left[j][i]))
            
            # Below constraint considering rotations
            solver.Add(y[i] + 
                      h * (1 - rotate[i]) + w * rotate[i] <= 
                      y[j] + M * (1 - below[i][j]))
            solver.Add(y[j] + 
                      hj * (1 - rotate[j]) + wj * rotate[j] <= 
                      y[i] + M * (1 - below[j][i]))
            
            # At least one condition must hold
            solver.Add(left[i][j] + left[j][i] + below[i][j] + below[j][i] >= 1)
            
            # ===== SYMMETRY BREAKING CONSTRAINTS =====
            
            # 1. Large Rectangles Constraint (C2 - large)
            # If two rectangles can't fit side by side even with rotation
            min_width_i = min(w, h)
            min_width_j = min(wj, hj)
            
            if min_width_i + min_width_j > strip_width:
                # They must be stacked vertically
                solver.Add(left[i][j] == 0)
                solver.Add(left[j][i] == 0)
            
            # 2. Large Rectangles Vertical Constraint
            min_height_i = min(w, h)
            min_height_j = min(wj, hj)
            
            # If they can't be stacked vertically
            if min_height_i + min_height_j > ub:
                # They must be placed side by side
                solver.Add(below[i][j] == 0)
                solver.Add(below[j][i] == 0)
            
            # 3. Same Rectangles Symmetry Breaking
            if w == wj and h == hj:  # Same dimensions
                # Force a fixed positional relationship to break symmetry
                solver.Add(left[j][i] == 0)  # Rectangle j cannot be to the left of rectangle i
                # Additional constraint: either i is left of j or j is below i
                solver.Add(left[i][j] + below[j][i] >= 1)
    
    # 4. One Pair Constraint (C2 - pair)
    if n >= 2:
        # Rectangle 1 cannot be to the left of rectangle 0
        solver.Add(left[1][0] == 0)
        # Rectangle 1 cannot be below rectangle 0
        solver.Add(below[1][0] == 0)
    
    # 5. Find max area rectangle for domain reduction
    max_area_idx = 0
    max_area = items[0][0] * items[0][1]
    for i in range(1, n):
        area = items[i][0] * items[i][1]
        if area > max_area:
            max_area = area
            max_area_idx = i
    
    # Domain Reduction constraint (C2 - domain)
    # Limit x position of largest rectangle to first half of strip
    w_max, h_max = items[max_area_idx]
    # Considering possible rotation for domain reduction
    solver.Add(x[max_area_idx] <= (strip_width - min(w_max, h_max)) / 2)

    # Objective: minimize height
    solver.Minimize(height)

    # Save checkpoint before solving
    save_checkpoint(instance_id, best_height if best_height != float('inf') else upper_bound)
    
    # Solve
    status = solver.Solve()
    
    if status == pywraplp.Solver.OPTIMAL or status == pywraplp.Solver.FEASIBLE:
        solution_height = height.solution_value()
        positions = [(x[i].solution_value(), y[i].solution_value()) for i in range(n)]
        rotations = [rotate[i].solution_value() for i in range(n)]
        
        # Update best solution
        if solution_height < best_height:
            best_height = solution_height
            best_positions = positions
            best_rotations = rotations
            
            # Save checkpoint after finding better solution
            save_checkpoint(instance_id, best_height)
        
        solution = {
            'height': solution_height,
            'positions': positions,
            'rotations': rotations,
            'status': 'Optimal' if status == pywraplp.Solver.OPTIMAL else 'Feasible',
            'objective_value': solver.Objective().Value(),
            'solve_time': solver.WallTime() / 1000.0  # Convert to seconds
        }
        return solution
    else:
        print("No optimal solution found")
        return None

def display_solution(strip_width, items, solution, instance_name):
    """Visualize the packing solution with rotations and save to file"""
    height = solution['height']
    positions = solution['positions']
    rotations = solution['rotations']
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.title(f"Strip Packing Solution for {instance_name} (Width: {strip_width}, Height: {height:.2f})")
    
    colors = cm.viridis(np.linspace(0, 1, len(items)))
    
    for i, (pos_x, pos_y) in enumerate(positions):
        w, h = items[i]
        if rotations[i]:
            w, h = h, w  # Apply rotation
            
        rect = plt.Rectangle((pos_x, pos_y), w, h, 
                           edgecolor='black', facecolor=colors[i], alpha=0.7)
        ax.add_patch(rect)
        
        # Add rectangle number at center with rotation info
        cx, cy = pos_x + w/2, pos_y + h/2
        
        # Calculate text color based on background brightness
        rgb = mcolors.to_rgb(colors[i])
        brightness = 0.299*rgb[0] + 0.587*rgb[1] + 0.114*rgb[2]
        text_color = 'white' if brightness < 0.6 else 'black'
        
        rot_info = 'R' if rotations[i] else 'NR'
        ax.text(cx, cy, f'{i}\n{rot_info}', ha='center', va='center', 
                color=text_color, fontweight='bold')
    
    ax.set_xlim(0, strip_width)
    ax.set_ylim(0, height)
    ax.set_xticks(range(0, strip_width+1, max(1, strip_width//20)))
    ax.set_yticks(range(0, int(height)+1, max(1, int(height)//20)))
    ax.set_xlabel('Width')
    ax.set_ylabel('Height')
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Save the plot to output folder
    plt.savefig(f'OR-TOOLS_MIP_R_SB/{instance_name}.png')
    plt.close()

if __name__ == "__main__":
    # Phần controller mode
    if len(sys.argv) == 1:
        # This is the controller mode - running without arguments
        # Create output folder if it doesn't exist
        if not os.path.exists('OR-TOOLS_MIP_R_SB'):
            os.makedirs('OR-TOOLS_MIP_R_SB')
        
        # Đọc file Excel hiện có để kiểm tra instances đã hoàn thành
        excel_file = 'OR-TOOLS_MIP_R_SB.xlsx'
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
        TIMEOUT = 1800  # 30 minutes timeout
        
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
            command = f"./runlim --real-time-limit={TIMEOUT} python3 OR-TOOLS_MIP_R_SB.py {instance_id}"
            
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
            items = []
            
            for i in range(2, 2+n_rec):
                if i < len(input_data):
                    rect = [int(val) for val in input_data[i].split()]
                    items.append(rect)
                else:
                    raise IndexError(f"Missing rectangle data at line {i}")
            
            # Calculate initial bounds
            total_area = sum(w * h for w, h in items)
            max_height = max(min(h, w) for w, h in items)  # Consider rotation
            area_lb = math.ceil(total_area / strip_width)
            lower_bound = max(area_lb, max_height)
            
            # Upper bound considering rotation
            upper_bound = min(sum(max(h, w) for w, h in items), 
                            calculate_first_fit_upper_bound(strip_width, items, allow_rotation=True))

            print(f"Solving 2D Strip Packing with OR-Tools MIP (with rotation) for instance {instance_name}")
            print(f"Width: {strip_width}")
            print(f"Number of rectangles: {n_rec}")
            print(f"Lower bound: {lower_bound}")
            print(f"Upper bound: {upper_bound}")
            
            # Solve with MIP
            solution = solve_strip_packing(strip_width, items, time_limit=1800) # Leave 60 seconds for other operations
            
            stop = timeit.default_timer()
            runtime = stop - start

            # Display and save the solution if we found one
            if solution:
                display_solution(strip_width, items, solution, instance_name)
                solution_status = solution['status']
                solution_height = solution['height']
            else:
                print("No feasible solution found.")
                solution_status = 'ERROR'
                solution_height = best_height if best_height != float('inf') else upper_bound

            # Tạo result object
            result_data = {
                'Instance': instance_name,
                'Runtime': runtime,
                'Optimal_Height': solution_height,
                'Status': 'COMPLETE' if solution and solution_status == 'Optimal' else 'SUBOPTIMAL',
                'Allow_Rotation': 'Yes'
            }
            
            # Ghi kết quả vào Excel trực tiếp
            excel_file = 'OR-TOOLS_MIP_R_SB.xlsx'
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
            
            print(f"Instance {instance_name} completed - Runtime: {runtime:.2f}s, Height: {solution_height}")

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
            excel_file = 'OR-TOOLS_MIP_R_SB.xlsx'
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