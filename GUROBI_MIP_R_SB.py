from gurobipy import Model, GRB, quicksum
import matplotlib.pyplot as plt
import numpy as np
import time
import math
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as mcolors
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
if not os.path.exists('GUROBI_MIP_R_SB'):
    os.makedirs('GUROBI_MIP_R_SB')

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

def solve_strip_packing(widths, heights, strip_width, time_limit=1800, allow_rotation=True):
    """
    Solve the Strip Packing Problem using Gurobi MIP with time limit and rotation option.
    Returns the optimal height and positions of rectangles.
    """
    global best_height, best_positions, best_rotations
    
    # Create model
    mdl = Model("StripPacking")
    mdl.setParam('OutputFlag', 1)  # Enable output logs
    mdl.setParam('TimeLimit', time_limit)  # Set time limit in seconds
    mdl.setParam('MIPGap', 0.001)  # 0.1% relative gap
    mdl.setParam('MIPGapAbs', 0.01)  # 0.01 absolute gap
    
    # Number of rectangles
    n = len(widths)
    
    # Calculate bounds
    total_area = sum(widths[i] * heights[i] for i in range(n))
    
    if allow_rotation:
        # With rotation: lower bound uses max of minimum dimensions
        max_single_height = max(min(widths[i], heights[i]) for i in range(n))
        area_bound = total_area / strip_width
        lower_bound = max(max_single_height, area_bound)
        upper_bound_h = min(sum(max(widths[i], heights[i]) for i in range(n)), 
                           calculate_first_fit_upper_bound(strip_width, list(zip(widths, heights)), allow_rotation))
    else:
        max_single_height = max(heights)
        area_bound = total_area / strip_width
        lower_bound = max(max_single_height, area_bound)
        upper_bound_h = min(sum(heights), 
                           calculate_first_fit_upper_bound(strip_width, list(zip(widths, heights)), allow_rotation))
    
    # Variables
    x = mdl.addVars(n, lb=0, ub=strip_width, vtype=GRB.CONTINUOUS, name="x")
    y = mdl.addVars(n, lb=0, vtype=GRB.CONTINUOUS, name="y")
    H = mdl.addVar(lb=lower_bound, ub=upper_bound_h, vtype=GRB.CONTINUOUS, name="H")
    
    # Rotation variables
    if allow_rotation:
        rotate = mdl.addVars(n, vtype=GRB.BINARY, name="rotate")
    else:
        rotate = {i: 0 for i in range(n)}
    
    # Non-overlapping variables
    lr = mdl.addVars([(i,j) for i in range(n) for j in range(i+1, n)], vtype=GRB.BINARY, name="lr")
    rl = mdl.addVars([(i,j) for i in range(n) for j in range(i+1, n)], vtype=GRB.BINARY, name="rl")
    ab = mdl.addVars([(i,j) for i in range(n) for j in range(i+1, n)], vtype=GRB.BINARY, name="ab")
    ba = mdl.addVars([(i,j) for i in range(n) for j in range(i+1, n)], vtype=GRB.BINARY, name="ba")
    
    # Big-M values
    M_width = strip_width
    M_height = upper_bound_h
    
    # Effective width and height variables
    w = mdl.addVars(n, vtype=GRB.CONTINUOUS, name="w")
    h = mdl.addVars(n, vtype=GRB.CONTINUOUS, name="h")
    
    # Constraints
    for i in range(n):
        if allow_rotation:
            mdl.addConstr(w[i] == rotate[i] * heights[i] + (1 - rotate[i]) * widths[i], f"w_def_{i}")
            mdl.addConstr(h[i] == rotate[i] * widths[i] + (1 - rotate[i]) * heights[i], f"h_def_{i}")
        else:
            mdl.addConstr(w[i] == widths[i], f"w_def_{i}")
            mdl.addConstr(h[i] == heights[i], f"h_def_{i}")
        
        mdl.addConstr(x[i] + w[i] <= strip_width, f"width_bound_{i}")
        mdl.addConstr(y[i] + h[i] <= H, f"height_bound_{i}")
    
    # Non-overlapping constraints
    for i in range(n):
        for j in range(i + 1, n):
            mdl.addConstr(lr[i,j] + rl[i,j] + ab[i,j] + ba[i,j] >= 1, f"disjunction_{i}_{j}")
            mdl.addConstr(x[i] + w[i] <= x[j] + M_width * (1 - lr[i,j]), f"left_{i}_{j}")
            mdl.addConstr(x[j] + w[j] <= x[i] + M_width * (1 - rl[i,j]), f"right_{i}_{j}")
            mdl.addConstr(y[i] + h[i] <= y[j] + M_height * (1 - ab[i,j]), f"above_{i}_{j}")
            mdl.addConstr(y[j] + h[j] <= y[i] + M_height * (1 - ba[i,j]), f"below_{i}_{j}")
            
            # ===== SYMMETRY BREAKING CONSTRAINTS =====
            
            # 1. Large Rectangles Constraint (C2 - large)
            # Nếu hai hình chữ nhật quá lớn để đặt cạnh nhau, buộc chúng xếp chồng
            min_width_i = min(widths[i], heights[i])
            min_width_j = min(widths[j], heights[j])
            
            if min_width_i + min_width_j > strip_width:
                # Nếu không thể đặt cạnh nhau, buộc phải xếp chồng (vô hiệu hóa vị trí ngang)
                mdl.addConstr(lr[i,j] == 0, f"large_rect_h1_{i}_{j}")
                mdl.addConstr(rl[i,j] == 0, f"large_rect_h2_{i}_{j}")
            
            # 2. Large Rectangles Vertical Constraint
            min_height_i = min(widths[i], heights[i])
            min_height_j = min(widths[j], heights[j])
            
            # Sử dụng upper_bound thay cho height cố định
            # Nếu không thể xếp chồng, buộc phải đặt cạnh nhau
            if min_height_i + min_height_j > upper_bound_h:
                mdl.addConstr(ab[i,j] == 0, f"large_rect_v1_{i}_{j}")
                mdl.addConstr(ba[i,j] == 0, f"large_rect_v2_{i}_{j}")
                
            # 3. Same Rectangles Constraint (Mới thêm)
            # Nếu hai hình chữ nhật có cùng kích thước, cố định mối quan hệ vị trí
            if widths[i] == widths[j] and heights[i] == heights[j]:
                # Hình j không thể nằm bên trái hình i
                mdl.addConstr(rl[i,j] == 0, f"same_rect_left_{i}_{j}")
                
                # Thêm ràng buộc: hoặc i nằm bên trái j, hoặc j nằm dưới i
                mdl.addConstr(lr[i,j] + ba[i,j] >= 1, f"same_rect_pos_{i}_{j}")
    
    # 4. Domain Reduction Constraint (Mới thêm)
    # Tìm hình chữ nhật có diện tích lớn nhất
    max_area_idx = 0
    max_area = widths[0] * heights[0]
    for i in range(1, n):
        area = widths[i] * heights[i]
        if area > max_area:
            max_area = area
            max_area_idx = i
    
    # Giới hạn vị trí x của hình chữ nhật lớn nhất đến nửa đầu của strip
    # Chỉ áp dụng nếu không xoay
    max_width = widths[max_area_idx]
    if allow_rotation:
        # Khi cho phép xoay, chỉ áp dụng khi không xoay
        mdl.addConstr(x[max_area_idx] <= (strip_width - max_width)/2 + 
                     M_width * rotate[max_area_idx], f"domain_reduction_{max_area_idx}")
    else:
        mdl.addConstr(x[max_area_idx] <= (strip_width - max_width)/2, f"domain_reduction_{max_area_idx}")
    
    # 5. One Pair Constraint (C2 - pair) - Đã có sẵn
    if n >= 2:
        # Hình chữ nhật 1 không thể nằm bên trái hình chữ nhật 0
        mdl.addConstr(rl[0,1] == 0, "pair_left")
        # Hình chữ nhật 1 không thể nằm dưới hình chữ nhật 0
        mdl.addConstr(ba[0,1] == 0, "pair_below")
    
    # Objective: minimize height
    mdl.setObjective(H, GRB.MINIMIZE)
    
    # Save checkpoint before solving
    save_checkpoint(instance_id, upper_bound_h)
    
    # Optimize
    mdl.optimize()
    
    # Process results
    if mdl.status == GRB.OPTIMAL or mdl.status == GRB.TIME_LIMIT:
        if mdl.solCount > 0:
            current_height = H.X
            positions = [(x[i].X, y[i].X) for i in range(n)]
            
            if allow_rotation:
                rotations = [rotate[i].X for i in range(n)]
            else:
                rotations = [0] * n
                
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
                'objective_value': mdl.objVal,
                'solve_time': mdl.Runtime,
                'status': 'Optimal' if mdl.status == GRB.OPTIMAL else 'Time Limit',
                'lower_bound': lower_bound,
                'upper_bound': upper_bound_h
            }
    else:
        print(f"No feasible solution found. Status: {mdl.status}")
        return None

def display_solution(strip, rectangles, pos_circuits, rotations, instance_name):
    """Visualize the packing solution with improved styling and rotation"""
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.title(f"Strip Packing Solution for {instance_name} (Width: {strip[0]}, Height: {strip[1]:.2f})")

    if len(pos_circuits) > 0:
        # Use colormap to make visualization more pleasing
        colors = cm.viridis(np.linspace(0, 1, len(rectangles)))
        
        for i in range(len(rectangles)):
            # Apply rotation if needed
            w, h = rectangles[i]
            if rotations[i] >= 0.5:
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
            rot_info = 'R' if rotations[i] >= 0.5 else 'NR'
            ax.annotate(f"{i}\n{rot_info}", (cx, cy), color=text_color, 
                        ha='center', va='center', fontsize=10, fontweight='bold')

    # Use tight layout to maximize figure space
    ax.set_xlim(0, strip[0])
    ax.set_ylim(0, strip[1] + 1)
    
    # Set grid
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Add axis ticks with appropriate spacing
    max_ticks = 20
    if strip[0] <= max_ticks:
        ax.set_xticks(range(strip[0] + 1))
    else:
        ax.set_xticks(np.linspace(0, strip[0], max_ticks).astype(int))
    
    if int(strip[1]) <= max_ticks:
        ax.set_yticks(range(int(strip[1]) + 2))
    else:
        ax.set_yticks(np.linspace(0, int(strip[1]) + 1, max_ticks).astype(int))
    
    # Add labels
    ax.set_xlabel('Width', fontsize=12)
    ax.set_ylabel('Height', fontsize=12)
    
    # Save the plot to the output folder
    plt.savefig(f'GUROBI_MIP_R_SB/{instance_name}.png')
    plt.close()

if __name__ == "__main__":
    # Phần controller mode
    if len(sys.argv) == 1:
        # This is the controller mode - running without arguments
        # Create output folder if it doesn't exist
        if not os.path.exists('GUROBI_MIP_R_SB'):
            os.makedirs('GUROBI_MIP_R_SB')
        
        # Đọc file Excel hiện có để kiểm tra instances đã hoàn thành
        excel_file = 'GUROBI_MIP_R_SB.xlsx'
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
            command = f"./runlim --real-time-limit={TIMEOUT} python3 GUROBI_MIP_R_SB.py {instance_id}"
            
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

            print(f"Solving 2D Strip Packing with Gurobi MIP (with rotation) for instance {instance_name}")
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
            excel_file = 'GUROBI_MIP_R_SB.xlsx'
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
            excel_file = 'GUROBI_MIP_R_SB.xlsx'
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