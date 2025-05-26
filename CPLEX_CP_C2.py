from docplex.cp.model import CpoModel
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

# Create SPP folder if it doesn't exist
if not os.path.exists('CPLEX_CP_C2'):
    os.makedirs('CPLEX_CP_C2')

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

def solve_2OPP(strip_width, items, fixed_height, time_limit=1800):
    """
    Solves the 2D Orthogonal Packing Problem with fixed height using CP.
    
    Args:
        strip_width: Width of the strip
        items: List of (width, height) tuples for each rectangle
        fixed_height: Fixed height of the strip
        time_limit: Time limit in seconds
        
    Returns:
        None if no solution, or dict with positions if feasible
    """
    # Create the CP model
    model = CpoModel(name="2OPP")
    
    n = len(items)
    
    # Variables: coordinates of the bottom-left corner of each rectangle
    x = [model.integer_var(0, strip_width - items[i][0], f'x_{i}') for i in range(n)]
    y = [model.integer_var(0, fixed_height - items[i][1], f'y_{i}') for i in range(n)]
    
    # Create interval variables for each rectangle on x-axis and y-axis
    interval_x = [model.interval_var(size=items[i][0], name=f'interval_x_{i}') for i in range(n)]
    interval_y = [model.interval_var(size=items[i][1], name=f'interval_y_{i}') for i in range(n)]
    
    # Link interval variables with position variables
    for i in range(n):
        model.add(model.start_of(interval_x[i]) == x[i])
        model.add(model.start_of(interval_y[i]) == y[i])
    
    # Define binary variables for relative positions
    lr = {}  # i is left of j
    rl = {}  # j is left of i
    ab = {}  # i is above j
    ba = {}  # j is above i
    
    # Non-overlapping constraints
    for i in range(n):
        for j in range(i + 1, n):
            # Create binary decision variables
            lr[(i,j)] = model.binary_var(name=f'lr_{i}_{j}')
            rl[(i,j)] = model.binary_var(name=f'rl_{i}_{j}')
            ab[(i,j)] = model.binary_var(name=f'ab_{i}_{j}')
            ba[(i,j)] = model.binary_var(name=f'ba_{i}_{j}')
            
            # Add constraints to enforce the meaning of decision variables
            model.add(model.if_then(lr[(i,j)] == 1, 
                                    model.end_of(interval_x[i]) <= model.start_of(interval_x[j])))
            model.add(model.if_then(rl[(i,j)] == 1, 
                                    model.end_of(interval_x[j]) <= model.start_of(interval_x[i])))
            model.add(model.if_then(ab[(i,j)] == 1, 
                                    model.end_of(interval_y[i]) <= model.start_of(interval_y[j])))
            model.add(model.if_then(ba[(i,j)] == 1, 
                                    model.end_of(interval_y[j]) <= model.start_of(interval_y[i])))
            
            # At least one of the relative positions must be true
            model.add(lr[(i,j)] + rl[(i,j)] + ab[(i,j)] + ba[(i,j)] >= 1)
            
            # C2 Symmetry Breaking - Large rectangles
            if items[i][0] + items[j][0] > strip_width:
                # If two rectangles can't fit side by side, they must be stacked vertically
                # Disable horizontal placements
                model.add(lr[(i,j)] == 0)
                model.add(rl[(i,j)] == 0)
    
    # C2 Symmetry Breaking - One pair of rectangles
    if n >= 2:
        # Fix the positional relationship between rectangles 0 and 1
        # Rectangle 1 cannot be to the left of rectangle 0
        model.add(rl[(0,1)] == 0)
        # Rectangle 1 cannot be below rectangle 0
        model.add(ba[(0,1)] == 0)
    
    # Solve the model - TimeLimit is passed directly to solve()
    solution = model.solve(TimeLimit=time_limit)
    
    if solution and solution.is_solution():
        # Extract solution
        positions = [(solution.get_value(x[i]), solution.get_value(y[i])) for i in range(n)]
        return {
            'positions': positions,
        }
    else:
        return None

def solve_2SPP(strip_width, items, time_limit=1800):
    """
    Solves the 2D Strip Packing Problem using bisection search.
    
    Args:
        strip_width: Width of the strip
        items: List of (width, height) tuples for each rectangle
        
    Returns:
        dict with optimal height and positions
    """
    global best_height, best_positions, upper_bound
    
    # Calculate lower and upper bounds
    total_area = sum(w * h for w, h in items)
    max_height = max(h for _, h in items)
    area_lb = math.ceil(total_area / strip_width)
    lb = max(area_lb, max_height)
    
    # Upper bound can be sum of all heights (worst case stacking)
    ub = sum(h for _, h in items)
    upper_bound = ub
    
    # Store the best solution found
    best_solution = None
    optimal_height = ub
    best_height = ub
    
    print(f"Initial bounds: LB={lb}, UB={ub}")
    iter_start_time = time.time()
    
    # Binary search for optimal height
    while lb <= ub:
        mid = (lb + ub) // 2
        print(f"Trying height: {mid}")
        
        # Lưu checkpoint trước khi giải
        save_checkpoint(instance_id, best_height)
            
        # Solve the 2OPP with fixed height = mid
        result = solve_2OPP(strip_width, items, mid, time_limit=time_limit)
        
        if result:  # Feasible solution found
            # Update best solution
            best_solution = result
            optimal_height = mid
            best_height = mid
            best_positions = result['positions']
            
            # Lưu checkpoint sau khi tìm thấy giải pháp tốt hơn
            save_checkpoint(instance_id, best_height)
            
            # Continue searching for better solutions
            ub = mid - 1
            print(f"  Solution found. New UB={ub}")
        else:  # No solution with this height
            # Increase lower bound
            lb = mid + 1
            print(f"  No solution. New LB={lb}")
    
    if best_solution:
        return {
            'height': optimal_height,
            'positions': best_solution['positions'],
            'solve_time': time.time() - iter_start_time
        }
    else:
        return None

def display_solution(strip, rectangles, pos_circuits, instance_name):
    """Visualize the packing solution with improved styling"""
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.title(f"Strip Packing Solution (Width: {strip[0]}, Height: {strip[1]:.2f})")

    if len(pos_circuits) > 0:
        # Use colormap to make visualization more pleasing
        colors = cm.viridis(np.linspace(0, 1, len(rectangles)))
        
        for i in range(len(rectangles)):
            # Create rectangle with nice coloring
            rect = plt.Rectangle(pos_circuits[i], rectangles[i][0], rectangles[i][1], 
                               edgecolor="#333", facecolor=colors[i], alpha=0.7)
            ax.add_patch(rect)
            
            # Add rectangle number at center
            rx, ry = pos_circuits[i]
            cx, cy = rx + rectangles[i][0]/2, ry + rectangles[i][1]/2
            
            # Calculate text color based on background brightness
            rgb = mcolors.to_rgb(colors[i])
            brightness = 0.299*rgb[0] + 0.587*rgb[1] + 0.114*rgb[2]
            text_color = 'white' if brightness < 0.6 else 'black'
            
            ax.annotate(str(i), (cx, cy), color=text_color, 
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
    plt.savefig(f'CPLEX_CP_C2/{instance_name}.png')
    plt.close()

def calculate_first_fit_upper_bound(width, rectangles):
    # Sort rectangles by height (non-increasing)
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

if __name__ == "__main__":
    # Phần controller mode
    if len(sys.argv) == 1:
        # This is the controller mode - running without arguments
        # Create output folder if it doesn't exist
        if not os.path.exists('CPLEX_CP_C2'):
            os.makedirs('CPLEX_CP_C2')
        
        # Đọc file Excel hiện có để kiểm tra instances đã hoàn thành
        excel_file = 'CPLEX_CP_C2.xlsx'
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
            command = f"./runlim --real-time-limit={TIMEOUT} python3 CPLEX_CP_C2.py {instance_id}"
            
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
            
            # Calculate initial bounds
            total_area = sum(w * h for w, h in rectangles)
            max_height = max(h for _, h in rectangles)
            area_lb = math.ceil(total_area / strip_width)
            lower_bound = max(area_lb, max_height)
            
            # Upper bound can be sum of all heights (worst case stacking)
            upper_bound = min(sum(h for _, h in rectangles), calculate_first_fit_upper_bound(strip_width, rectangles))

            print(f"Solving 2D Strip Packing with CPLEX CP for instance {instance_name}")
            print(f"Width: {strip_width}")
            print(f"Number of rectangles: {n_rec}")
            print(f"Lower bound: {lower_bound}")
            print(f"Upper bound: {upper_bound}")
            
            # Solve with CP
            result = solve_2SPP(strip_width, rectangles, time_limit=1800) # Leave 60 seconds for other operations
            
            stop = timeit.default_timer()
            runtime = stop - start

            # Display and save the solution if we found one
            if result:
                optimal_height = result['height']
                display_solution((strip_width, optimal_height), rectangles, result['positions'], instance_name)
            else:
                optimal_height = best_height if best_height != float('inf') else upper_bound
                print(f"No solution found, using best height: {optimal_height}")

            # Tạo result object
            result_data = {
                'Instance': instance_name,
                'Runtime': runtime,
                'Optimal_Height': optimal_height,
                'Status': 'COMPLETE' if result else 'ERROR'
            }
            
            # Ghi kết quả vào Excel trực tiếp
            excel_file = 'CPLEX_CP_C2.xlsx'
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
                'Status': 'ERROR'
            }
            
            # Ghi kết quả lỗi vào Excel
            excel_file = 'CPLEX_CP_C2.xlsx'
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