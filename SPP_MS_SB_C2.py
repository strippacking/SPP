import math
import fileinput
import matplotlib.pyplot as plt
import timeit
import sys
import signal
import pandas as pd
import os
import tempfile
import subprocess
import traceback
import json
import time

from pysat.formula import CNF
from pysat.solvers import Glucose42

# Global variables to track best solution found so far
best_height = float('inf')
best_positions = []
variables_length = 0
clauses_length = 0
upper_bound = 0  # Biến toàn cục để lưu trữ upper_bound

# Signal handler for graceful interruption (e.g., by runlim)
def handle_interrupt(signum, frame):
    print(f"\nReceived interrupt signal {signum}. Saving current best solution.")
    
    # Lấy chiều cao tốt nhất (hoặc là giá trị tìm được, hoặc là upper_bound)
    current_height = best_height if best_height != float('inf') else upper_bound
    print(f"Best height found before interrupt: {current_height}")
    
    # Save result as JSON for the controller to pick up
    result = {
        'Instance': instances[instance_id],  # Thêm tên instance
        'Variables': variables_length,
        'Clauses': clauses_length,
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
if not os.path.exists('SPP_MS_SB_C2'):
    os.makedirs('SPP_MS_SB_C2')

def read_file_instance(n_instance):
    s = ''
    filepath = "inputs/ins-{}.txt".format(n_instance)
    for line in fileinput.input(files=filepath):
        s += line
    return s.splitlines()

instances= [ "",
  "HT01(c1p1)", "HT02(c1p2)", "HT03(c1p3)", "HT04(c2p1)", "HT05(c2p2)", "HT06(c2p3)", 
  "HT07(c3p1)", "HT08(c3p2)", "HT09(c3p3)", 
  "CGCUT01", "CGCUT02", "CGCUT03", 
  "GCUT01", "GCUT02", "GCUT03", "GCUT04", 
  "NGCUT01", "NGCUT02", "NGCUT03", "NGCUT04", "NGCUT05", "NGCUT06", "NGCUT07", 
  "NGCUT08", "NGCUT09", "NGCUT10", "NGCUT11", "NGCUT12", 
  "BENG01", "BENG02", "BENG03", "BENG04", "BENG05", "BENG06", "BENG07", "BENG08", "BENG09", "BENG10", "HT10(c4p1)", "HT11(c4p2)", "HT12(c4p3)"
]

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

def display_solution(strip, rectangles, pos_circuits, instance_name):
    # define Matplotlib figure and axis
    fig, ax = plt.subplots()
    ax = plt.gca()
    plt.title(strip)

    if len(pos_circuits) > 0:
        for i in range(len(rectangles)):
            rect = plt.Rectangle(pos_circuits[i], *rectangles[i], edgecolor="#333")
            ax.add_patch(rect)

    ax.set_xlim(0, strip[0])
    ax.set_ylim(0, strip[1] + 1)
    ax.set_xticks(range(strip[0] + 1))
    ax.set_yticks(range(strip[1] + 1))
    ax.set_xlabel('width')
    ax.set_ylabel('height')
    
    # Save the plot to SPP folder
    plt.savefig(f'SPP_MS_SB_C2/{instance_name}.png')
    plt.close()

def positive_range(end):
    if end < 0:
        return []
    return range(end)

# Thêm hàm save_checkpoint để lưu tiến trình giải
def save_checkpoint(instance_id, variables, clauses, height, status="IN_PROGRESS"):
    checkpoint = {
        'Variables': variables,
        'Clauses': clauses,
        'Runtime': timeit.default_timer() - start,
        'Optimal_Height': height if height != float('inf') else upper_bound,
        'Status': status
    }
    
    # Ghi ra file checkpoint
    with open(f'checkpoint_{instance_id}.json', 'w') as f:
        json.dump(checkpoint, f)

def SPP_MaxSAT(width, rectangles, lower_bound, upper_bound):
    """Solve 2SPP using Max-SAT approach with tt-open-wbo-inc"""
    global variables_length, clauses_length, best_height, best_positions
    n_rectangles = len(rectangles)
    variables = {}
    counter = 1

    # Create a temporary file for the Max-SAT formula
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.wcnf') as file:
        wcnf_file = file.name
        
        # Add comments for clarity (optional)
        file.write(f"c Strip Packing Problem, W={width}, n={n_rectangles}\n")
        file.write(f"c Lower bound={lower_bound}, Upper bound={upper_bound}\n")
        
        # Define variables for rectangle positions and relations
        for i in range(n_rectangles):
            for j in range(n_rectangles):
                if i != j:
                    variables[f"lr{i + 1},{j + 1}"] = counter  # lri,rj
                    counter += 1
                    variables[f"ud{i + 1},{j + 1}"] = counter  # uri,rj
                    counter += 1
            for e in range(width - rectangles[i][0] + 2):  # Position variables for x-axis
                variables[f"px{i + 1},{e}"] = counter  # pxi,e
                counter += 1
            for h in range(upper_bound - rectangles[i][1] + 2):  # Position variables for y-axis
                variables[f"py{i + 1},{h}"] = counter  # pyi,h
                counter += 1
        
        # Height constraint variables - ph_h means "can pack with height ≤ h"
        for h in range(lower_bound, upper_bound + 1):
            variables[f"ph_{h}"] = counter
            counter += 1
            
        # Prepare hard clauses (basic packing constraints)
        hard_clauses = []
        
        # Order encoding axioms
        for i in range(n_rectangles):
            for e in range(width - rectangles[i][0] + 1):
                hard_clauses.append([-variables[f"px{i + 1},{e}"], variables[f"px{i + 1},{e + 1}"]])
            for h in range(upper_bound - rectangles[i][1] + 1):
                hard_clauses.append([-variables[f"py{i + 1},{h}"], variables[f"py{i + 1},{h + 1}"]])
        
        # Height variable ordering - this enforces that ph_h implies ph_{h+1}
        for h in range(lower_bound, upper_bound):
            hard_clauses.append([-variables[f"ph_{h}"], variables[f"ph_{h+1}"]])
        
        # Non-overlapping constraints function remains the same
        def non_overlapping(i, j, h1, h2, v1, v2):
            i_width = rectangles[i][0]
            i_height = rectangles[i][1]
            j_width = rectangles[j][0]
            j_height = rectangles[j][1]
        
            # lri,j v lrj,i v udi,j v udj,i
            four_literal = []
            if h1: four_literal.append(variables[f"lr{i + 1},{j + 1}"])
            if h2: four_literal.append(variables[f"lr{j + 1},{i + 1}"])
            if v1: four_literal.append(variables[f"ud{i + 1},{j + 1}"])
            if v2: four_literal.append(variables[f"ud{j + 1},{i + 1}"])
            hard_clauses.append(four_literal)
        
            # First type of constraints - prevent rectangle j's left edge from being in first i_width positions
            if h1:
                for e in range(i_width):
                    if f"px{j + 1},{e}" in variables:
                        hard_clauses.append([-variables[f"lr{i + 1},{j + 1}"], -variables[f"px{j + 1},{e}"]])
            
            # First type for h2 - prevent rectangle i's left edge from being in first j_width positions
            if h2:
                for e in range(j_width):
                    if f"px{i + 1},{e}" in variables:
                        hard_clauses.append([-variables[f"lr{j + 1},{i + 1}"], -variables[f"px{i + 1},{e}"]])
        
            # First type for v1 - prevent rectangle j's bottom edge from being in first i_height positions
            if v1:
                for y_pos in range(i_height):
                    if f"py{j + 1},{y_pos}" in variables:
                        hard_clauses.append([-variables[f"ud{i + 1},{j + 1}"], -variables[f"py{j + 1},{y_pos}"]])
            
            # First type for v2 - prevent rectangle i's bottom edge from being in first j_height positions
            if v2:
                for y_pos in range(j_height):
                    if f"py{i + 1},{y_pos}" in variables:
                        hard_clauses.append([-variables[f"ud{j + 1},{i + 1}"], -variables[f"py{i + 1},{y_pos}"]])
        
            # Second type of constraints - enforce relative positions when certain relations hold
            # For h1: if rectangle i is to the left of j, then i's left edge at e implies j can't be at e+i_width
            if h1:
                for e in positive_range(width - i_width):
                    if f"px{i + 1},{e}" in variables and f"px{j + 1},{e + i_width}" in variables:
                        hard_clauses.append([-variables[f"lr{i + 1},{j + 1}"],
                                      variables[f"px{i + 1},{e}"], 
                                      -variables[f"px{j + 1},{e + i_width}"]])
            
            # For h2: if rectangle j is to the left of i, then j's left edge at e implies i can't be at e+j_width
            if h2:
                for e in positive_range(width - j_width):
                    if f"px{j + 1},{e}" in variables and f"px{i + 1},{e + j_width}" in variables:
                        hard_clauses.append([-variables[f"lr{j + 1},{i + 1}"],
                                      variables[f"px{j + 1},{e}"], 
                                      -variables[f"px{i + 1},{e + j_width}"]])
        
            # For v1: if rectangle i is below j, then i's bottom edge at f implies j can't be at f+i_height
            if v1:
                for y_pos in positive_range(upper_bound - i_height):
                    if f"py{i + 1},{y_pos}" in variables and f"py{j + 1},{y_pos + i_height}" in variables:
                        hard_clauses.append([-variables[f"ud{i + 1},{j + 1}"],
                                      variables[f"py{i + 1},{y_pos}"], 
                                      -variables[f"py{j + 1},{y_pos + i_height}"]])
            
            # For v2: if rectangle j is below i, then j's bottom edge at f implies i can't be at f+j_height
            if v2:
                for y_pos in positive_range(upper_bound - j_height):
                    if f"py{j + 1},{y_pos}" in variables and f"py{i + 1},{y_pos + j_height}" in variables:
                        hard_clauses.append([-variables[f"ud{j + 1},{i + 1}"],
                                      variables[f"py{j + 1},{y_pos}"], 
                                      -variables[f"py{i + 1},{y_pos + j_height}"]])
        
        # Find max height and width for symmetry breaking
        max_height = max([int(rectangle[1]) for rectangle in rectangles])
        max_width = max([int(rectangle[0]) for rectangle in rectangles])
        second_max_width = max([int(rectangle[0]) for rectangle in rectangles if int(rectangle[0]) != max_width])

        # Symmetry Breaking - Config 2
        for i in range(n_rectangles):
            for j in range(i + 1, n_rectangles):
                #Fix the position of the largest rectangle and the second largest rectangle
                if rectangles[i][0] == max_width and rectangles[j][0] == second_max_width:
                    non_overlapping(i, j, True, False, True, False)
                # Large-rectangles horizontal
                elif rectangles[i][0] + rectangles[j][0] > width:
                    non_overlapping(i, j, False, False, True, True)
                # Large-rectangles vertical
                elif rectangles[i][1] + rectangles[j][1] > upper_bound:
                    non_overlapping(i, j, True, True, False, False)
                # Same-sized rectangles
                elif rectangles[i] == rectangles[j]:
                    non_overlapping(i, j, True, False, True, True)
                else:
                    non_overlapping(i, j, True, True, True, True)
        
        # Domain encoding to ensure rectangles are placed within strip boundaries
        for i in range(n_rectangles):
            # Right edge within strip width
            hard_clauses.append([variables[f"px{i + 1},{width - rectangles[i][0]}"]])
            
            # Top edge within strip height - integrated with height constraint variables
            for h in range(lower_bound, upper_bound + 1):
                if h >= rectangles[i][1]:  # Rectangle must fit below height h
                    hard_clauses.append([-variables[f"ph_{h}"], variables[f"py{i + 1},{h - rectangles[i][1]}"]])
        
        # Prepare soft clauses with weights for height minimization
        soft_clauses = []
        
        # Use weight 1 for all height variables
        for h in range(lower_bound, upper_bound + 1):
            soft_clauses.append((1, [variables[f"ph_{h}"]]))  # We want ph_h to be FALSE when possible
        
        # Require at least one ph_h to be true (ensures a solution exists)
        all_ph_vars = [variables[f"ph_{h}"] for h in range(lower_bound, upper_bound + 1)]
        hard_clauses.append(all_ph_vars)
        
        # Write hard clauses with 'h' prefix
        for clause in hard_clauses:
            file.write(f"h {' '.join(map(str, clause))} 0\n")
        
        # Write soft clauses with their weights
        for weight, clause in soft_clauses:
            file.write(f"{weight} {' '.join(map(str, clause))} 0\n")
            
        # For debugging, print details about the WCNF file
        print(f"Created WCNF file with {len(hard_clauses)} hard clauses and {len(soft_clauses)} soft clauses")
        print(f"Variable count: {counter-1}")
        print(f"Sample variables: ph_{lower_bound}={variables[f'ph_{lower_bound}']}, " +
              f"px{1},{0}={variables.get(f'px{1},{0}', 'N/A')}")
        
        # Flush and close the file
        file.flush()

    variables_length = len(variables)
    clauses_length = len(hard_clauses) + len(soft_clauses)
    
    # Lưu checkpoint trước khi giải
    save_checkpoint(instance_id, variables_length, clauses_length, 
                   best_height if best_height != float('inf') else upper_bound)
            
    # Call tt-open-wbo-inc solver
    try:
        print(f"Running tt-open-wbo-inc on {wcnf_file}...")
        result = subprocess.run(
            ["./tt-open-wbo-inc-Glucose4_1_static", wcnf_file], 
            capture_output=True, 
            text=True
        )
        
        output = result.stdout
        print(f"Solver output preview: {output[:200]}...")  # Debug: Show beginning of output
        
        # Parse the output to find the model
        optimal_height = upper_bound
        positions = [[0, 0] for _ in range(n_rectangles)]
        
        if "OPTIMUM FOUND" in output:
            print("Optimal solution found!")
            
            # Extract the model line (starts with "v ")
            for line in output.split('\n'):
                if line.startswith('v '):
                    print(f"Found solution line: {line[:50]}...")  # Debug output
                    
                    # Remove the 'v ' prefix
                    binary_string = line[2:].strip()
                    
                    # Convert binary values to true variable set
                    true_vars = set()
                    
                    # Process the solution string - tt-open-wbo-inc can output in different formats
                    # Try to interpret as space-separated list of integers first
                    if " " in binary_string:
                        try:
                            for val in binary_string.split():
                                val_int = int(val)
                                if val_int > 0:  # Positive literals represent true variables
                                    true_vars.add(val_int)
                        except ValueError:
                            # Not integers, try as space-separated binary values
                            for i, val in enumerate(binary_string.split()):
                                if val == '1':
                                    true_vars.add(i + 1)  # 1-indexed
                    else:
                        # No spaces - treat as continuous binary string
                        for i, val in enumerate(binary_string):
                            if val == '1':
                                true_vars.add(i + 1)  # 1-indexed
                    
                    # Extract height variables and find minimum height where ph_h is true
                    ph_true_heights = []
                    for h in range(lower_bound, upper_bound + 1):
                        var_key = f"ph_{h}"
                        if var_key in variables and variables[var_key] in true_vars:
                            ph_true_heights.append(h)
                    
                    if ph_true_heights:
                        optimal_height = min(ph_true_heights)
                        print(f"Heights where ph_h is true: {sorted(ph_true_heights)}")
                        # Update best_height global variable
                        if optimal_height < best_height:
                            best_height = optimal_height
                            save_checkpoint(instance_id, variables_length, clauses_length, best_height)
                    else:
                        print("WARNING: No ph_h variables are true! This may indicate a parsing issue.")
                        # Check if we're within bounds - the solution string might not include all variables
                        height_var_indices = [variables[f"ph_{h}"] for h in range(lower_bound, upper_bound + 1)]
                        min_height_var = min(height_var_indices)
                        max_height_var = max(height_var_indices)
                        
                        if len(binary_string) < min_height_var:
                            print(f"Binary string length ({len(binary_string)}) is less than smallest height variable index ({min_height_var}).")
                            print("This suggests the output format needs to be interpreted differently.")
                            # Assume lowest possible height when uncertain
                            optimal_height = lower_bound
                            best_height = lower_bound

                    # If we couldn't parse any variables but the solver found a solution,
                    # use the lower bound as a fallback
                    if not true_vars:
                        print("WARNING: Solution parsing failed. Using lower bound height as fallback.")
                        optimal_height = lower_bound
                        best_height = lower_bound
                        
                        # Set default positions - simple greedy left-bottom placement
                        x_pos = 0
                        y_pos = 0
                        max_height = 0
                        for i in range(n_rectangles):
                            # If current rectangle doesn't fit in the current row, move to next row
                            if x_pos + rectangles[i][0] > width:
                                x_pos = 0
                                y_pos = max_height
                            
                            positions[i][0] = x_pos
                            positions[i][1] = y_pos
                            
                            # Update position for next rectangle
                            x_pos += rectangles[i][0]
                            max_height = max(max_height, y_pos + rectangles[i][1])
                    
                    # Extract positions - Find the exact transition point for each rectangle
                    for i in range(n_rectangles):
                        # Find x position (first position where px is true)
                        found_x = False
                        for e in range(width - rectangles[i][0] + 1):
                            var_key = f"px{i + 1},{e}"
                            if var_key in variables and variables[var_key] in true_vars:
                                if e == 0 or f"px{i + 1},{e-1}" not in variables or variables[f"px{i + 1},{e-1}"] not in true_vars:
                                    positions[i][0] = e
                                    found_x = True
                                    break
                                    
                        # Find y position (first position where py is true)
                        found_y = False
                        for y_pos in range(upper_bound - rectangles[i][1] + 1):
                            var_key = f"py{i + 1},{y_pos}"
                            if var_key in variables and variables[var_key] in true_vars:
                                if y_pos == 0 or f"py{i + 1},{y_pos-1}" not in variables or variables[f"py{i + 1},{y_pos-1}"] not in true_vars:
                                    positions[i][1] = y_pos
                                    found_y = True
                                    break
                
                # Save best positions
                best_positions = positions
        else:
            print("No optimal solution found.")
            print(f"Solver output: {output}")
        
        # Clean up the temporary file
        os.unlink(wcnf_file)
        return optimal_height, positions
    
    except Exception as e:
        print(f"Error running Max-SAT solver: {e}")
        traceback.print_exc()  # Print traceback for more detailed error information
        if os.path.exists(wcnf_file):
            os.unlink(wcnf_file)
        return None, None

if __name__ == "__main__":
    # Phần controller mode
    if len(sys.argv) == 1:
        # This is the controller mode - running without arguments
        # Create SPP folder if it doesn't exist
        if not os.path.exists('SPP_MS_SB_C2'):
            os.makedirs('SPP_MS_SB_C2')
        
        # Đọc file Excel hiện có để kiểm tra instances đã hoàn thành
        excel_file = 'SPP_MS_SB_C2.xlsx'
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
        TIMEOUT = 900  # 30 minutes timeout
        
        for instance_id in range(39, 42):
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
            command = f"./runlim --time-limit={TIMEOUT} python3 SPP_MS_SB_C2.py {instance_id}"
            
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
            input = read_file_instance(instance_id)
            width = int(input[0])
            n_rec = int(input[1])
            rectangles = []
            for i in range(2, 2 + n_rec):
                if i < len(input):
                    rect = [int(val) for val in input[i].split()]
                    rectangles.append(rect)
                else:
                    raise IndexError(f"Missing rectangle data at line {i}")
            
            # Calculate initial bounds
            heights = [int(rectangle[1]) for rectangle in rectangles]
            upper_bound = min(sum(heights), calculate_first_fit_upper_bound(width, rectangles))
            lower_bound = max(math.ceil(sum([int(rectangle[0] * rectangle[1]) for rectangle in rectangles]) / width), max(heights))

            print(f"Solving 2D Strip Packing with MaxSAT for instance {instance_name}")
            print(f"Width: {width}")
            print(f"Number of rectangles: {n_rec}")
            print(f"Lower bound: {lower_bound}")
            print(f"Upper bound: {upper_bound}")
            
            # Solve with MaxSAT
            optimal_height, optimal_pos = SPP_MaxSAT(width, rectangles, lower_bound, upper_bound)
            
            stop = timeit.default_timer()
            runtime = stop - start

            # Display and save the solution if we found one
            if optimal_height is not None and optimal_pos is not None:
                display_solution((width, optimal_height), rectangles, optimal_pos, instance_name)

            # Tạo result object
            result = {
                'Instance': instance_name,
                'Variables': variables_length,
                'Clauses': clauses_length,
                'Runtime': runtime,
                'Optimal_Height': optimal_height if optimal_height is not None else upper_bound,
                'Status': 'COMPLETE' if optimal_height is not None else 'ERROR'
            }
            
            # Ghi kết quả vào Excel trực tiếp
            excel_file = 'SPP_MS_SB_C2.xlsx'
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
            
            # Save result to a JSON file that the controller will read
            with open(f'results_{instance_id}.json', 'w') as f:
                json.dump(result, f)
            
            print(f"Instance {instance_name} completed - Runtime: {runtime:.2f}s, Height: {optimal_height}")

        except Exception as e:
            print(f"Error in instance {instance_name}: {str(e)}")
            traceback.print_exc()  # Print the traceback for the error
            
            # Save error result - use upper_bound if no best_height
            current_height = best_height if best_height != float('inf') else upper_bound
            result = {
                'Instance': instance_name,
                'Variables': variables_length,
                'Clauses': clauses_length,
                'Runtime': timeit.default_timer() - start,
                'Optimal_Height': current_height,
                'Status': 'ERROR'
            }
            
            # Ghi kết quả lỗi vào Excel
            excel_file = 'SPP_MS_SB_C2.xlsx'
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
                except Exception as ex:
                    print(f"Lỗi khi đọc file Excel hiện có: {str(ex)}")
                    existing_df = pd.DataFrame([result])
            else:
                # Tạo DataFrame mới nếu chưa có file Excel
                existing_df = pd.DataFrame([result])
            
            # Lưu DataFrame vào Excel
            existing_df.to_excel(excel_file, index=False)
            print(f"Error results saved to {excel_file}")
            
            # Save result to a JSON file that the controller will read
            with open(f'results_{instance_id}.json', 'w') as f:
                json.dump(result, f)