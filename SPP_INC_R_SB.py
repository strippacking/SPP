import math
import fileinput
import matplotlib.pyplot as plt
import timeit
import sys
import signal
import pandas as pd
import os
import json
import subprocess
import time

from pysat.formula import CNF
from pysat.solvers import Glucose42

# Global variables to track best solution found so far
best_height = float('inf')
best_positions = []
best_rotations = []
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
if not os.path.exists('SPP_INC_R_SB'):
    os.makedirs('SPP_INC_R_SB')

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

def display_solution(strip, rectangles, pos_circuits, rotations, instance_name):
    fig, ax = plt.subplots()
    plt.title(f"Strip Packing Solution (Width: {strip[0]}, Height: {strip[1]})")

    if len(pos_circuits) > 0:
        for i in range(len(rectangles)):
            w = rectangles[i][1] if rotations[i] else rectangles[i][0]
            h = rectangles[i][0] if rotations[i] else rectangles[i][1]
            rect = plt.Rectangle(pos_circuits[i], w, h, 
                               edgecolor="#333", facecolor="lightblue", alpha=0.6)
            ax.add_patch(rect)
            rx, ry = pos_circuits[i]
            cx, cy = rx + w/2, ry + h/2
            ax.annotate(str(i), (cx, cy), color='black', ha='center', va='center')

    ax.set_xlim(0, strip[0])
    ax.set_ylim(0, strip[1] + 1)
    ax.set_xticks(range(strip[0] + 1))
    ax.set_yticks(range(int(strip[1]) + 2))
    ax.set_xlabel('Width')
    ax.set_ylabel('Height')
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Save the plot to SPP folder
    plt.savefig(f'SPP_INC_R_SB/{instance_name}.png')
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

def SPP_Incremental_Rotation(rectangles, strip_width, lower_bound, upper_bound):
    """
    Solve 2D Strip Packing with rotation using incremental SAT solving.
    Returns the optimal height, rectangle positions, and rotation states.
    """
    global variables_length, clauses_length, best_height, best_positions, best_rotations
    n_rectangles = len(rectangles)
    
    # Initialize the CNF formula and variables
    cnf = CNF()
    variables = {}
    counter = 1
    
    # Create variables for rectangle positions and relations
    # lr (left-right) and ud (up-down) variables
    for i in range(n_rectangles):
        for j in range(n_rectangles):
            if i != j:
                variables[f"lr{i+1},{j+1}"] = counter  # lri,rj
                counter += 1
                variables[f"ud{i+1},{j+1}"] = counter  # udi,rj
                counter += 1
        
        # Position variables with order encoding
        for e in range(strip_width):
            variables[f"px{i+1},{e}"] = counter  # pxi,e
            counter += 1
            
        for h in range(upper_bound):
            variables[f"py{i+1},{h}"] = counter  # pyi,h
            counter += 1
        
        # Rotation variables - True means the rectangle is rotated (width↔height)
        variables[f"r{i+1}"] = counter
        counter += 1
    
    # Height variables - ph_h means "can pack with height ≤ h"
    for h in range(lower_bound, upper_bound + 1):
        variables[f"ph_{h}"] = counter
        counter += 1
    
    # Add order encoding axiom clauses
    for i in range(n_rectangles):
        # ¬pxi,e ∨ pxi,e+1
        for e in range(strip_width - 1):
            cnf.append([-variables[f"px{i+1},{e}"], variables[f"px{i+1},{e+1}"]])
        
        # ¬pyi,h ∨ pyi,h+1
        for h in range(upper_bound - 1):
            cnf.append([-variables[f"py{i+1},{h}"], variables[f"py{i+1},{h+1}"]])
    
    # Add height variable ordering constraints
    # If ph_h is true, then ph_{h+1} must also be true
    for h in range(lower_bound, upper_bound):
        cnf.append([-variables[f"ph_{h}"], variables[f"ph_{h+1}"]])
    
    # Non-overlapping constraints function
    def add_non_overlapping(rotated, i, j, h1, h2, v1, v2):
        # Get dimensions based on rotation
        if not rotated:
            i_width = rectangles[i][0]
            i_height = rectangles[i][1]
            j_width = rectangles[j][0]
            j_height = rectangles[j][1]
            i_rotation = variables[f"r{i+1}"]  # Not rotated means r{i+1} is False
            j_rotation = variables[f"r{j+1}"]  # Not rotated means r{j+1} is False
        else:
            i_width = rectangles[i][1]
            i_height = rectangles[i][0]
            j_width = rectangles[j][1]
            j_height = rectangles[j][0]
            i_rotation = -variables[f"r{i+1}"]   # Rotated means r{i+1} is True
            j_rotation = -variables[f"r{j+1}"]   # Rotated means r{j+1} is True

        # lri,j v lrj,i v udi,j v udj,i v rotation literals
        four_literal = []
        if h1: four_literal.append(variables[f"lr{i+1},{j+1}"])
        if h2: four_literal.append(variables[f"lr{j+1},{i+1}"])
        if v1: four_literal.append(variables[f"ud{i+1},{j+1}"])
        if v2: four_literal.append(variables[f"ud{j+1},{i+1}"])
        
        # Make clause only apply if rotation is correct
        cnf.append(four_literal + [i_rotation])
        cnf.append(four_literal + [j_rotation])

        # First type constraints - prevent rectangle j's left edge from being too close to i
        if h1:
            for e in range(min(strip_width, i_width)):
                cnf.append([i_rotation, -variables[f"lr{i+1},{j+1}"], -variables[f"px{j+1},{e}"]])
        
            for e in positive_range(strip_width - i_width):
                cnf.append([i_rotation, -variables[f"lr{i+1},{j+1}"],
                          variables[f"px{i+1},{e}"], 
                          -variables[f"px{j+1},{e + i_width}"]])
        
        if h2:
            for e in range(min(strip_width, j_width)):
                cnf.append([j_rotation, -variables[f"lr{j+1},{i+1}"], -variables[f"px{i+1},{e}"]])
            
            for e in positive_range(strip_width - j_width):
                cnf.append([j_rotation, -variables[f"lr{j+1},{i+1}"],
                          variables[f"px{j+1},{e}"], 
                          -variables[f"px{i+1},{e + j_width}"]])

        if v1:
            for h in range(min(upper_bound, i_height)):
                cnf.append([i_rotation, -variables[f"ud{i+1},{j+1}"], -variables[f"py{j+1},{h}"]])
            
            for h in positive_range(upper_bound - i_height):
                cnf.append([i_rotation, -variables[f"ud{i+1},{j+1}"],
                          variables[f"py{i+1},{h}"], 
                          -variables[f"py{j+1},{h + i_height}"]])
        
        if v2:
            for h in range(min(upper_bound, j_height)):
                cnf.append([j_rotation, -variables[f"ud{j+1},{i+1}"], -variables[f"py{i+1},{h}"]])
            
            for h in positive_range(upper_bound - j_height):
                cnf.append([j_rotation, -variables[f"ud{j+1},{i+1}"],
                          variables[f"py{j+1},{h}"], 
                          -variables[f"py{i+1},{h + j_height}"]])
    
    max_width = max([int(rectangle[0]) for rectangle in rectangles])
    second_max_width = max([int(rectangle[0]) for rectangle in rectangles if int(rectangle[0]) != max_width])

    #Symmetry Breaking
    for i in range(len(rectangles)):
        for j in range(i + 1, len(rectangles)):
            #Fix the position of the largest rectangle and the second largest rectangle
            if rectangles[i][0] == max_width and rectangles[j][0] == second_max_width:
                add_non_overlapping(False, i, j, False, False, True, True)
                add_non_overlapping(True, i, j, False, False, True, True)
            # Large-rectangles horizontal
            if min(rectangles[i][0], rectangles[i][1]) + min(rectangles[j][0], rectangles[j][1]) > strip_width:
                add_non_overlapping(False, i, j, False, False, True, True)
                add_non_overlapping(True, i, j, False, False, True, True)
            # Large rectangles vertical
            elif min(rectangles[i][0], rectangles[i][1]) + min(rectangles[j][0], rectangles[j][1]) > upper_bound:
                add_non_overlapping(False, i, j, True, True, False, False)
                add_non_overlapping(True, i, j, True, True, False, False)
            # Normal rectangles
            else:
                add_non_overlapping(False, i, j, True, True, True, True)
                add_non_overlapping(True, i, j, True, True, True, True)
    
    # Domain encoding to ensure every rectangle stays inside strip's boundary
    for i in range(n_rectangles):
        # Width checks - normal orientation
        if rectangles[i][0] > strip_width:
            cnf.append([variables[f"r{i+1}"]])  # Must be rotated if width > strip_width
        else:
            for e in range(strip_width - rectangles[i][0], strip_width):
                cnf.append([variables[f"r{i+1}"], variables[f"px{i+1},{e}"]])
        
        # Height checks - normal orientation
        if rectangles[i][1] > upper_bound:
            cnf.append([variables[f"r{i+1}"]])  # Must be rotated if height > upper_bound
        else:
            for h in range(upper_bound - rectangles[i][1], upper_bound):
                cnf.append([variables[f"r{i+1}"], variables[f"py{i+1},{h}"]])

        # Width checks - rotated orientation
        if rectangles[i][1] > strip_width:
            cnf.append([-variables[f"r{i+1}"]])  # Cannot be rotated if rotated_width > strip_width
        else:
            for e in range(strip_width - rectangles[i][1], strip_width):
                cnf.append([-variables[f"r{i+1}"], variables[f"px{i+1},{e}"]])
        
        # Height checks - rotated orientation
        if rectangles[i][0] > upper_bound:
            cnf.append([-variables[f"r{i+1}"]])  # Cannot be rotated if rotated_height > upper_bound
        else:
            for h in range(upper_bound - rectangles[i][0], upper_bound):
                cnf.append([-variables[f"r{i+1}"], variables[f"py{i+1},{h}"]])
    
    # Height constraints - connecting ph_h with rectangle positions
    for h in range(lower_bound, upper_bound + 1):
        for i in range(n_rectangles):
            # Normal orientation
            rect_height = rectangles[i][1]
            if h >= rect_height:
                cnf.append([-variables[f"ph_{h}"], variables[f"r{i+1}"], 
                          variables[f"py{i+1},{h - rect_height}"]])
            
            # Rotated orientation
            rotated_height = rectangles[i][0]
            if h >= rotated_height:
                cnf.append([-variables[f"ph_{h}"], -variables[f"r{i+1}"], 
                          variables[f"py{i+1},{h - rotated_height}"]])
                
    variables_length = len(variables)
    clauses_length = len(cnf.clauses)
    
    # Lưu checkpoint ban đầu
    save_checkpoint(instance_id, variables_length, clauses_length, upper_bound)
    
    # Initialize the incremental SAT solver with the CNF formula
    with Glucose42(bootstrap_with=cnf) as solver:
        optimal_height = upper_bound
        positions = None
        rotations = [False] * n_rectangles
        
        # For model reuse (as described in the paper)
        best_model = None
        
        # Binary search with incremental solving
        current_lb = lower_bound
        current_ub = upper_bound
        
        while current_lb <= current_ub:
            mid = (current_lb + current_ub) // 2
            print(f"Trying height: {mid} (lower={current_lb}, upper={current_ub})")
            
            # Set up assumptions for this iteration - test if we can pack with height ≤ mid
            assumptions = [variables[f"ph_{mid}"]]
            
            # Lưu checkpoint trước khi giải
            save_checkpoint(instance_id, variables_length, clauses_length, 
                         best_height if best_height != float('inf') else upper_bound)
            
            # Solve with assumptions
            is_sat = solver.solve(assumptions=assumptions)
            
            if is_sat:
                # We found a solution with height ≤ mid
                optimal_height = mid
                best_height = mid

                save_checkpoint(instance_id, variables_length, clauses_length, best_height)
                
                # Save the model for reuse
                best_model = solver.get_model()
                
                # Extract positions and rotations from the model
                positions = [[0, 0] for _ in range(n_rectangles)]
                rotations = [False for _ in range(n_rectangles)]
                
                # Convert model to dictionary for faster lookup
                true_vars = set(var for var in best_model if var > 0)
                
                # Extract rotation variables
                for i in range(n_rectangles):
                    rotations[i] = variables[f"r{i+1}"] in true_vars
                
                # Extract positions
                for i in range(n_rectangles):
                    # Find x position (first position where px is true)
                    found_x = False
                    for e in range(strip_width):
                        var = variables.get(f"px{i+1},{e}", None)
                        if var in true_vars:
                            if e == 0 or variables[f"px{i+1},{e-1}"] not in true_vars:
                                positions[i][0] = e
                                found_x = True
                                break
                    if not found_x:
                        print(f"WARNING: Could not determine x-position for rectangle {i}!")
                    
                    # Find y position (first position where py is true)
                    found_y = False
                    for y in range(upper_bound):
                        var = variables.get(f"py{i+1},{y}", None)
                        if var in true_vars:
                            if y == 0 or variables[f"py{i+1},{y-1}"] not in true_vars:
                                positions[i][1] = y
                                found_y = True
                                break
                    if not found_y:
                        print(f"WARNING: Could not determine y-position for rectangle {i}!")
                
                # Save the best positions and rotations
                best_positions = positions
                best_rotations = rotations
                
                # Update search range - try lower height
                current_ub = mid - 1
            
            else:
                # No solution with height ≤ mid
                # Update search range - try higher height
                current_lb = mid + 1

        variables_length = len(variables)
        clauses_length = len(cnf.clauses)

        # Final validation of the solution
        if positions is None:
            return None, None, None
        
        actual_max_height = 0
        for i in range(n_rectangles):
            rect_height = rectangles[i][0] if rotations[i] else rectangles[i][1]
            top_edge = positions[i][1] + rect_height
            actual_max_height = max(actual_max_height, top_edge)
        
        if actual_max_height > optimal_height:
            print(f"WARNING: Actual max height {actual_max_height} exceeds optimal_height {optimal_height}")
            optimal_height = actual_max_height
        
        print(f"Final optimal height: {optimal_height}")
        return optimal_height, positions, rotations

if __name__ == "__main__":
    # Phần controller mode
    if len(sys.argv) == 1:
        # This is the controller mode - running without arguments
        # Create SPP folder if it doesn't exist
        if not os.path.exists('SPP_INC_R_SB'):
            os.makedirs('SPP_INC_R_SB')
        
        # Đọc file Excel hiện có để kiểm tra instances đã hoàn thành
        excel_file = 'SPP_INC_R_SB.xlsx'
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
            command = f"./runlim --time-limit={TIMEOUT} python3 SPP_INC_R_SB.py {instance_id}"
            
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
                    
                    # Nếu status là TIMEOUT, ghi vào Excel
                    if result['Status'] == 'TIMEOUT':
                        # Thêm tên instance vào kết quả nếu chưa có
                        result['Runtime'] = "TIMEOUT"
                        if 'Instance' not in result:
                            result['Instance'] = instance_name
                        
                        # Các phần code xử lý Excel
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
                        print(f"Timeout results saved to {excel_file}")
                        
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
            
            # Calculate initial bounds - cập nhật biến toàn cục
            heights = [int(rectangle[1]) for rectangle in rectangles]
            upper_bound = min(sum(heights), calculate_first_fit_upper_bound(width, rectangles))
            lower_bound = max(math.ceil(sum([int(rectangle[0] * rectangle[1]) for rectangle in rectangles]) / width), max(heights))

            print(f"Solving 2D Strip Packing with rotation and incremental SAT for instance {instance_name}")
            print(f"Width: {width}")
            print(f"Number of rectangles: {n_rec}")
            print(f"Lower bound: {lower_bound}")
            print(f"Upper bound: {upper_bound}")
            
            # Solve with incremental SAT and rotation
            optimal_height, optimal_pos, optimal_rot = SPP_Incremental_Rotation(rectangles, width, lower_bound, upper_bound)
            
            stop = timeit.default_timer()
            runtime = stop - start

            # Display and save the solution if we found one
            if optimal_height is not None:
                display_solution((width, optimal_height), rectangles, optimal_pos, optimal_rot, instance_name)

            # Tạo result object
            result = {
                'Instance': instance_name,
                'Variables': variables_length,
                'Clauses': clauses_length,
                'Runtime': runtime,
                'Optimal_Height': optimal_height if optimal_height is not None else upper_bound,
                'Status': 'COMPLETE'
            }
            
            # Ghi kết quả vào Excel trực tiếp
            excel_file = 'SPP_INC_R_SB.xlsx'
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
            # Save error result - sử dụng upper_bound nếu không có best_height
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
            excel_file = 'SPP_INC_R_SB.xlsx'
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