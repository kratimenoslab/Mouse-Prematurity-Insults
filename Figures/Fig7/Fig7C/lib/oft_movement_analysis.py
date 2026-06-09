"""
OFT Analysis - Movement Analysis Module
Functions for analyzing mouse movement patterns
"""

import numpy as np
import pandas as pd
import os
from scipy.ndimage import gaussian_filter1d


def load_boundary_overrides(experiment_name):
    """
    Load boundary overrides from CSV file
    
    Parameters:
    - experiment_name: name of the experiment to look up
    
    Returns:
    - dict with override values or None if no overrides found
    """
    if not experiment_name:
        return None
        
    override_file = r"C:\pt-social-project-data\logs\dlc_boundaries.csv"
    
    if not os.path.exists(override_file):
        print(f"Debug - No boundary override file found at {override_file}")
        return None
    
    try:
        df = pd.read_csv(override_file)
        
        # Find the row for this experiment
        experiment_row = df[df['experiment'] == experiment_name]
        
        if len(experiment_row) == 0:
            print(f"Debug - No boundary overrides found for experiment '{experiment_name}'")
            return None
        
        row = experiment_row.iloc[0]
        overrides = {}
        
        # Check each boundary value
        for col in ['min_x_dlc', 'max_x_dlc', 'min_y_dlc', 'max_y_dlc']:
            if col in row:
                value = row[col]
                # Check if value is valid (not NaN, not 0, not "NA", not empty)
                if pd.isna(value) or value == 0 or str(value).upper() == 'NA' or str(value).strip() == '':
                    overrides[col] = None  # Will trigger calculation
                else:
                    try:
                        overrides[col] = float(value)
                        print(f"Debug - Found override for {col}: {overrides[col]:.1f}")
                    except ValueError:
                        overrides[col] = None  # Invalid value, will trigger calculation
                        print(f"Debug - Invalid override value for {col}: '{value}', will calculate")
        
        return overrides
        
    except Exception as e:
        print(f"Debug - Error loading boundary overrides: {e}")
        return None


def filter_empty_dlc_rows(dlc_df):
    """
    Filter out rows where all tracked points have zero coordinates (empty frames)
    
    Parameters:
    - dlc_df: DataFrame with DLC data
    
    Returns:
    - Filtered DataFrame with empty rows removed
    """
    # Get all x and y columns (assuming they follow pattern x, y, x.1, y.1, etc.)
    x_columns = [col for col in dlc_df.columns if col.startswith('x')]
    y_columns = [col for col in dlc_df.columns if col.startswith('y')]
    
    # Combine all position columns
    position_columns = x_columns + y_columns
    
    # Find rows where ALL position values are 0
    # This indicates frames with no tracking data
    empty_rows = (dlc_df[position_columns] == 0).all(axis=1)
    
    # Count and report empty rows
    num_empty = empty_rows.sum()
    total_rows = len(dlc_df)
    
    if num_empty > 0:
        print(f"Found {num_empty} empty rows out of {total_rows} total rows ({num_empty/total_rows*100:.1f}%)")
        
        # Filter out empty rows
        dlc_df_filtered = dlc_df[~empty_rows].copy()
        
        # Reset index to ensure continuous frame numbers
        dlc_df_filtered.reset_index(drop=True, inplace=True)
        
        print(f"Filtered data contains {len(dlc_df_filtered)} valid frames")
        
        return dlc_df_filtered
    else:
        print(f"No empty rows found in {total_rows} total rows")
        return dlc_df


def calculate_min_x_boundary(x, left_x, right_x, threshold_percent):
    """Calculate min_x_dlc using density-based approach"""
    x_in_chamber = x[(x >= left_x) & (x <= right_x)]
    
    if len(x_in_chamber) == 0:
        print("Warning: No data found within chamber X boundaries!")
        return left_x
    
    x_range = np.max(x_in_chamber) - np.min(x_in_chamber)
    n_bins_x = max(100, int(x_range / 5))
    x_hist, x_edges = np.histogram(x_in_chamber, bins=n_bins_x)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    
    x_density = gaussian_filter1d(x_hist.astype(float), sigma=1)
    density_threshold = np.max(x_density) * (threshold_percent / 100.0)
    significant_bins = np.where(x_density > density_threshold)[0]
    
    if len(significant_bins) > 0:
        min_x_dlc = x_centers[significant_bins[0]]
        min_x_dlc = max(min_x_dlc, left_x)  # Ensure within chamber limits
        print(f"Debug - Calculated min_x_dlc using {threshold_percent}% density: {min_x_dlc:.1f}")
        return min_x_dlc
    else:
        print(f"Debug - No significant X density found, using left boundary: {left_x}")
        return left_x


def calculate_max_x_boundary(x, left_x, right_x, threshold_percent, min_x_dlc):
    """Calculate max_x_dlc using density-based approach"""
    x_in_chamber = x[(x >= left_x) & (x <= right_x)]
    
    if len(x_in_chamber) == 0:
        return right_x
    
    x_range = np.max(x_in_chamber) - np.min(x_in_chamber)
    n_bins_x = max(100, int(x_range / 5))
    x_hist, x_edges = np.histogram(x_in_chamber, bins=n_bins_x)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    
    x_density = gaussian_filter1d(x_hist.astype(float), sigma=1)
    density_threshold = np.max(x_density) * (threshold_percent / 100.0)
    significant_bins = np.where(x_density > density_threshold)[0]
    
    if len(significant_bins) > 0:
        max_x_dlc = x_centers[significant_bins[-1]]
        max_x_dlc = min(max_x_dlc, right_x)  # Ensure within chamber limits
        print(f"Debug - Calculated max_x_dlc using {threshold_percent}% density: {max_x_dlc:.1f}")
        return max_x_dlc
    else:
        print(f"Debug - No significant X density found, using right boundary: {right_x}")
        return right_x


def calculate_min_y_boundary(y, low_y):
    """Calculate min_y_dlc using constraint-based approach"""
    y_range = np.max(y) - np.min(y)
    n_bins_y = max(100, int(y_range / 5))
    y_hist, y_edges = np.histogram(y, bins=n_bins_y)
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    
    y_density = gaussian_filter1d(y_hist.astype(float), sigma=1)
    min_y_threshold = np.max(y_density) * 0.005  # 0.5% threshold
    min_y_candidates = np.where(y_density > min_y_threshold)[0]
    
    if len(min_y_candidates) > 0:
        min_y_density_based = y_centers[min_y_candidates[0]]
        raw_min_y = np.min(y)
        
        if min_y_density_based > low_y:
            min_y_dlc = raw_min_y + (low_y - raw_min_y) * 0.3
            print(f"Debug - Calculated min_y_dlc (interpolated): {min_y_dlc:.1f}")
        elif min_y_density_based < (raw_min_y + 10):
            min_y_dlc = raw_min_y + (low_y - raw_min_y) * 0.05
            print(f"Debug - Calculated min_y_dlc (adjusted): {min_y_dlc:.1f}")
        else:
            min_y_dlc = min_y_density_based
            print(f"Debug - Calculated min_y_dlc (density-based): {min_y_dlc:.1f}")
        return min_y_dlc
    else:
        raw_min_y = np.min(y)
        min_y_dlc = raw_min_y + (low_y - raw_min_y) * 0.25
        print(f"Debug - Calculated min_y_dlc (fallback): {min_y_dlc:.1f}")
        return min_y_dlc


def calculate_max_y_boundary(y, high_y):
    """Calculate max_y_dlc using density-based approach"""
    y_range = np.max(y) - np.min(y)
    n_bins_y = max(100, int(y_range / 5))
    y_hist, y_edges = np.histogram(y, bins=n_bins_y)
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    
    y_density = gaussian_filter1d(y_hist.astype(float), sigma=1)
    max_y_threshold = np.max(y_density) * 0.01  # 1% threshold
    max_y_candidates = np.where(y_density > max_y_threshold)[0]
    
    if len(max_y_candidates) > 0:
        max_y_density_based = y_centers[max_y_candidates[-1]]
        raw_max_y = np.max(y)
        
        if max_y_density_based < high_y:
            max_y_dlc = high_y + (raw_max_y - high_y) * 0.3
            print(f"Debug - Calculated max_y_dlc (interpolated): {max_y_dlc:.1f}")
        elif max_y_density_based > (raw_max_y - 5):
            max_y_dlc = raw_max_y - 5
            print(f"Debug - Calculated max_y_dlc (adjusted): {max_y_dlc:.1f}")
        else:
            max_y_dlc = max_y_density_based
            print(f"Debug - Calculated max_y_dlc (density-based): {max_y_dlc:.1f}")
        return max_y_dlc
    else:
        raw_max_y = np.max(y)
        max_y_dlc = high_y + (raw_max_y - high_y) * 0.85
        print(f"Debug - Calculated max_y_dlc (fallback): {max_y_dlc:.1f}")
        return max_y_dlc


def calculate_velocity(x, y, fps):
    """Calculate velocity in pixels per second based on head position"""
    # Calculate displacement between consecutive frames
    dx = np.diff(x)
    dy = np.diff(y)
    displacement = np.sqrt(dx**2 + dy**2)
    
    # Calculate velocity in pixels per second
    velocity = displacement * fps
    
    # Handle NaN or Inf values
    velocity = np.nan_to_num(velocity, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Add a 0 at the beginning to make velocity array same length as position arrays
    velocity = np.insert(velocity, 0, 0)
    
    return velocity


def detect_gate_time(x, y, boundaries, gate_radius=30):
    """
    Detect time spent at chamber gates
    
    Parameters:
    - x, y: position arrays
    - boundaries: dict with high_y, low_y, mid_x
    - gate_radius: radius of circle defining the gate area in pixels
    
    Returns:
    - mask arrays for top and bottom gates
    """
    # Define gate centers
    top_gate_center = (boundaries['mid_x'], boundaries['low_y'])
    bottom_gate_center = (boundaries['mid_x'], boundaries['high_y'])
    
    # Calculate distance from each point to gate centers
    top_distances = np.sqrt((x - top_gate_center[0])**2 + (y - top_gate_center[1])**2)
    bottom_distances = np.sqrt((x - bottom_gate_center[0])**2 + (y - bottom_gate_center[1])**2)
    
    # Create masks for points within gate radius
    top_gate_mask = top_distances <= gate_radius
    bottom_gate_mask = bottom_distances <= gate_radius
    
    return {
        'top_gate_mask': top_gate_mask,
        'bottom_gate_mask': bottom_gate_mask,
        'top_gate_center': top_gate_center,
        'bottom_gate_center': bottom_gate_center,
        'gate_radius': gate_radius
    }


def detect_corner_time(x, y, boundaries, corner_radius=60, experiment_name=None):
    """
    Detect time spent in corners using quarter circles at chamber corners
    
    Parameters:
    - x, y: position arrays (already cleaned)
    - boundaries: dict with high_y, low_y, left_x, right_x
    - corner_radius: radius of quarter circle defining the corner area in pixels
    - experiment_name: name of the experiment for boundary override lookup
    
    Returns:
    - mask arrays for all corners
    """
    
    # Get boundary constraints
    left_x = boundaries['left_x']
    right_x = boundaries['right_x']
    low_y = boundaries['low_y']
    high_y = boundaries['high_y']
    
    # Load boundary overrides if experiment_name is provided
    overrides = load_boundary_overrides(experiment_name) if experiment_name else None
    
    # Simple density-based boundary detection
    # Find the actual borders of DLC data based on density of datapoints
    # Boundaries must be within chamber constraints (left_x to right_x)
    
    # Adjustable threshold parameter (can be modified easily)
    DENSITY_THRESHOLD_PERCENT = 1.0  # 1% of max density - adjust this value as needed
    
    # Initialize boundary values
    min_x_dlc = None
    max_x_dlc = None
    min_y_dlc = None
    max_y_dlc = None
    
    # Check for overrides and calculate missing values
    if overrides:
        min_x_dlc = overrides.get('min_x_dlc')
        max_x_dlc = overrides.get('max_x_dlc')
        min_y_dlc = overrides.get('min_y_dlc')
        max_y_dlc = overrides.get('max_y_dlc')
    
    # Calculate min_x_dlc if not overridden or invalid
    if min_x_dlc is None:
        min_x_dlc = calculate_min_x_boundary(x, left_x, right_x, DENSITY_THRESHOLD_PERCENT)
    
    # Calculate max_x_dlc if not overridden or invalid
    if max_x_dlc is None:
        max_x_dlc = calculate_max_x_boundary(x, left_x, right_x, DENSITY_THRESHOLD_PERCENT, min_x_dlc)
    
    # Calculate min_y_dlc if not overridden or invalid
    if min_y_dlc is None:
        min_y_dlc = calculate_min_y_boundary(y, low_y)
    
    # Calculate max_y_dlc if not overridden or invalid
    if max_y_dlc is None:
        max_y_dlc = calculate_max_y_boundary(y, high_y)
    
    print(f"DLC data range (density-based) - X: {min_x_dlc:.1f} to {max_x_dlc:.1f}, Y: {min_y_dlc:.1f} to {max_y_dlc:.1f}")
    print(f"Chamber boundaries - X: {left_x} to {right_x}, Y: {low_y} to {high_y}")
    
    # Store the final boundary values that will be used (after any corrections)
    # Update Y values as well
    final_min_y = min_y_dlc
    final_max_y = max_y_dlc
    
    # Define corners using quarter circles at DLC data extremes
    # All 12 corners (4 per chamber) with proper coordinates
    # Remember: min_y_dlc = TOP, max_y_dlc = BOTTOM, low_y = TOP boundary, high_y = BOTTOM boundary
    corners = []
    
    # Chamber 1 (Top chamber) - 4 corners
    # Top-left corner - actual TOP edge of data
    corners.append({
        'center': (min_x_dlc, min_y_dlc),  # min_y_dlc is the TOP
        'type': 'top_left',
        'name': 'Chamber1_TopLeft',
        'chamber': 1
    })
    
    # Top-right corner - actual TOP edge of data
    corners.append({
        'center': (max_x_dlc, min_y_dlc),  # min_y_dlc is the TOP
        'type': 'top_right',
        'name': 'Chamber1_TopRight',
        'chamber': 1
    })
    
    # Bottom-left corner - at TOP chamber boundary (low_y)
    corners.append({
        'center': (min_x_dlc, low_y),  # low_y is the TOP boundary
        'type': 'bottom_left',
        'name': 'Chamber1_BottomLeft',
        'chamber': 1
    })
    
    # Bottom-right corner - at TOP chamber boundary (low_y)
    corners.append({
        'center': (max_x_dlc, low_y),  # low_y is the TOP boundary
        'type': 'bottom_right',
        'name': 'Chamber1_BottomRight',
        'chamber': 1
    })
    
    # Chamber 2 (Middle chamber) - 4 corners
    # Top-left corner - slightly below TOP boundary
    corners.append({
        'center': (min_x_dlc, low_y + 5),  # 5px below TOP boundary
        'type': 'top_left',
        'name': 'Chamber2_TopLeft',
        'chamber': 2
    })
    
    # Top-right corner - slightly below TOP boundary
    corners.append({
        'center': (max_x_dlc, low_y + 5),  # 5px below TOP boundary
        'type': 'top_right',
        'name': 'Chamber2_TopRight',
        'chamber': 2
    })
    
    # Bottom-left corner - slightly above BOTTOM boundary
    corners.append({
        'center': (min_x_dlc, high_y - 5),  # 5px above BOTTOM boundary
        'type': 'bottom_left',
        'name': 'Chamber2_BottomLeft',
        'chamber': 2
    })
    
    # Bottom-right corner - slightly above BOTTOM boundary
    corners.append({
        'center': (max_x_dlc, high_y - 5),  # 5px above BOTTOM boundary
        'type': 'bottom_right',
        'name': 'Chamber2_BottomRight',
        'chamber': 2
    })
    
    # Chamber 3 (Bottom chamber) - 4 corners
    # Top-left corner - at BOTTOM chamber boundary (high_y)
    corners.append({
        'center': (min_x_dlc, high_y),  # high_y is the BOTTOM boundary
        'type': 'top_left',
        'name': 'Chamber3_TopLeft',
        'chamber': 3
    })
    
    # Top-right corner - at BOTTOM chamber boundary (high_y)
    corners.append({
        'center': (max_x_dlc, high_y),  # high_y is the BOTTOM boundary
        'type': 'top_right',
        'name': 'Chamber3_TopRight',
        'chamber': 3
    })
    
    # Bottom-left corner - actual BOTTOM edge of data
    corners.append({
        'center': (min_x_dlc, max_y_dlc),  # max_y_dlc is the BOTTOM
        'type': 'bottom_left',
        'name': 'Chamber3_BottomLeft',
        'chamber': 3
    })
    
    # Bottom-right corner - actual BOTTOM edge of data
    corners.append({
        'center': (max_x_dlc, max_y_dlc),  # max_y_dlc is the BOTTOM
        'type': 'bottom_right',
        'name': 'Chamber3_BottomRight',
        'chamber': 3
    })
    
    # Calculate masks for each corner
    corner_masks = []
    for corner in corners:
        corner_x, corner_y = corner['center']
        
        # Calculate distance from each point to corner center
        distances = np.sqrt((x - corner_x)**2 + (y - corner_y)**2)
        
        # Create mask for points within corner radius
        within_radius = distances <= corner_radius
        
        # For quarter circles, check quadrant based on corner type
        if corner['type'] == 'top_left':
            # Top-left quarter: x >= center_x AND y >= center_y
            in_quadrant = (x >= corner_x) & (y >= corner_y)
        elif corner['type'] == 'top_right':
            # Top-right quarter: x <= center_x AND y >= center_y
            in_quadrant = (x <= corner_x) & (y >= corner_y)
        elif corner['type'] == 'bottom_left':
            # Bottom-left quarter: x >= center_x AND y <= center_y
            in_quadrant = (x >= corner_x) & (y <= corner_y)
        else:  # bottom_right
            # Bottom-right quarter: x <= center_x AND y <= center_y
            in_quadrant = (x <= corner_x) & (y <= corner_y)
        
        # Combine conditions
        corner_mask = within_radius & in_quadrant
        
        corner_masks.append({
            'mask': corner_mask,
            'center': (corner_x, corner_y),
            'radius': corner_radius,
            'name': corner['name'],
            'type': corner['type'],
            'chamber': corner['chamber']
        })
        
        print(f"{corner['name']}: center=({corner_x:.1f}, {corner_y:.1f}), "
              f"points_in_corner={np.sum(corner_mask)}")
    
    return corner_masks


def analyze_mouse_movement(dlc_df, boundaries, fps, experiment_name=None):
    """
    Analyze mouse movement including velocity, gate time, and corner time
    
    Parameters:
    - dlc_df: DataFrame with DLC data
    - boundaries: dict with boundary information
    - fps: frames per second
    - experiment_name: name of the experiment for boundary override lookup
    
    Returns:
    - dict with analysis results
    """
    # FIRST: Filter out empty rows (rows with all zeros)
    dlc_df_filtered = filter_empty_dlc_rows(dlc_df)
    
    # Extract head position from filtered data
    x_raw = dlc_df_filtered['x.1'].values
    y_raw = dlc_df_filtered['y.1'].values
    
    # Clean the data - remove obvious tracking errors
    # Method 1: Remove negative values (though after filtering zeros, this might be less common)
    valid_mask = (x_raw >= 0) & (y_raw >= 0)
    
    # Method 2: Remove NaN and infinite values
    valid_mask = valid_mask & np.isfinite(x_raw) & np.isfinite(y_raw)
    
    # Apply cleaning
    x = x_raw.copy()
    y = y_raw.copy()
    
    # For invalid points, interpolate if possible
    if np.sum(~valid_mask) > 0 and np.sum(valid_mask) > 10:
        valid_indices = np.where(valid_mask)[0]
        invalid_indices = np.where(~valid_mask)[0]
        
        # Only interpolate for small gaps (less than 10 frames)
        for idx in invalid_indices:
            # Find nearest valid points
            before_idx = valid_indices[valid_indices < idx]
            after_idx = valid_indices[valid_indices > idx]
            
            if len(before_idx) > 0 and len(after_idx) > 0:
                nearest_before = before_idx[-1]
                nearest_after = after_idx[0]
                
                # Only interpolate for small gaps
                if (nearest_after - nearest_before) < 10:
                    # Linear interpolation
                    alpha = (idx - nearest_before) / (nearest_after - nearest_before)
                    x[idx] = x[nearest_before] * (1 - alpha) + x[nearest_after] * alpha
                    y[idx] = y[nearest_before] * (1 - alpha) + y[nearest_after] * alpha
                else:
                    # For large gaps, use nearest valid value
                    if (idx - nearest_before) < (nearest_after - idx):
                        x[idx] = x[nearest_before]
                        y[idx] = y[nearest_before]
                    else:
                        x[idx] = x[nearest_after]
                        y[idx] = y[nearest_after]
            elif len(before_idx) > 0:
                x[idx] = x[before_idx[-1]]
                y[idx] = y[before_idx[-1]]
            elif len(after_idx) > 0:
                x[idx] = x[after_idx[0]]
                y[idx] = y[after_idx[0]]
    
    # Final check - ensure no negative values remain
    x = np.maximum(x, 0)
    y = np.maximum(y, 0)
    
    print(f"Data cleaning summary: {np.sum(~valid_mask)} points cleaned/interpolated out of {len(x)} total valid points")
    
    # Calculate velocity
    velocity = calculate_velocity(x, y, fps)
    avg_velocity = np.mean(velocity)
    
    # Gate time analysis (with 30px radius)
    gate_results = detect_gate_time(x, y, boundaries, gate_radius=30)
    top_gate_time = np.sum(gate_results['top_gate_mask']) / fps
    bottom_gate_time = np.sum(gate_results['bottom_gate_mask']) / fps
    
    # Corner time analysis
    corner_results = detect_corner_time(x, y, boundaries, corner_radius=70, experiment_name=experiment_name)
    corner_times = []
    for i, corner in enumerate(corner_results):
        corner_time = np.sum(corner['mask']) / fps
        corner_times.append(corner_time)
    
    # Total recording time (based on filtered data)
    total_time = len(x) / fps
    
    # Define final boundary values for visualization consistency
    # These represent the actual bounds used in corner detection after any corrections
    raw_min_x = np.min(x)
    raw_max_x = np.max(x)
    raw_min_y = np.min(y)
    raw_max_y = np.max(y)
    
    # Extract the boundary values that were actually used in corner detection
    # by getting the min/max from the corner centers
    if corner_results and len(corner_results) > 0:
        corner_x_positions = [corner['center'][0] for corner in corner_results]
        corner_y_positions = [corner['center'][1] for corner in corner_results]
        final_min_x = min(corner_x_positions)
        final_max_x = max(corner_x_positions)
        final_min_y = min(corner_y_positions)
        final_max_y = max(corner_y_positions)
    else:
        # Fallback to raw data if no corners detected
        final_min_x = raw_min_x
        final_max_x = raw_max_x
        final_min_y = raw_min_y
        final_max_y = raw_max_y
    
    # Calculate number of rows removed
    num_rows_removed = len(dlc_df) - len(dlc_df_filtered)
    
    return {
        'x': x,
        'y': y,
        'velocity': velocity,
        'avg_velocity': avg_velocity,
        'gate_results': gate_results,
        'top_gate_time': top_gate_time,
        'bottom_gate_time': bottom_gate_time,
        'corner_results': corner_results,
        'corner_times': corner_times,
        'total_time': total_time,
        'fps': fps,
        'num_cleaned_points': np.sum(~valid_mask),
        'percent_cleaned': (np.sum(~valid_mask) / len(x)) * 100,
        'num_empty_rows_removed': num_rows_removed,
        'original_frame_count': len(dlc_df),
        'valid_frame_count': len(dlc_df_filtered),
        # Add the actual boundary values used for corner detection
        'actual_bounds': {
            'min_x': final_min_x,
            'max_x': final_max_x, 
            'min_y': final_min_y,
            'max_y': final_max_y,
            'raw_min_x': raw_min_x,
            'raw_max_x': raw_max_x,
            'raw_min_y': raw_min_y,
            'raw_max_y': raw_max_y
        }
    }