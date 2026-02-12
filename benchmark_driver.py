import json
import subprocess
import argparse
import sys
import os
import time
import csv
from datetime import datetime

def load_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

def compile_rebound(profile=False):
    print(f"Compiling REBOUND (Profile mode: {profile})...")
    opt_flags = "-march=native -O3"
    if profile:
        opt_flags += " -DPROF"
    
    os.environ['OPT'] = opt_flags
    
    try:
        subprocess.check_call(['make', 'clean'], cwd='examples/whfast512_solar_system_jac', stdout=subprocess.DEVNULL)
        subprocess.check_call(['make'], cwd='examples/whfast512_solar_system_jac', stdout=subprocess.DEVNULL)
        print("Compilation successful.")
        return True
    except subprocess.CalledProcessError:
        print("Compilation failed!")
        return False

def run_benchmark(config, executable='./rebound'):
    cmd = [executable]
    for key, value in config['args'].items():
        cmd.append(f"--{key}")
        cmd.append(str(value))
    
    try:
        result = subprocess.run(
            cmd, 
            cwd='examples/whfast512_solar_system_jac',
            capture_output=True,
            text=True,
            check=True
        )
        
        if "PROFILING_START" in result.stdout:
            print("\n" + "-"*40)
            print(f"Profiling Output for {config['name']}:")
            start = result.stdout.find("PROFILING_START")
            end = result.stdout.find("PROFILING_END") + 13
            print(result.stdout[start:end])
            print("-"*40 + "\n")
            
        lines = result.stdout.strip().split('\n')
        last_line = lines[-1]
        
        if "PROFILING" in last_line: 
             pass

        parts = last_line.split(',')
        walltime = float(parts[0])
        verification = float(parts[1])
        energy_err = float(parts[2]) if len(parts) > 2 else None
        return walltime, verification, energy_err
    except subprocess.CalledProcessError as e:
        print(f"Error running configuration {config['name']}: {e}")
        return None, None
    except ValueError:
        print(f"Error parsing output for {config['name']}: {result.stdout}")
        return None, None

def write_csv(results, output_file, profile_mode):
    """Write results to CSV file"""
    # Use WHFast512 Jacobi (Baseline) as reference for speedup
    baseline = next((r for r in results if 'Jacobi (Baseline)' in r['name']), results[0])
    reference_time = baseline['time'] if baseline['time'] is not None else 0
    reference_x = results[0]['verify'] if results else 0
    
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['configuration', 'time_s', 'speedup', 'position_error', 'energy_error', 'profile_mode', 'timestamp']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        
        for res in results:
            name = res['name']
            t = res['time']
            
            if t is None:
                writer.writerow({
                    'configuration': name,
                    'time_s': 'FAILED',
                    'speedup': 'N/A',
                    'position_error': 'N/A',
                    'energy_error': 'N/A',
                    'profile_mode': profile_mode,
                    'timestamp': timestamp
                })
                continue
            
            speedup = reference_time / t if t > 0 else 0
            
            # Position Error vs reference
            pos_error = 0.0
            if reference_x != 0:
                pos_error = abs(res['verify'] - reference_x) / abs(reference_x)
                
            # Energy Error (from simulation)
            energy_error = res.get('energy_err', 0.0)
            if energy_error is None: energy_error = 0.0
            
            writer.writerow({
                'configuration': name,
                'time_s': f"{t:.6f}",
                'speedup': f"{speedup:.4f}",
                'position_error': f"{pos_error:.6e}",
                'energy_error': f"{energy_error:.6e}",
                'profile_mode': profile_mode,
                'timestamp': timestamp
            })
    
    print(f"Results written to: {output_file}")

def print_table(results):
    # Use WHFast512 Jacobi (Baseline) as reference for speedup
    baseline = next((r for r in results if 'Jacobi (Baseline)' in r['name']), results[0])
    reference_time = baseline['time'] if baseline['time'] is not None else 0
    reference_x = results[0]['verify'] if results else 0
    
    print(f"\n" + "="*100)
    print(f"{'Method/Configuration':<30} | {'Time (s)':<10} | {'Speedup':<10} | {'Pos Error':<12} | {'Energy Err':<12}")
    print("-" * 100)
    
    for res in results:
        name = res['name']
        t = res['time']
        
        if t is None:
            print(f"{name:<30} | {'FAILED':<10} | {'-':<10} | {'-':<12} | {'-':<12}")
            continue
            
        speedup = reference_time / t if t > 0 else 0
        
        # Position Error vs reference
        pos_error = 0.0
        if reference_x != 0:
            pos_error = abs(res['verify'] - reference_x) / abs(reference_x)
            
        # Energy Error (from simulation)
        energy_error = res.get('energy_err', 0.0)
        if energy_error is None: energy_error = 0.0

        print(f"{name:<30} | {t:<10.4f} | {speedup:<10.2f}x | {pos_error:<12.2e} | {energy_error:<12.2e}")
    print("=" * 100 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Run REBOUND Ablation Benchmarks')
    parser.add_argument('--config', default='benchmarks.json', help='Path to configuration JSON')
    parser.add_argument('--profile', action='store_true', help='Enable profiling (-DPROF)')
    parser.add_argument('--filter', help='Run only configs containing this string')
    parser.add_argument('--output', default=None, help='Custom output CSV filename')
    args = parser.parse_args()
    
    configs = load_config(args.config)
    
    if args.filter:
        configs = [c for c in configs if args.filter in c['name']]
    
    if not compile_rebound(args.profile):
        sys.exit(1)
        
    results = []
    print(f"\nRunning {len(configs)} configurations...\n")
    
    for i, cfg in enumerate(configs):
        print(f"[{i+1}/{len(configs)}] Running: {cfg['name']}...", end='', flush=True)
        t, v, e = run_benchmark(cfg)
        if t is not None:
             print(f" Done ({t:.4f}s)")
        else:
             print(f" FAILED")
             
        results.append({
            'name': cfg['name'],
            'time': t,
            'verify': v,
            'energy_err': e
        })
        
    print_table(results)
    
    # Determine output filename
    if args.output:
        output_file = args.output
    else:
        profile_suffix = "_profiled" if args.profile else "_noprofile"
        timestamp_suffix = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_file = f"benchmark_results{profile_suffix}_{timestamp_suffix}.csv"
    
    # Write CSV output
    profile_mode = "profiled" if args.profile else "no_profile"
    write_csv(results, output_file, profile_mode)

if __name__ == "__main__":
    main()
