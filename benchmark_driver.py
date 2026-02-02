import json
import subprocess
import argparse
import sys
import os
import time

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

        walltime, verification = map(float, last_line.split(','))
        return walltime, verification
    except subprocess.CalledProcessError as e:
        print(f"Error running configuration {config['name']}: {e}")
        return None, None
    except ValueError:
        print(f"Error parsing output for {config['name']}: {result.stdout}")
        return None, None

def print_table(results):
    reference_time = results[0]['time'] if results else 0
    reference_x = results[0]['verify'] if results else 0
    
    print("\n" + "="*85)
    print(f"{'Method/Configuration':<30} | {'Time (s)':<10} | {'Speedup':<10} | {'Rel Error':<15}")
    print("-" * 85)
    
    for res in results:
        name = res['name']
        t = res['time']
        
        if t is None:
            print(f"{name:<30} | {'FAILED':<10} | {'-':<10} | {'-':<15}")
            continue
            
        speedup = reference_time / t if t > 0 else 0
        
        # Relative error vs reference
        rel_error = 0.0
        if reference_x != 0:
            rel_error = abs(res['verify'] - reference_x) / abs(reference_x)
            
        print(f"{name:<30} | {t:<10.4f} | {speedup:<10.2f}x | {rel_error:<15.2e}")
    print("=" * 85 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Run REBOUND Ablation Benchmarks')
    parser.add_argument('--config', default='benchmarks.json', help='Path to configuration JSON')
    parser.add_argument('--profile', action='store_true', help='Enable profiling (-DPROF)')
    parser.add_argument('--filter', help='Run only configs containing this string')
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
        t, v = run_benchmark(cfg)
        if t is not None:
             print(f" Done ({t:.4f}s)")
        else:
             print(f" FAILED")
             
        results.append({
            'name': cfg['name'],
            'time': t,
            'verify': v
        })
        
    print_table(results)

if __name__ == "__main__":
    main()
