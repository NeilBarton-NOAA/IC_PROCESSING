#!/usr/bin/env python3
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import xarray as xr


def process_member_recenter(args):
    """Worker function to re-center a single ensemble member."""
    mem_num, origin_dir, work_dir, pert_filename, mean_nc = args
    t0_mem = time.perf_counter()

    mem_str = f"mem{mem_num:03d}"
    in_nc = os.path.join(origin_dir, mem_str, "analysis", "ocean", pert_filename)
    
    out_dir = os.path.join(work_dir, mem_str, "analysis", "ocean")
    os.makedirs(out_dir, exist_ok=True)
    out_nc = os.path.join(out_dir, pert_filename)

    if not os.path.exists(in_nc):
        return False, f"[ERROR] Input member file missing: {in_nc}", 0.0, 0.0

    try:
        ds_mean = xr.open_dataset(mean_nc)

        with xr.open_dataset(in_nc) as ds_mem:
            ds_recentered = ds_mem.copy(deep=True)

            # Re-center t_pert and s_pert (Subtract Mean)
            if "t_pert" in ds_recentered and "t_pert" in ds_mean:
                ds_recentered["t_pert"] = (ds_mem["t_pert"] - ds_mean["t_pert"]).astype(np.float32)
            
            if "s_pert" in ds_recentered and "s_pert" in ds_mean:
                ds_recentered["s_pert"] = (ds_mem["s_pert"] - ds_mean["s_pert"]).astype(np.float32)

            for var in ds_recentered.data_vars:
                ds_recentered[var].attrs = ds_mem[var].attrs
            ds_recentered.attrs = ds_mem.attrs

            encoding = {var: {"zlib": True, "complevel": 1} for var in ds_recentered.data_vars}
            ds_recentered.to_netcdf(out_nc, encoding=encoding)

        ds_mean.close()
        elapsed = time.perf_counter() - t0_mem
        out_size_mb = os.path.getsize(out_nc) / (1024 * 1024)

        return True, f"✓ Re-centered {mem_str} in {elapsed:5.2f}s | Output: {out_size_mb:6.1f} MB", elapsed, out_size_mb

    except Exception as e:
        return False, f"[ERROR] Failed processing {mem_str}: {str(e)}", 0.0, 0.0


def main():
    total_start = time.perf_counter()

    if len(sys.argv) < 7:
        print(
            "[FATAL ERROR] Usage: python recenter_parallel.py "
            "<origin_dir> <work_dir> <mem_start> <mem_end> <pert_filename> <mean_nc> <num_workers>"
        )
        sys.exit(1)

    origin_dir = sys.argv[1]
    work_dir = sys.argv[2]
    mem_start = int(sys.argv[3])
    mem_end = int(sys.argv[4])
    pert_filename = sys.argv[5]
    mean_nc = sys.argv[6]
    num_workers = int(sys.argv[7]) if len(sys.argv) > 7 else (mem_end - mem_start + 1)

    n_mems = mem_end - mem_start + 1

    print("=" * 75)
    print("STEP 2: PARALLEL RE-CENTERING (SUBTRACT ENSEMBLE MEAN)")
    print("=" * 75)
    print(f" Source Directory : {origin_dir}")
    print(f" Work Directory   : {work_dir}")
    print(f" Target Members   : mem{mem_start:03d} to mem{mem_end:03d} ({n_mems} members)")
    print(f" Ensemble Mean NC : {mean_nc}")
    print(f" Parallel Workers : {num_workers}")
    print("-" * 75)

    if not os.path.exists(mean_nc):
        print(f"[FATAL ERROR] Ensemble mean file not found: {mean_nc}")
        sys.exit(1)

    tasks = [
        (m, origin_dir, work_dir, pert_filename, mean_nc)
        for m in range(mem_start, mem_end + 1)
    ]

    print(f"Launching {n_mems} parallel tasks across {num_workers} worker(s)...\n")
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = list(executor.map(process_member_recenter, tasks))

    failures = 0
    total_mb = 0.0
    for success, msg, _, size_mb in results:
        print(msg)
        total_mb += size_mb
        if not success:
            failures += 1

    total_elapsed = time.perf_counter() - total_start

    print("-" * 75)
    if failures > 0:
        print(f"[FATAL ERROR] Step 2 failed for {failures} member(s)!")
        sys.exit(1)
    else:
        print(f"✓ RE-CENTERING COMPLETED SUCCESSFULLY!")
        print(f"  Processed {n_mems} members ({total_mb:.1f} MB total) in {total_elapsed:.2f} seconds.")
        print("=" * 75)


if __name__ == "__main__":
    main()
