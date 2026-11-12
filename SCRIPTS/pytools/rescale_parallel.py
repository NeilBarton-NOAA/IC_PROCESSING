#!/usr/bin/env python3
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import xarray as xr


def process_member_rescale(args):
    """Worker function to scale a single recentered ensemble member using depth interpolation."""
    mem_num, work_dir, out_root, yyyymmdd, hh, pert_filename, scale_csv, max_cap = args
    t0_mem = time.perf_counter()

    mem_str = f"mem{mem_num:03d}"
    
    in_nc = os.path.join(work_dir, mem_str, "analysis", "ocean", pert_filename)
    out_dir = os.path.join(out_root, f"sfs.{yyyymmdd}", hh, mem_str, "analysis", "ocean")
    os.makedirs(out_dir, exist_ok=True)
    out_nc = os.path.join(out_dir, pert_filename)

    if not os.path.exists(in_nc):
        return False, f"[ERROR] Recentered input file missing: {in_nc}", 0.0, 0.0

    try:
        # 1. Load CSV scaling profile (32 Z-levels)
        df_scale = pd.read_csv(scale_csv)
        df_scale.columns = df_scale.columns.str.strip()

        depth_col = [c for c in df_scale.columns if 'depth' in c.lower()][0]
        t_col = [c for c in df_scale.columns if 'ratio_t' in c.lower() or 'temp' in c.lower()][0]
        s_col = [c for c in df_scale.columns if 'ratio_s' in c.lower() or 'sal' in c.lower() or 'so' in c.lower()][0]

        csv_depths = df_scale[depth_col].values.astype(np.float64)
        csv_t_ratios = df_scale[t_col].values.astype(np.float64)
        csv_s_ratios = df_scale[s_col].values.astype(np.float64)

        # 2. Apply upper bound cap across entire column to prevent numerical instability
        csv_t_ratios = np.minimum(csv_t_ratios, max_cap)
        csv_s_ratios = np.minimum(csv_s_ratios, max_cap)

        with xr.open_dataset(in_nc) as ds:
            ds_scaled = ds.copy(deep=True)

            # Identify vertical dimension/coordinate
            z_dim = None
            for dim in ["Layer", "z_l", "z_i"]:
                if dim in ds.dims or dim in ds.coords:
                    z_dim = dim
                    break

            if z_dim is None:
                return False, f"[ERROR] Could not identify vertical dimension in {in_nc}", 0.0, 0.0

            # Extract target depths (75 layers) from NetCDF 'Layer' variable
            if z_dim in ds.coords:
                target_depths = ds[z_dim].values.astype(np.float64)
            else:
                return False, f"[ERROR] Vertical coordinate '{z_dim}' missing depth values in {in_nc}", 0.0, 0.0

            # 3. Perform 1D Depth Interpolation
            # left=csv_ratios[0] holds shallowest value, right=csv_ratios[-1] holds constant for deep layers below 467.4m
            t_factors = np.interp(
                target_depths,
                csv_depths,
                csv_t_ratios,
                left=csv_t_ratios[0],
                right=csv_t_ratios[-1]
            ).astype(np.float32)

            s_factors = np.interp(
                target_depths,
                csv_depths,
                csv_s_ratios,
                left=csv_s_ratios[0],
                right=csv_s_ratios[-1]
            ).astype(np.float32)

            # Construct 1D DataArray for 3D spatial broadcasting across (Layer, lath, lonh)
            t_scale_da = xr.DataArray(t_factors, dims=[z_dim])
            s_scale_da = xr.DataArray(s_factors, dims=[z_dim])

            # Apply vertical scaling to perturbations
            if "t_pert" in ds_scaled:
                ds_scaled["t_pert"] = (ds_scaled["t_pert"] * t_scale_da).astype(np.float32)
            if "s_pert" in ds_scaled:
                ds_scaled["s_pert"] = (ds_scaled["s_pert"] * s_scale_da).astype(np.float32)

            # Preserve metadata attributes
            for var in ds_scaled.data_vars:
                ds_scaled[var].attrs = ds[var].attrs
            ds_scaled.attrs = ds.attrs

            # Set zlib compression encoding
            encoding = {var: {"zlib": True, "complevel": 1} for var in ds_scaled.data_vars}
            ds_scaled.to_netcdf(out_nc, encoding=encoding)

        elapsed = time.perf_counter() - t0_mem
        out_size_mb = os.path.getsize(out_nc) / (1024 * 1024)

        return True, f"✓ Rescaled {mem_str} in {elapsed:5.2f}s | Output: {out_size_mb:6.1f} MB", elapsed, out_size_mb

    except Exception as e:
        return False, f"[ERROR] Failed scaling {mem_str}: {str(e)}", 0.0, 0.0


def main():
    total_start = time.perf_counter()

    # Allow 9 or 10 command-line arguments safely
    if len(sys.argv) < 10:
        print(
            "[FATAL ERROR] Usage: python rescale_parallel.py "
            "<work_dir> <out_root> <yyyymmdd> <hh> <mem_start> <mem_end> <pert_filename> <scale_csv> <max_cap> [num_workers]"
        )
        sys.exit(1)

    work_dir = sys.argv[1]
    out_root = sys.argv[2]
    yyyymmdd = sys.argv[3]
    hh = sys.argv[4]
    mem_start = int(sys.argv[5])
    mem_end = int(sys.argv[6])
    pert_filename = sys.argv[7]
    scale_csv = sys.argv[8]
    max_cap = float(sys.argv[9])

    # Safely handle optional num_workers argument
    if len(sys.argv) > 10:
        num_workers = int(sys.argv[10])
    else:
        num_workers = mem_end - mem_start + 1

    n_mems = mem_end - mem_start + 1

    print("=" * 75)
    print("STEP 3: PARALLEL VERTICAL PROFILE RESCALING & OUTPUT DELIVERY")
    print("=" * 75)
    print(f" Work Directory   : {work_dir}")
    print(f" Target Output    : {out_root}/sfs.{yyyymmdd}/{hh}/")
    print(f" Target Members   : mem{mem_start:03d} to mem{mem_end:03d} ({n_mems} members)")
    print(f" Scale Profile    : {scale_csv}")
    print(f" Max Ratio Cap    : {max_cap:.2f}x")
    print(f" Parallel Workers : {num_workers}")
    print("-" * 75)

    if not os.path.exists(scale_csv):
        print(f"[FATAL ERROR] Scale CSV file not found: {scale_csv}")
        sys.exit(1)

    tasks = [
        (m, work_dir, out_root, yyyymmdd, hh, pert_filename, scale_csv, max_cap)
        for m in range(mem_start, mem_end + 1)
    ]

    print(f"Launching {n_mems} parallel tasks across {num_workers} worker(s)...\n")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = list(executor.map(process_member_rescale, tasks))

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
        print(f"[FATAL ERROR] Step 3 failed for {failures} member(s)!")
        sys.exit(1)
    else:
        print(f"✓ RESCALING COMPLETED SUCCESSFULLY!")
        print(f"  Delivered {n_mems} members ({total_mb:.1f} MB total) in {total_elapsed:.2f} seconds.")
        print("=" * 75)


if __name__ == "__main__":
    main()
