#!/usr/bin/env python3
import os
import sys
import time
import numpy as np
import xarray as xr


def compute_ensemble_mean():
    total_start_time = time.perf_counter()

    # 1. Parse Command Line Arguments
    if len(sys.argv) < 6:
        print(
            "[FATAL ERROR] Usage: python compute_ensemble_mean.py "
            "<origin_dir> <mem_start> <mem_end> <pert_filename> <out_mean_nc>"
        )
        sys.exit(1)

    origin_dir = sys.argv[1]
    mem_start = int(sys.argv[2])
    mem_end = int(sys.argv[3])
    pert_filename = sys.argv[4]
    out_mean_nc = sys.argv[5]

    n_mems = mem_end - mem_start + 1

    print("=" * 75)
    print("STEP 1: COMPUTING CPC ENSEMBLE MEAN (t_pert & s_pert)")
    print("=" * 75)
    print(f" Source Directory : {origin_dir}")
    print(f" Target Members   : mem{mem_start:03d} to mem{mem_end:03d} ({n_mems} total)")
    print(f" Input Filename   : {pert_filename}")
    print(f" Output Mean File : {out_mean_nc}")
    print("-" * 75)

    # 2. Build and Validate File List
    file_list = []
    for m in range(mem_start, mem_end + 1):
        fpath = os.path.join(
            origin_dir, f"mem{m:03d}", "analysis", "ocean", pert_filename
        )
        if not os.path.exists(fpath):
            print(f"[FATAL ERROR] Member perturbation file missing: {fpath}")
            sys.exit(1)
        file_list.append(fpath)

    print(f"✓ Validated presence of all {n_mems} ensemble member files.")

    # 3. Read Template Grid & Metadata from First Member
    t0_read = time.perf_counter()
    first_file = file_list[0]
    print(f"\n[1/3] Extracting metadata template from: mem{mem_start:03d}")

    ds_template = xr.open_dataset(first_file)

    # Validate required variables
    for var in ["t_pert", "s_pert"]:
        if var not in ds_template.data_vars:
            print(
                f"[FATAL ERROR] Variable '{var}' missing from template file!"
            )
            ds_template.close()
            sys.exit(1)

    t_dims = ds_template["t_pert"].dims
    s_dims = ds_template["s_pert"].dims
    t_attrs = ds_template["t_pert"].attrs
    s_attrs = ds_template["s_attrs"] if "s_attrs" in ds_template else ds_template["s_pert"].attrs

    # 4. Accumulate Sum Across All Members
    print(f"\n[2/3] Accumulation loop across {n_mems} member(s)...")

    t_sum = None
    s_sum = None

    for idx, fpath in enumerate(file_list, start=1):
        t0_mem = time.perf_counter()
        mem_name = f"mem{mem_start + idx - 1:03d}"

        with xr.open_dataset(fpath) as ds:
            # Extract raw numpy arrays for fast accumulation
            t_data = ds["t_pert"].values
            s_data = ds["s_pert"].values

            if t_sum is None:
                t_sum = np.zeros_like(t_data, dtype=np.float64)
                s_sum = np.zeros_like(s_data, dtype=np.float64)

            # Ignore NaN values during summation if present
            np.add(t_sum, np.nan_to_num(t_data, nan=0.0), out=t_sum)
            np.add(s_sum, np.nan_to_num(s_data, nan=0.0), out=s_sum)

        mem_elapsed = time.perf_counter() - t0_mem
        print(f"    [{idx:02d}/{n_mems:02d}] Accumulated {mem_name} ({mem_elapsed:.2f}s)")

    # 5. Calculate Average and Construct Dataset
    print("\n[3/3] Calculating mean and exporting NetCDF output...")
    t0_write = time.perf_counter()

    t_mean = (t_sum / n_mems).astype(np.float32)
    s_mean = (s_sum / n_mems).astype(np.float32)

    # Build clean output dataset maintaining coordinates
    ds_mean = xr.Dataset(
        data_vars={
            "t_pert": (t_dims, t_mean, t_attrs),
            "s_pert": (s_dims, s_mean, s_attrs),
        },
        coords={
            k: ds_template.coords[k]
            for k in ds_template.coords
            if k in t_dims or k in s_dims
        },
        attrs=ds_template.attrs,
    )

    # NetCDF zlib compression encoding
    encoding = {
        "t_pert": {"zlib": True, "complevel": 1},
        "s_pert": {"zlib": True, "complevel": 1},
    }

    # Ensure output directory exists
    os.makedirs(os.path.dirname(out_mean_nc), exist_ok=True)

    ds_mean.to_netcdf(out_mean_nc, encoding=encoding)
    ds_template.close()
    ds_mean.close()

    write_elapsed = time.perf_counter() - t0_write
    out_size_mb = os.path.getsize(out_mean_nc) / (1024 * 1024)

    print(f"    Exported: {out_mean_nc} ({out_size_mb:.2f} MB in {write_elapsed:.2f}s)")

    total_elapsed = time.perf_counter() - total_start_time
    print("=" * 75)
    print(f"✓ ENSEMBLE MEAN COMPUTED SUCCESSFULLY in {total_elapsed:.2f} seconds!")
    print("=" * 75)


if __name__ == "__main__":
    compute_ensemble_mean()
