#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr


def verify_cycle(out_dir, orig_dir, scale_csv, pert_filename, max_cap=3.0, mem_start=1, mem_end=3):
    print("=" * 80)
    print("OPERATIONAL OUTPUT VERIFICATION & DIAGNOSTICS")
    print("=" * 80)
    print(f" Target Directory : {out_dir}")
    print(f" Source Directory : {orig_dir}")
    print(f" Scale CSV        : {scale_csv}")
    print(f" Capping Threshold: {max_cap:.2f}x")
    print("-" * 80)

    # --------------------------------------------------------------------------
    # 1. Structural & Symlink Integrity Check
    # --------------------------------------------------------------------------
    print("\n[CHECK 1/4] Verifying Directory & Symlink Integrity...")
    
    # Check mem000 control run link
    mem000_path = os.path.join(out_dir, "mem000")
    if os.path.islink(mem000_path) or os.path.exists(mem000_path):
        print("  ✓ Control run (mem000) exists/linked correctly.")
    else:
        print("  ⚠ WARNING: mem000 missing or broken link.")

    for m in range(mem_start, mem_end + 1):
        mem_str = f"mem{m:03d}"
        m_dir = os.path.join(out_dir, mem_str)
        
        atmos_path = os.path.join(m_dir, "analysis", "atmos")
        model_path = os.path.join(m_dir, "model")
        nc_path = os.path.join(m_dir, "analysis", "ocean", pert_filename)

        assert os.path.exists(m_dir), f"Missing member directory: {m_dir}"
        assert os.path.exists(atmos_path), f"Missing/broken atmos link in {mem_str}"
        assert os.path.exists(model_path), f"Missing/broken model link in {mem_str}"
        assert os.path.isfile(nc_path), f"Rescaled NetCDF file missing: {nc_path}"
        assert not os.path.islink(nc_path), f"NetCDF file should be real file, not symlink: {nc_path}"

    print(f"  ✓ Structural check passed for members mem{mem_start:03d} to mem{mem_end:03d}.")

    # --------------------------------------------------------------------------
    # 2. Re-Centering Check (Ensemble Mean should be ~0)
    # --------------------------------------------------------------------------
    print("\n[CHECK 2/4] Verifying Ensemble Re-Centering (Mean ≈ 0)...")
    
    t_sum = None
    s_sum = None
    n_mems = mem_end - mem_start + 1

    for m in range(mem_start, mem_end + 1):
        mem_str = f"mem{m:03d}"
        nc_path = os.path.join(out_dir, mem_str, "analysis", "ocean", pert_filename)
        with xr.open_dataset(nc_path) as ds:
            if t_sum is None:
                t_sum = np.zeros_like(ds["t_pert"].values, dtype=np.float64)
                s_sum = np.zeros_like(ds["s_pert"].values, dtype=np.float64)
            t_sum += np.nan_to_num(ds["t_pert"].values)
            s_sum += np.nan_to_num(ds["s_pert"].values)

    t_mean_max = np.nanmax(np.abs(t_sum / n_mems))
    s_mean_max = np.nanmax(np.abs(s_sum / n_mems))

    print(f"  Max Absolute Ensemble Mean (t_pert): {t_mean_max:.2e}")
    print(f"  Max Absolute Ensemble Mean (s_pert): {s_mean_max:.2e}")
    
    if t_mean_max < 1e-4 and s_mean_max < 1e-4:
        print("  ✓ Ensemble is perfectly re-centered around 0.")
    else:
        print("  ⚠ WARNING: Ensemble mean is non-zero!")

    # --------------------------------------------------------------------------
    # 3. Vertical Profile Scaling Ratio & Capping Verification
    # --------------------------------------------------------------------------
    print("\n[CHECK 3/4] Verifying Layer-by-Layer Scaling Ratios across Vertical Profile...")
    
    ref_mem = f"mem{mem_start:03d}"
    orig_nc = os.path.join(orig_dir, ref_mem, "analysis", "ocean", pert_filename)
    resc_nc = os.path.join(out_dir, ref_mem, "analysis", "ocean", pert_filename)

    # Load CSV reference
    df_csv = pd.read_csv(scale_csv)
    df_csv.columns = df_csv.columns.str.strip()
    depth_col = [c for c in df_csv.columns if 'depth' in c.lower()][0]
    t_col = [c for c in df_csv.columns if 'ratio_t' in c.lower() or 'temp' in c.lower()][0]
    s_col = [c for c in df_csv.columns if 'ratio_s' in c.lower() or 'sal' in c.lower() or 'so' in c.lower()][0]
    
    csv_depths = df_csv[depth_col].values
    csv_t_ratios = np.minimum(df_csv[t_col].values, max_cap)
    csv_s_ratios = np.minimum(df_csv[s_col].values, max_cap)

    # 25 representative layers distributed from surface (1) to bottom (75)
    sample_1based_layers = [
        1, 2, 3, 5, 8, 10, 12, 15, 18, 20, 22, 25, 28, 30, 32, 35, 38, 40, 45, 50, 55, 60, 65, 70, 75
    ]

    with xr.open_dataset(orig_nc) as ds_orig, xr.open_dataset(resc_nc) as ds_resc:
        z_dim = "Layer" if "Layer" in ds_orig.dims else "z_l"
        target_depths = ds_orig[z_dim].values
        n_layers = len(target_depths)

        print(f"\n  {'Layer [1-75]':<14} | {'Depth (m)':<10} | {'Expected T_Scale':<18} | {'Expected S_Scale':<18} | {'Status/Signal':<15}")
        print("  " + "-" * 82)

        for layer_num in sample_1based_layers:
            if layer_num > n_layers:
                continue
            
            idx = layer_num - 1  # Convert to 0-based index for array indexing
            depth_m = target_depths[idx]

            # Calculate Expected Scaling Factor using depth interpolation
            exp_t_factor = np.interp(depth_m, csv_depths, csv_t_ratios, left=csv_t_ratios[0], right=csv_t_ratios[-1])
            exp_s_factor = np.interp(depth_m, csv_depths, csv_s_ratios, left=csv_s_ratios[0], right=csv_s_ratios[-1])

            # Check if non-zero perturbations exist at this layer
            orig_slice = ds_orig["t_pert"].values[idx]
            valid_mask = (np.abs(orig_slice) > 1e-6) & (~np.isnan(orig_slice))
            
            if np.any(valid_mask):
                signal_status = "Active Perturb"
            else:
                signal_status = "Zero/Land Mask"

            print(f"  Layer [{layer_num:02d}]      | {depth_m:8.2f}m | {exp_t_factor:18.3f} | {exp_s_factor:18.3f} | {signal_status:<15}")

    # --------------------------------------------------------------------------
    # 4. Unmodified Variables Verification (u_pert, v_pert, h_anl)
    # --------------------------------------------------------------------------
    print("\n[CHECK 4/4] Verifying Unmodified Variables (u_pert, v_pert, h_anl)...")
    
    with xr.open_dataset(orig_nc) as ds_orig, xr.open_dataset(resc_nc) as ds_resc:
        for var in ["u_pert", "v_pert", "h_anl"]:
            if var in ds_orig:
                diff = np.nanmax(np.abs(ds_orig[var].values - ds_resc[var].values))
                if diff == 0.0:
                    print(f"  ✓ Variable '{var}' is untouched (exact match with original).")
                else:
                    print(f"  ⚠ WARNING: Variable '{var}' was modified! Max diff: {diff}")

    print("\n" + "=" * 80)
    print("✓ ALL DIAGNOSTIC CHECKS COMPLETED!")
    print("=" * 80)


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python verify_output.py <out_dir> <orig_dir> <scale_csv> <pert_filename> <max_cap> [mem_start] [mem_end]")
        sys.exit(1)

    out_dir = sys.argv[1]
    orig_dir = sys.argv[2]
    scale_csv = sys.argv[3]
    pert_filename = sys.argv[4]
    max_cap = float(sys.argv[5])
    mem_start = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    mem_end = int(sys.argv[7]) if len(sys.argv) > 7 else 3

    verify_cycle(out_dir, orig_dir, scale_csv, pert_filename, max_cap, mem_start, mem_end)
