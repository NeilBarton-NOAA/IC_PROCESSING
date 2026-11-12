#!/usr/bin/env python3
import os
import sys
import numpy as np
import xarray as xr


def compare_recentered_vs_rescaled(recenter_nc, rescale_nc, target_lat=0.0, target_lon=180.0):
    print("=" * 115)
    print("STEP 3 SCALING CHECK: RE-CENTERED INTERMEDIATE vs. FINAL RESCALED")
    print("=" * 115)
    print(f" Re-Centered NC : {recenter_nc}")
    print(f" Final Rescaled : {rescale_nc}")
    print(f" Target Position : Lat {target_lat:.2f}°, Lon {target_lon:.2f}°")
    print("-" * 115)

    if not os.path.exists(recenter_nc):
        print(f"[FATAL ERROR] Re-centered file missing: {recenter_nc}")
        sys.exit(1)
    if not os.path.exists(rescale_nc):
        print(f"[FATAL ERROR] Rescaled file missing: {rescale_nc}")
        sys.exit(1)

    with xr.open_dataset(recenter_nc) as ds_rec, xr.open_dataset(rescale_nc) as ds_res:
        lats = ds_rec["lat"].values
        lons = ds_rec["lon"].values

        target_lon_adj = target_lon + 360.0 if (np.nanmax(lons) > 180 and target_lon < 0) else target_lon

        dist = (lats - target_lat) ** 2 + (lons - target_lon_adj) ** 2
        j_idx, i_idx = np.unravel_index(np.nanargmin(dist), dist.shape)

        actual_lat = lats[j_idx, i_idx]
        actual_lon = lons[j_idx, i_idx]

        print(f" Grid Cell Index: (j={j_idx}, i={i_idx}) -> Actual Lat {actual_lat:.2f}°, Lon {actual_lon:.2f}°")
        print("-" * 115)

        depths = ds_rec["Layer"].values
        rec_t = ds_rec["t_pert"].values[:, j_idx, i_idx]
        res_t = ds_res["t_pert"].values[:, j_idx, i_idx]
        rec_s = ds_rec["s_pert"].values[:, j_idx, i_idx]
        res_s = ds_res["s_pert"].values[:, j_idx, i_idx]

        header = f"{'Layer [1-75]':<13} | {'Depth (m)':<10} | {'Rec_t_pert':<13} | {'Res_t_pert':<13} | {'t_Scale Factor':<15} | {'s_Scale Factor':<15}"
        print(header)
        print("-" * 115)

        for k in range(len(depths)):
            layer_num = k + 1
            z_m = depths[k]
            
            rt, st = rec_t[k], res_t[k]
            rs_val, ss_val = rec_s[k], res_s[k]

            t_factor = (st / rt) if (abs(rt) > 1e-7 and not np.isnan(rt)) else 0.0
            s_factor = (ss_val / rs_val) if (abs(rs_val) > 1e-7 and not np.isnan(rs_val)) else 0.0

            rt_str = f"{rt:13.6f}" if not np.isnan(rt) else f"{'NaN':^13}"
            st_str = f"{st:13.6f}" if not np.isnan(st) else f"{'NaN':^13}"
            tf_str = f"{t_factor:15.3f}" if t_factor != 0.0 else f"{'N/A (0)':^15}"
            sf_str = f"{s_factor:15.3f}" if s_factor != 0.0 else f"{'N/A (0)':^15}"

            print(f"Layer [{layer_num:02d}]      | {z_m:8.2f}m | {rt_str} | {st_str} | {tf_str} | {sf_str}")

        print("=" * 115)


if __name__ == "__main__":
    # Updated paths matching your exact directory layout
    recenter_nc = "/scratch3/NCEPDEV/climate/Jiande.Wang/SFS-work/re-center-git/work-recenter-ed/20260701/00/mem005/analysis/ocean/20260701.000000.mom6_perturbation.nc"
    rescale_nc = "/scratch3/NCEPDEV/climate/Jiande.Wang/SFS-work/re-center-git/ICs/CPC_land-rescaled/C192mx025/sfs.20260701/00/mem005/analysis/ocean/20260701.000000.mom6_perturbation.nc"

    lat_val = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
    lon_val = float(sys.argv[2]) if len(sys.argv) > 2 else 180.0

    if len(sys.argv) > 3:
        recenter_nc = sys.argv[3]
    if len(sys.argv) > 4:
        rescale_nc = sys.argv[4]

    compare_recentered_vs_rescaled(recenter_nc, rescale_nc, lat_val, lon_val)
