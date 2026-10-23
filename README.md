Scripts to Process GFS ICs, GFS_DA ICs CPC_land ICs, CPC_landice ICs, and Replay ICs for NOAA Seasonal Forecast System (SFS)
Latest updates were made on October 2026

To get the IC_PROCESSING package:
(1) Git clone from https://github.com/YangxingZheng-NOAA/IC_PROCESSING.git or https://github.com/NeilBarton-NOAA/IC_PROCESSING.git
(2) Run chgres_compile.sh to get UFS_UTILS


To Generate the ICs:

FOR GFS ICS: run in the following order
    (1) GET_ICS.sh ${DTG} GFS
    (2) RUN_CHGRES.sh ${DTG} GFS
    (3) RUN_REGRID.sh ${DTG} GFS

FOR CPC ICS update with LAND SPIN-UP (CPC_land): run in the following order
    (1) GET_ICS.sh ${DTG} CPC_land
    (2) GET_PERTIRBATIONS.sh ${DTG} CPC_land
    (3) RUN_CHNGRES.sh ${DTG} CPC_land
    (4) LAND_REPLACE.sh ${DTG} CPC_land
    (5) LINK.sh ${DTG} CPC_land

FOR CPC ICS update with LAND SPIN-UP and SNOW CLIMATOLOGY in ice ICs (CPC_landice): run in the following order
    (1) GET_ICS.sh ${DTG} CPC_landice
    (2) GET_PERTIRBATIONS.sh ${DTG} CPC_landice
    (3) RUN_CHNGRES.sh ${DTG} CPC_landice
    (4) LAND_REPLACE.sh ${DTG} CPC_landice
    (5) ICE_REPLACE.sh ${DTG} CPC_landice
    (5) LINK.sh ${DTG} CPC_landice

FOR GFS ICS WITH UPDATED MARINE ICs (GFS_DA): run in the following order
    (1) GET_ICS.sh ${DTG} GFS_DA
    (2) RUN_CHGRES.sh ${DTG} GFS_DA
    (3) GET_PERTIRBATIONS.sh ${DTG} GFS_DA
    (4) LINK.sh ${DTG} GFS_DA
    (5) FINALIZE_GFS_DA_DIR.sh ${DTG} GFS_DA

FOR REPLAY ICS: run in the following order
    (1) GET_ICS.sh ${DTG} REPLAY
    (1) GET_PERTIRBATIONS.sh ${DTG} REPLAY
    (2) RUN_CHNGRES.sh ${DTG} REPLAY
    (3) LINK.sh ${DTG} REPLAY
