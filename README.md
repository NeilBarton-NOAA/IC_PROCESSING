Scripts to Process GFS ICs, GFS_DA ICs CPC_land ICs, CPC_landice ICs, and Replay ICs for NOAA Seasonal Forecast System (SFS)
Latest updates were made on October 2026

To get the IC_PROCESSING package:
(1) Git clone from https://github.com/YangxingZheng-NOAA/IC_PROCESSING.git or https://github.com/NeilBarton-NOAA/IC_PROCESSING.git
(2) Run chgres_compile.sh to get UFS_UTILS

#############################################################################################
!!! SPECIAL NOTES TO GET the GLORE ICS on GaeaC5:

Option 1: Here we use Globus to get GLORE ICs, users need to login in to the Globus CLI in the terminal and perform a manual tranfer before using the package. Below is an example on Ursa to transfer from GaeaC5
(A) Login to Globus CLI
    1. Load the module: module load globus-cli
    2. Login in to Globus: globus login
(B) Perform a Manual Tranfer
    Perform a manual transfer using the endpoints in your script immediately after logging in. This will activate your specific source and destination endpoints (UUID_C5 and UUID_DES) and cache the required credentials in ~/.globus/:

    # Get the folder paths by stripping the file names
    UUID_C5="57da978b-4434-4ac0-ae5e-3de5161305e1"
    UUID_DES="8bbc822a-a625-4058-b68e-65edacd73828" #(this is for Ursa)
    dir_in=$(dirname "${file_in}")
    dir_out=$(dirname "${file_out}")
    globus transfer "${UUID_C5}:${dir_in}/" "${UUID_DES}:${dir_out}/" --recursive --dry-run
    where you need to specify full path of file_in (i.e., dir_in) and file_out (i.e., dir_out) before running the above globus transfer
    
    When you run your globus transfer command with --dry-run, Globus will go through all the preparatory steps of your transfer without actually moving a single file

Option 2: You can manually transfer the GLORE data using scp or rsync or on-line Globus freely.
#################################################################################################

To Generate the ICs:

FOR GFS ICS: run in the following order
    (1) GET_ICS.sh ${DTG} GFS
    (2) RUN_CHGRES.sh ${DTG} GFS
    (3) RUN_REGRID.sh ${DTG} GFS

FOR CPC ICS update with LAND SPIN-UP (CPC_land): run in the following order
    (1) GET_ICS.sh ${DTG} CPC_land
    (2) RUN_CHNGRES.sh ${DTG} CPC_land
    (3) LAND_REPLACE.sh ${DTG} CPC_land
    (4) LINK.sh ${DTG} CPC_land

FOR CPC ICS update with LAND SPIN-UP and SNOW CLIMATOLOGY in ice ICs (CPC_landice): run in the following order
    (1) GET_ICS.sh ${DTG} CPC_landice
    (2) RUN_CHNGRES.sh ${DTG} CPC_landice
    (3) LAND_ICE_REPLACE.sh ${DTG} CPC_landice
    (4) LINK.sh ${DTG} CPC_landice

FOR GFS ICS WITH UPDATED MARINE ICs (GFS_DA): run in the following order
    (1) GET_ICS.sh ${DTG} GFS_DA
    (2) RUN_CHGRES.sh ${DTG} GFS_DA
    (3) LINK.sh ${DTG} GFS_DA

FOR REPLAY ICS: run in the following order
    (1) GET_ICS.sh ${DTG} REPLAY
    (2) RUN_CHNGRES.sh ${DTG} REPLAY
    (3) LINK.sh ${DTG} REPLAY
