# Copy raw data to data_clean
cp data_raw/kittel_et_al_2021/MARcst-AN35km-176x148.cdf data_clean/mar_model/

# Generate mean-MAR_ACCESS1.3-1980-1989_ams.nc, which creates the grid_imask variable
python schneida_tools/mar3.py

# Write grid_imask variable to MAR data
ncks -A -C -v grid_imask data_clean/mar_model/mean-MAR_ACCESS1.3-1980-1989_ams.nc data_clean/mar_model/MARcst-AN35km-176x148.cdf

# infer grid file file from MAR data 
ncremap -d data_clean/mar_model/MARcst-AN35km-176x148.cdf -g data_clean/mar_model/MARcst-AN35km-176x148_grd.nc

# Create conservative mapping file
ncremap -a conserve -s data_clean/mar_model/MARcst-AN35km-176x148_grd.nc -g data_clean/mar_model/MARcst-AN35km-176x148_grd.nc -m data_clean/map_marv3.11_to_marv3.11_aave.20220818.nc

# Apply mapping to regrid MAR data
ncremap --sgs_frc=solid_mask --sgs_msk=grid_imask -m data_clean/map_marv3.11_to_marv3.11_aave.20220818.nc data_clean/mar_model/mean-MAR_ACCESS1.3-1980-1989_ams.nc data_clean/mar_model/mean-MAR_ACCESS1.3-1980-1989_amsrgr.nc