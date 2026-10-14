#!/bin/bash

#generate the masks
#python


file_full="science/GL388/2539802o_pp_e2dsff_tcorr_AB.fits"
recon_file="science/GL388/2539802o_pp_e2dsff_recon_AB.fits"
wave_night_file="calib/6420F52BD3a_pp_e2dsff_AB_wave_night_AB.fits"
sky_mask_file="masks/3334061o_pp_e2dsff_AB_sky_mask_3sigma.fits"
file_full="science/GL388/2539802o_pp_e2dsff_tcorr_AB.fits"

#mask using a 3 sigma sky mask
dest_file="masked_3sigma_sky/GL388/2539802o_pp_e2dsff_tcorr_AB.fits"
python cp_mask_sky.py $file_full $wave_night_file $sky_mask_file $dest_file 

#mask using a 3 sigma sky mask and 50% absorption
dest_file="masked_trans50_3sigma_sky/GL388/2539802o_pp_e2dsff_tcorr_AB.fits"
python cp_mask_recon_sky.py  $file_full $wave_night_file $sky_mask_file $recon_file 0.5 $dest_file


