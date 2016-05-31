PRO power_spectra_generation

;Code in idl to generate power spectra of the EoR signal in an automated manner.
;Editing is still going on

;Editing code
;   <----------->: Indicates that the user should fill in the values applying to his code (!!!!!! should be made into a parameter file)
; !!!!!!!!!!!!!!!!: Something I have to or want to improve on.

;Computing parameters   <----------->
CPU, TPOOL_NTHREADS = 16           ;number of threads IDL may use


;Simulation parameters   <----------->
simulation_box_pixels = 300        ;Number of pixels of the box on a side
dimensions = 3                     ;Number of dimensions of the box
simulation_box_cmpc = 200          ;Size of a side of the box in comoving Mpc
simulation_name = 'EoRSIM_pix'+STR(simulation_box_pixels)+'_cmpc'+STR(simulation_box_cmpc)

;PS extraction parameters  <----------->
ps_bins_num = 100
ps_bins_type = 'linear'      ;linear or logarithmic
; Decide which spectra to compute
; which_spectra_to_compute=0: all, which_spectra_to_compute=1 Pd, =2 Px, =3 Pxd, =4 Cxd, =5 Cxdd, =6 Cxxd, =7 PdTb, =8 Pd, Px and PdTb
which_spectra_to_compute = 0


;File locations    <----------->
loc_load_folder = '/net/machine/data/users/user/simulated_data/'
loc_save_folder = '/net/machine/data/users/user/processed_data/'

loc_save_folder_PS = loc_save_folder + 'PS/'
loc_save_folder_neutral_fraction = loc_save_folder +' neutral_fraction/'

;Subfolders    <----------->
loc_load_d = loc_load_folder       ;Location of density files
loc_load_v = loc_load_folder       ;Location of velocity files
loc_load_x = loc_load_folder       ;Location of neutral fraction files
loc_load_T = loc_load_folder       ;Location of spin temperature files
loc_load_dTb = loc_load_folder     ;Location of 21cm field temperature files

;In this example I only use a folder in which the density and neutral fields are available.
;The density field files look like : '/net/machine/users/user/simulated_data/6.905n_all.dat'
;The neutral field files look like : /net/machine/users/user/simulated_data/xfrac3d_6.905.bin'

;Structure of the file names    <----------->
filename_d_start = loc_load_d
filename_d_end = 'n_all.dat'

filename_x_start = loc_load_x + 'xfrac3d_'
filename_x_end = '.bin'

strlengt_filename_d_start = STRLEN(filename_d_start)
strlengt_filename_x_start = STRLEN(filename_x_start)

;We assume that the EoR happens between z of 0 and 100.    <----------->
number_of_decimals_z = 3


;Find all redshifts present
;In my case the redshift slices in the density folder do not always match the slices in the neutral field folder
;The following is to find matching files
density_files=FILE_SEARCH(filename_d_start+'*'+filename_d_end)
nslices_d=SIZE(density_files,/n_elements)

neutral_files=FILE_SEARCH(filename_x_start+'*'+filename_x_end)
nslices_n=SIZE(neutral_files,/n_elements)

;Create arrays for redshift information, both as float and string. 
zz=FLTARR(nslices)
zzstr=FLTARR(nslices)
count=0

;Loop searches for all z slices for which both density and neutral field files are present
;The for loop can be removed using set operations. !!!!!!!!!!!!!!!!!
FOR i=0,nslices-1 DO BEGIN
    z_d = FLOAT(STRMID(density_files(i),strlengt_filename_d_start,2+number_of_decimals_z))
    count2=0

    FOR j=0,nslices2-1 DO BEGIN
        IF ABS(AA - FLOAT(STRMID(neutral_files(j),strlengt_filename_x_start,2+number_of_decimals_z))) LE 0.01 THEN BEGIN
        count2=1
        BREAK
    ENDFOR
    
    IF count2 EQ 1 THEN count = count +1
    IF AA LE 9.999 AND count2 EQ 1 THEN BEGIN
        zzstr(count) = STRMID(density_files(i),strlengt_filename_d_start,2+number_of_decimals_z)
        zz(count) = FLOAT(STRMID(density_files(i),strlengt_filename_d_start,2+number_of_decimals_z))
    ENDIF
    IF AA GE 10. and count2 EQ 1 THEN BEGIN
        zzstr(count) = STRMID(density_files(i),strlengt_filename_d_start,3+number_of_decimals_z)
        zz(count) = FLOAT(STRMID(density_files(i),strlengt_filename_d_start,3+number_of_decimals_z))
    ENDIF
ENDFOR

;Throw away empty bins in redshift arrays
zz=zz(1:count)
zzstr=zzstr(1:count)

;Create arrays for neutral fraction
neutral_volume_fraction = FLTARR(count)
neutral_mass_fraction =  FLTARR(count)

;Put files in order of increasing z
order=SORT(zz)
zz=zz(order)
zzstr=zzstr(order)
delvar, order

for i=0,count do begin
    density_file = FILE_SEARCH(filename_d_start+zzstr(i)+filename_d_end)
    neutral_file = FILE_SEARCH(filename_x_start+zzstr(i)+filename_x_end)
    
    density_array = FLTARR(simulation_box_pixels,simulation_box_pixels,simulation_box_pixels)
    neutral_array = FLTARR(simulation_box_pixels,simulation_box_pixels,simulation_box_pixels)
    
    density_array = readdens, density_file, simulation_box_pixels
    neutral_array = readneutral, neutral_file, simulation_box_pixels
    
    ;compute neutral fraction (volume and mass weighted)
    neutral_volume_fraction[i],neutral_mass_fractio[i] = mean_neutral_fraction, neutral_field_array, density_field_array
    
    dx = simulation_box_mpc / simulation_box_pixels  ;Size of a pixel in cMpc
    k=FLTARR(ps_bins_num)
     
    ;Make selected power spectra of all files
    IF single EQ 0 OR single EQ 7 OR single EQ 8 THEN BEGIN
        PdTb=FLTARR(ps_bins_num)
        save_loc_PdTb = loc_save_folder_PS + 'PdTb_z'+zzstr[i]+'_'+simulation_name+'.dat'
        dTb = 24*((1+zz)/10.)^0.5*((1+density_array)*neutral_array
        powerspec, dTb, dTb, dx, k, PdTb
        DELVAR, dTb
        save_array_as_txt, PdTb, save_loc_PdTb
        DELVAR, PdTb, save_loc_PdTb
    ENDIF
    IF single EQ 0 OR single EQ 1 OR single EQ 8 THEN BEGIN
        Pd=FLTARR(ps_bins_num)
        save_loc_Pd = loc_save_folder_PS + 'Pd_bins'+ps_bins_num+ps_bins_type+'_z'+zzstr[i]+'_'+simulation_name+'.dat'
        powerspec, density_array, density_array, dx, k, Pd
        save_array_as_txt, Pd, save_loc_Pd
        DELVAR, Pd, save_loc_Pd
    ENDIF     
    IF single EQ 0 OR single EQ 2 OR single EQ 8 THEN BEGIN
        Px=FLTARR(ps_bins_num)
        save_loc_Px= loc_save_folder_PS + 'Px_bins'+ps_bins_num+ps_bins_type+'_z'+zzstr[i]+'_'+simulation_name+'.dat'
        powerspec, neutral_array, neutral_array, dx, k, Px
        save_array_as_txt, Px, save_loc_Px
        DELVAR, Px, save_loc_Px
    ENDIF
    IF single EQ 0 OR single EQ 3 THEN BEGIN
        Pxd=FLTARR(ps_bins_num)
        save_loc_Pxd= loc_save_folder_PS + 'Pxd_bins'+ps_bins_num+ps_bins_type+'_z'+zzstr[i]+'_'+simulation_name+'.dat'
        powerspec, neutral_array*(1+density_array),neutral_array*(1+density_array), dx, k, Pxd
        save_array_as_txt, Pxd, save_loc_Pxd
        DELVAR, Pxd, save_loc_Pxd
    ENDIF
    IF single EQ 0 OR single EQ 4 THEN BEGIN
        Cxd = FLTARR(ps_bins_num)
        save_loc_Cxd = loc_save_folder_PS + 'Cxd_bins'+ps_bins_num+ps_bins_type+'_z'+zzstr[i]+'_'+simulation_name+'.dat'
        powerspec, neutral_array, density_array, dx, k, Cxd
        save_array_as_txt, Cxd, save_loc_Cxd
        DELVAR, Cxd, save_loc_Cxd
    ENDIF     
    IF single EQ 0 OR single EQ 5 THEN BEGIN
        Cxdd = FLTARR(ps_bins_num)
        save_loc_Cxdd = loc_save_folder_PS + 'Cxdd_bins'+ps_bins_num+ps_bins_type+'_z'+zzstr[i]+'_'+simulation_name+'.dat'
        powerspec, neutral_array*(1+density_array), density_array, dx, k, Cxdd
        save_array_as_txt, Cxdd, save_loc_Cxdd
        DELVAR, Cxdd, save_loc_Cxdd
    ENDIF          
    IF single EQ 0 OR single EQ 6 THEN BEGIN
        Cxxd = FLTARR(ps_bins_num)
        save_loc_Cxxd = loc_save_folder_PS + 'Cxxd_bins'+ps_bins_num+ps_bins_type+'_z'+zzstr[i]+'_'+simulation_name+'.dat'
        powerspec, neutral_array*(1+density_array), neutral_array, dx, k, Cxxd
        save_array_as_txt, Cxxd, save_loc_Cxxd
        DELVAR, Cxxd, save_loc_Cxxd
    ENDIF          
ENDFOR

;Save resp. wavenumbers used in PS calculation, neutral fraction by volume and mass and the redshifts all slices.
save_array_as_txt, k, loc_save_folder_PS+'k.dat'
save_array_as_txt, neutral_volume_fraction, loc_save_folder_neutral_fraction+'neutral_volume_fractions_list.dat'
save_array_as_txt, neutral_mass_fraction, loc_save_folder_neutral_fraction+'neutral_mass_fractions_list.dat'
save_array_as_txt, zz, loc_save_folder_neutral_fraction+'redshifts_list.dat'

END
