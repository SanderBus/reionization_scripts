PRO power_spectra_generation

;Code in idl to generate power spectra of the EoR signal in an automated manner.
;Editing is still going on

;Define all necessary variables by going through the parameter file
;This parameter file is located at
parameter_file_location = '/net/machine/data/users/user/parameterfile.txt'


P=CREATE_STRUCT(
        'nthreads',16,
        'simulation_box_pixels, 300,
        'dimensions', 3,
        'simulation_box_cmpc', 200,
        'simulation_name', 'EoRSIM_pix'+STR(simulation_box_pixels)+'_cmpc'+STR(simulation_box_cmpc),
        'ps_bins_num', 100,
        'ps_bins_type','linear',
        'which_spectra_to_compute', 0,
        'loc_load_folder', '/net/machine/data/users/user/simulated_data/,'
        'loc_save_folder', '/net/machine/data/users/user/processed_data/',
        'loc_save_folder_PS', loc_save_folder + 'PS/',
        'loc_save_folder_neutral_fraction', loc_save_folder +' neutral_fraction/',
        'loc_load_d', loc_load_folder,
        'loc_load_v', loc_load_folder,
        'loc_load_x', loc_load_folder,
        'loc_load_T', loc_load_folder,
        'loc_load_dTb', loc_load_folder,
        'filename_d_start', loc_load_d,
        'filename_d_end','n_all.dat',
        'filename_x_start', loc_load_x + 'xfrac3d_',
        'filename_x_end','.bin',
        'number_of_decimals_z',3)

;Update the parameter structure from the parameter file
read_parameter_file, parameter_file_location, parameters


;Find all redshifts present
;In my case the redshift slices in the density folder do not always match the slices in the neutral field folder
;The following is to find matching files
density_files=FILE_SEARCH(P.filename_d_start+'*'+P.filename_d_end)
nslices_d=SIZE(density_files,/n_elements)

neutral_files=FILE_SEARCH(P.filename_x_start+'*'+P.filename_x_end)
nslices_n=SIZE(neutral_files,/n_elements)

;Create arrays for redshift information, both as float and string. 
zz=FLTARR(nslices_n)
zzstr=FLTARR(nslices_n)
count=0

;Loop searches for all z slices for which both density and neutral field files are present
;In order to cut the z information out of the file name we need to know at what position the redshift number starts
strlengt_filename_d_start = STRLEN(P.filename_d_start)
strlengt_filename_x_start = STRLEN(P.filename_x_start)

;Start loop
FOR i=0,nslices-1 DO BEGIN
    z_d = FLOAT(STRMID(density_files(i),P.strlengt_filename_d_start,2+P.number_of_decimals_z))
    count2=0

    FOR j=0,nslices2-1 DO BEGIN
        IF ABS(z_d - FLOAT(STRMID(neutral_files(j),P.strlengt_filename_x_start,2+P.number_of_decimals_z))) LE 0.01 THEN BEGIN
        count2=1
        BREAK
    ENDFOR
    
    IF count2 EQ 1 THEN count = count +1
    IF AA LE 9.999 AND count2 EQ 1 THEN BEGIN
        zzstr(count) = STRMID(density_files(i),P.strlengt_filename_d_start,2+P.number_of_decimals_z)
        zz(count) = FLOAT(STRMID(density_files(i),P.strlengt_filename_d_start,2+P.number_of_decimals_z))
    ENDIF
    IF AA GE 10. and count2 EQ 1 THEN BEGIN
        zzstr(count) = STRMID(density_files(i),P.strlengt_filename_d_start,3+P.number_of_decimals_z)
        zz(count) = FLOAT(STRMID(density_files(i),P.strlengt_filename_d_start,3+P.number_of_decimals_z))
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
    density_file = FILE_SEARCH(P.filename_d_start+zzstr(i)+P.filename_d_end)
    neutral_file = FILE_SEARCH(P.filename_x_start+zzstr(i)+P.filename_x_end)
    
    density_array = FLTARR(P.simulation_box_pixels,P.simulation_box_pixels,P.simulation_box_pixels)
    neutral_array = FLTARR(P.simulation_box_pixels,P.simulation_box_pixels,P.simulation_box_pixels)
    
    density_array = readdens, density_file, P.simulation_box_pixels
    neutral_array = readneutral, neutral_file, P.simulation_box_pixels
    
    ;compute neutral fraction (volume and mass weighted)
    neutral_volume_fraction[i],neutral_mass_fractio[i] = mean_neutral_fraction, neutral_field_array, density_field_array
    
    dx = P.simulation_box_mpc / P.simulation_box_pixels  ;Size of a pixel in cMpc
    k=FLTARR(P.ps_bins_num)
     
    ;Make selected power spectra of all files
    IF single EQ 0 OR single EQ 7 OR single EQ 8 THEN BEGIN
        PdTb=FLTARR(P.ps_bins_num)
        save_loc_PdTb = P.loc_save_folder_PS + 'PdTb_bins'+P.ps_bins_num+P.ps_bins_type+'_z'+zzstr[i]+'_'+P.simulation_name+'.dat'
        dTb = 24*((1+zz)/10.)^0.5*((1+density_array)*neutral_array
        powerspec, dTb, dTb, dx, k, PdTb
        DELVAR, dTb
        save_array_as_txt, PdTb, P.save_loc_PdTb
        DELVAR, PdTb, save_loc_PdTb
    ENDIF
    IF single EQ 0 OR single EQ 1 OR single EQ 8 THEN BEGIN
        Pd=FLTARR(P.ps_bins_num)
        save_loc_Pd = P.loc_save_folder_PS + 'Pd_bins'+P.ps_bins_num+P.ps_bins_type+'_z'+zzstr[i]+'_'+P.simulation_name+'.dat'
        powerspec, density_array, density_array, dx, k, Pd
        save_array_as_txt, Pd, P.save_loc_Pd
        DELVAR, Pd, save_loc_Pd
    ENDIF     
    IF single EQ 0 OR single EQ 2 OR single EQ 8 THEN BEGIN
        Px=FLTARR(P.ps_bins_num)
        save_loc_Px= P.loc_save_folder_PS + 'Px_bins'+P.ps_bins_num+P.ps_bins_type+'_z'+zzstr[i]+'_'+P.simulation_name+'.dat'
        powerspec, neutral_array, neutral_array, dx, k, Px
        save_array_as_txt, Px, P.save_loc_Px
        DELVAR, Px, save_loc_Px
    ENDIF
    IF single EQ 0 OR single EQ 3 THEN BEGIN
        Pxd=FLTARR(P.ps_bins_num)
        save_loc_Pxd= P.loc_save_folder_PS + 'Pxd_bins'+P.ps_bins_num+P.ps_bins_type+'_z'+zzstr[i]+'_'+P.simulation_name+'.dat'
        powerspec, neutral_array*(1+density_array),neutral_array*(1+density_array), dx, k, Pxd
        save_array_as_txt, Pxd, P.save_loc_Pxd
        DELVAR, Pxd, save_loc_Pxd
    ENDIF
    IF single EQ 0 OR single EQ 4 THEN BEGIN
        Cxd = FLTARR(P.ps_bins_num)
        save_loc_Cxd = P.loc_save_folder_PS + 'Cxd_bins'+P.ps_bins_num+P.ps_bins_type+'_z'+zzstr[i]+'_'+P.simulation_name+'.dat'
        powerspec, neutral_array, density_array, dx, k, Cxd
        save_array_as_txt, Cxd, P.save_loc_Cxd
        DELVAR, Cxd, save_loc_Cxd
    ENDIF     
    IF single EQ 0 OR single EQ 5 THEN BEGIN
        Cxdd = FLTARR(P.ps_bins_num)
        save_loc_Cxdd = P.loc_save_folder_PS + 'Cxdd_bins'+P.ps_bins_num+P.ps_bins_type+'_z'+zzstr[i]+'_'+P.simulation_name+'.dat'
        powerspec, neutral_array*(1+density_array), density_array, dx, k, Cxdd
        save_array_as_txt, Cxdd, P.save_loc_Cxdd
        DELVAR, Cxdd, save_loc_Cxdd
    ENDIF          
    IF single EQ 0 OR single EQ 6 THEN BEGIN
        Cxxd = FLTARR(P.ps_bins_num)
        save_loc_Cxxd = P.loc_save_folder_PS + 'Cxxd_bins'+P.ps_bins_num+P.ps_bins_type+'_z'+zzstr[i]+'_'+P.simulation_name+'.dat'
        powerspec, neutral_array*(1+density_array), neutral_array, dx, k, Cxxd
        save_array_as_txt, Cxxd, P.save_loc_Cxxd
        DELVAR, Cxxd, save_loc_Cxxd
    ENDIF          
ENDFOR


;Save resp. wavenumbers used in PS calculation, neutral fraction by volume and mass and the redshifts all slices.
save_array_as_txt, k, P.loc_save_folder_PS+'k_bins'+P.ps_bins_num+P.ps_bins_type+'_sim'+P.simulation_name+'.dat'
save_array_as_txt, neutral_volume_fraction, P.loc_save_folder_neutral_fraction+'neutral_volume_fractions_list.dat'
save_array_as_txt, neutral_mass_fraction, P.loc_save_folder_neutral_fraction+'neutral_mass_fractions_list.dat'
save_array_as_txt, zz, P.loc_save_folder_neutral_fraction+'redshifts_list.dat'


END
