pro power_spectra_generation, ii

;Code in idl to generate power spectra of the EoR signal in an automated manner.
;Editing is still going on

;Editing code
;   <----------->: Indicates that the user should fill in the values applying to his code (!!!!!! should be made into a parameter file)
; !!!!!!!!!!!!!!!!: Something I have to or want to improve on.

;Computing parameters   <----------->
cpu, TPOOL_NTHREADS = 16           ;number of threads IDL may use


;Simulation parameters   <----------->
simulation_box_pixels = 300        ;Number of pixels of the box on a side
dimensions = 3                     ;Number of dimensions of the box
simulation_box_cmpc = 200          ;Size of a side of the box in comoving Mpc


;File locations    <----------->
loc_load_folder = '/net/machine/users/user/simulated_data/'
loc_save_folder = '/net/machine/users/user/processed_data/'

;Subfolders    <----------->
loc_load_d = loc_load_folder       ;Location of density files
loc_load_v = loc_load_folder       ;Location of velocity files
loc_load_x = loc_load_folder       ;Location of neutral fraction files
loc_load_T = loc_load_folder       ;Location of spin temperature files
loc_load_dTb = loc_load_folder     ;Location of 21cm field temperature files

loc_save_d = loc_save_folder       ;Location of density files
loc_save_v = loc_save_folder       ;Location of velocity files
loc_save_x = loc_save_folder       ;Location of neutral fraction files
loc_save_T = loc_save_folder       ;Location of spin temperature files
loc_save_dTb = loc_save_folder     ;Location of 21cm field temperature files


;In this example I only use a folder in which the density and neutral fields are available.
;The density field files look like : '/net/machine/users/user/simulated_data/6.905n_all.dat'
;The neutral field files look like : /net/machine/users/user/simulated_data/xfrac3d_6.905.bin'

;Structure of the file names    <----------->
filename_d_start = loc_load_d
filename_d_end = 'n_all.dat'

filename_x_start = loc_load_x + 'xfrac3d_'
filename_x_end = '.bin'

strlengt_filename_d_start = strlen(filename_d_start)
strlengt_filename_x_start = strlen(filename_x_start)

;We assume that the EoR happens between z of 0 and 100.    <----------->
number_of_decimals_z = 3


;Find all redshifts present
;In my case the redshift slices in the density folder do not always match the slices in the neutral field folder
;The following is to find matching files
density_files=FILE_SEARCH(filename_d_start+'*'+filename_d_end)
nslices_d=size(density_files,/n_elements)

neutral_files=FILE_SEARCH(filename_x_start+'*'+filename_x_end)
nslices_n=size(neutral_files,/n_elements)

;Create arrays for redshift information, both as float and string. 
zz=fltarr(nslices)
zzstr=strarr(nslices)
count=0

;Loop searches for all z slices for which both density and neutral field files are present
;The for loop can be removed using set operations. !!!!!!!!!!!!!!!!!
for i=0,nslices-1 do begin
     z_d = float(strmid(density_files(i),strlengt_filename_d_start,2+number_of_decimals_z))
     count2=0

     for j=0,nslices2-1 do begin
         if abs(AA - float(strmid(neutral_files(j),strlengt_filename_x_start,2+number_of_decimals_z))) le 0.01 then begin
         count2=1
         break
     endfor
    
     if count2 eq 1 then count = count +1
     if AA le 9.999 and count2 eq 1 then begin
          zzstr(count) = strmid(density_files(i),strlengt_filename_d_start,2+number_of_decimals_z)
          zz(count) = float(strmid(density_files(i),strlengt_filename_d_start,2+number_of_decimals_z))
     endif
     if AA ge 10. and count2 eq 1 thne begin
          zzstr(count) = strmid(density_files(i),strlengt_filename_d_start,3+number_of_decimals_z)
          zz(count) = float(strmid(density_files(i),strlengt_filename_d_start,3+number_of_decimals_z))
     endif
endfor

;Throw away empty bins in redshift arrays
zz=zz(1:count)
zzstr=zzstr(1:count)

;Create arrays for neutral fraction
neutral_volume_fraction = fltarr(count)
neutral_mass_fraction =  fltarr(count)

;Put files in order of increasing z
order=sort(zz)
zz=zz(order)
zzstr=zzstr(order)
delvar, order

for i=0,count do begin
    density_file = FILE_SEARCH(filename_d_start+zzstr(i)+filename_d_end)
    neutral_file = FILE_SEARCH(filename_x_start+zzstr(i)+filename_x_end)
    
    density_array = fltarr(simulation_box_pixels,simulation_box_pixels,simulation_box_pixels)
    neutral_array = fltarr(simulation_box_pixels,simulation_box_pixels,simulation_box_pixels)
    
    density_array = readdens, density_file, simulation_box_pixels
    neutral_array = readneutral, neutral_file, simulation_box_pixels
    
    ;compute neutral fraction (volume and mass weighted)
    neutral_volume_fraction[i],neutral_mass_fractio[i] = mean_neutral_fraction, neutral_field_array, density_field_array
    
    ;!!!!!!!!!!!!!!!! Include different reionization models, which I called GIO, LIO, LOI (see also Watkinson & Pritchard 201?5?) 
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ; Single=0: all, Single=1 Pd, =2 Px, =3 Pxd, =4 Cxd, =5 Cxdd, =6 Cxxd, =7 PdTb
    single=7
;     type='_'
    runps, density_array, neutral_array,1, zzstr(i),zz(i), loc_real,type+'fullreionization', single
;     runps, z_density_array, z_neutral_array,double(scaling_due_to_velocity),zzstr(i),zz(i), loc_reds, '_RSD_'+type, single
;     stop
;     if i eq nslices-1 then stop
endfor


end
     
