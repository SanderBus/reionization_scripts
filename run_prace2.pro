pro run_prace2, ii

cpu, TPOOL_NTHREADS = 16           ;number of threads IDL may use

; File locations
loc_load_folder = '/net/machine/users/user/simulated_data/'
loc_save_folder = '/net/machine/users/user/processed_data/'

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
;The neutral field files look like : /net/machine/users/user/simulated_data/xfrac3d_6.905.bin/'

;Find all redshifts present
;In my case the redshift slices in the density folder do not always match the slices in the neutral field folder
;The following is to find matching files
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
density_files=FILE_SEARCH(loc_load_d+'*n_all.dat')
nslices=size(density_files,/n_elements)

neutral_files=FILE_SEARCH(loc_x+'xfrac3d_*')
nslices2=size(neutral_files,/n_elements)
zz=fltarr(nslices)
zzstr=strarr(nslices)
ll=0

for i=0,nslices-1 do begin
     AA = float(strmid(density_files(i),99,5))
     kk=0

     for j=0,nslices2-1 do begin
         if abs(AA - float(strmid(neutral_files(j),217,5))) le 0.01 then kk=1
         print, strmid(neutral_files(j),217,5)
     endfor
    
     if kk eq 1 then ll = ll +1
     if AA le 9.999 and kk eq 1 then zzstr(ll) = strmid(density_files(i),99,5)
     if AA le 9.999 and kk eq 1 then zz(ll) = float(strmid(density_files(i),99,5))
     if AA ge 10. and kk eq 1 then zzstr(ll) = strmid(density_files(i),99,6)
     if AA ge 10. and kk eq 1 then zz(ll) = float(strmid(density_files(i),99,6))
endfor

zz=zz(1:ll)
zzstr=zzstr(1:ll)

;Put files in order of increasing z
order=sort(zz)
zz=zz(order)
density_files=density_files(order)
zzstr=zzstr(order)
delvar, order


for i=0,ll do begin
;     i+=11
    density_file = FILE_SEARCH(loc_vd+zzstr(i)+'n_all.dat')
    neutral_file = FILE_SEARCH(loc_x+'xfrac3d_'+zzstr(i)+'.bin')
    velocity_file = FILE_SEARCH(loc_vd+zzstr(i)+'v_all.dat')
    neutral_frac_file = '/net/osiris/data/users/sanbus/PRACE4LOFAR2/real_space/GIO/neutralfracs/neutralfrac__'+zzstr(i)+'.dat'
    
    density_array = fltarr(300,300,300)
    density_array_full = fltarr(600,600,600)
    z_density_array = fltarr(300,300,300)
    machine_density_array = fltarr(300,300,300)
    neutral_array = fltarr(300,300,300)
    z_neutral_array = fltarr(300,300,300)
    velocity_array = fltarr(300,300,300)
    scaling_due_to_velocity = fltarr(300,300,300)
    
    zaxis=fltarr(300)
    z3D=fltarr([300,300,300])
    sizze=0
    
    print, density_file
    readdens_prace2, density_file, density_array_full, machine_density_array,sizze
    
    if sizze eq 600 then begin
        density_array = rebin(density_array_full,300L,300L,300L)
    endif else begin 
        density_array = density_array_full
    endelse
    delvar, density_array_full
    


    if type eq 'GIO' then readneutral, neutral_file, neutral_array else read_neutral_frac, neutral_frac_file, neutralfraction
    if type eq 'GOI' then begin
    loc_conjugate_neutralfrac_list=0
    goi_prace, zzstr(i),zzstr, loc_conjugate_neutralfrac
    neutral_frac_file = '/net/osiris/data/users/sanbus/PRACE4LOFAR2/real_space/GIO/neutralfracs/neutralfrac__z'+zzstr(loc_conjugate_neutralfrac)+'.dat'
    readneutral, neutral_file, neutral_array
;     neutral_array= 1 -neutral_array
    endif
    if type eq 'LIO' then begin
    lio_loi_reionization, density_array, neutralfraction, 1, neutral_array
    endif
    if type eq 'LOI' then begin
    lio_loi_reionization, density_array, neutralfraction, 0, neutral_array
    endif
    
    print, min(density_array), mean(density_array), max(density_array)
    print, min(neutral_array), mean(neutral_array), max(neutral_array)
    
; ; ; ; ; ; ; ; ;     ; Make redshift array
; ; ; ; ; ; ; ; ;     makez3d, zz(i), zaxis, z3D
; ; ; ; ; ; ; ; ;     readvel, velocity_file, machine_density_array, zz(i), velocity_array
; ; ; ; ; ; ; ; ;     delvar, machine_density_array
; ; ; ; ; ; ; ; ;     print, 'density params', min(density_array),mean(density_array),max(density_array)
; ; ; ; ; ; ; ; ;     print, 'neutral params',  min(neutral_array),mean(neutral_array),max(neutral_array)
; ; ; ; ; ; ; ; ;     print, 'velocity params', min(abs(velocity_array)),mean(abs(velocity_array)),max(abs(velocity_array))
; ; ; ; ; ; ; ; ;     
; ; ; ; ; ; ; ; ;     
; ; ; ; ; ; ; ; ;     vel_scale,velocity_array/100000.,zz(i),500./300.,scaling_due_to_velocity
; ; ; ; ; ; ; ; ;     
; ; ; ; ; ; ; ; ;     rsd, 1+density_array,neutral_array, velocity_array, zaxis, z_density_array,z_neutral_array
; ; ; ; ; ; ; ; ;     z_density_array -=1
; ; ; ; ; ; ; ; ;     delvar, velocity_array
; ; ; ; ; ; ; ; ;     
; ; ; ; ; ; ; ; ;     print, 'z density params', min(z_density_array),mean(z_density_array),max(z_density_array)
; ; ; ; ; ; ; ; ;     print, 'z neutral params',  min(z_neutral_array),mean(z_neutral_array),max(z_neutral_array)    
; ; ; ; ; ; ; ; ;     print, 'scaling due to vel params',  min(scaling_due_to_velocity),mean(scaling_due_to_velocity),max(scaling_due_to_velocity)
; ; ; ; ; ; ; ; ;     
; ; ; ; ; ; ; ; ;     
; ; ; ; ; ; ; ; ;     ;Save slices of all fields
; ; ; ; ; ; ; ; ;     run_saveslice, density_array, z_density_array, neutral_array, z_neutral_array, scaling_due_to_velocity,zzstr(i),zz(i),loc_real,loc_reds
; ; ; ; ; ; ; ; ;     
; ; ; ; ; ; ; ; ;     ;Sace redshifted cubes
; ; ; ; ; ; ; ; ;     save3darray_bin, z_density_array, '/net/osiris/data/users/sanbus/PRACE4LOFAR/redshift_space/'+type+'/Cubes/density'+type+'_cube_rsd_'+zzstr(i)+'.bin'
; ; ; ; ; ; ; ; ;     save3darray_bin, z_neutral_array, '/net/osiris/data/users/sanbus/PRACE4LOFAR/redshift_space/'+type+'/Cubes/neutral'+type+'_cube_rsd_'+zzstr(i)+'.bin'
; ; ; ; ; ; ; ; ;     save3darray_bin, double(24*((1+zz(i))/10.)^0.5*z_neutral_array*(1+z_density_array)*scaling_due_to_velocity), '/net/osiris/data/users/sanbus/PRACE4LOFAR/redshift_space/'+type+'/Cubes/neutral'+type+'_cube_rsd'+zzstr(i)+'.bin'
    
    ; Single=0: all, Single=1 Pd, =2 Px, =3 Pxd, =4 Cxd, =5 Cxdd, =6 Cxxd, =7 PdTb
    single=7
;     type='_'
    runps, density_array, neutral_array,1, zzstr(i),zz(i), loc_real,type+'fullreionization', single
;     runps, z_density_array, z_neutral_array,double(scaling_due_to_velocity),zzstr(i),zz(i), loc_reds, '_RSD_'+type, single
;     stop
;     if i eq nslices-1 then stop
endfor


end
     
