pro readdens, density_file,pixels,DATA_TYPE=4,ENDIAN=LITTLE,extra_pixels=0
; Read density file
; This depends very much on the  exact structure of the density field

density_original = READ_BINARY(density_file,DATA_TYPE=DATA_TYPE,ENDIAN=ENDIAN)
;IDL data types 2: 16-bit INT, 3: 32-bit INT/LONG, 4: 32-bit FLOAT, 5: 64-bit FLOAT/DOUBLE
;https://www.harrisgeospatial.com/docs/idl_data_types.html

;Normalize density field if necesarry and reshape
if ABS(MEAN(density_original[extra_pixels:-1]))>0.001 THEN BEGIN
  density_array = REFORM(density_original[extra_pixels:-1]/MEAN(densityo[extra_pixels:-1])-1,pixels,pixels,pixels)
ENDIF
ELSE BEGIN
  density_array = REFORM(density_original[extra_pixels:-1],pixels,pixels,pixels
ENDELSE
delvar, density_original
return, density_original
end
