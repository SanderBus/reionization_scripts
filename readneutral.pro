pro readneutral, neutral_file, pixels, neutral=1, DATA_TYPE=5, ENDIAN=LITTLE, extra_pixels=0
; Load neutral array

if neutral THEN BEGIN
    ionized_array=READ_BINARY(neutral_file,DATA_TYPE=DATA_TYPE,ENDIAN=ENDIAN)
    neutralarray=REFORM(1-ionizedfrac[extra_pixels:-1],pixels,pixels,pixels)
    delvar, ionized_array
ENDIF
ELSE BEGIN
    neutralarray_original = READ_BINARY(neutral_file,DATA_TYPE=DATA_TYPE,ENDIAN=ENDIAN)
    neutralarray = =REFORM(neutralarray_original[extra_pixels:-1],pixels,pixels,pixel)
    DELVAR, neutralarray_original
ENDELSE

RETURN, neutralarray
end
