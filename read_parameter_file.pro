PRO read_parameter_file, file_name, parameters

line = ''

OPENR, lun, file_name, /GET_LUN

WHILE NOT EOF(lun) DO BEGIN 
    READF, lun, line
    IF STRMID(line,0,1) NE '#' AND STRMID(line,0,1) NE ' ' THEN BEGIN
        void = EXECUTE('parameters.'+line)
    ENDIF
ENDWHILE
free_lun, lun

END
