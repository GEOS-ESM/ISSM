source issm_env_withpy
rm -rf Models
rm -f GreenlandGEOS.bin GreenlandGEOS.outbin GreenlandGEOS.outlog GreenlandGEOS.errlog
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Greenland_input_control.py)
/discover/nobackup/agstubbl/ISSM/GEOS-ISSM/ISSM/bin/issm.exe StressbalanceSolution /discover/nobackup/agstubbl/ISSM/projs/IRF-ISSM GreenlandGEOS 2>> GreenlandGEOS.errlog
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Greenland_input_transient.py)
