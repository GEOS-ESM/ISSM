source issm_env_py
rm -f ISSM_GreenlandGEOS.bin ISSM_GreenlandGEOS.outbin ISSM_GreenlandGEOS.outlog ISSM_GreenlandGEOS.errlog
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Greenland_input_control.py)
${ISSM_DIR}/bin/issm.exe StressbalanceSolution $(pwd) ISSM_GreenlandGEOS 2>> ISSM_GreenlandGEOS.errlog
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Greenland_input_transient.py)
