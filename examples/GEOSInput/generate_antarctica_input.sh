source issm_env_py
rm -f ISSM_AntarcticaGEOS.bin ISSM_AntarcticaGEOS.outbin ISSM_AntarcticaGEOS.outlog ISSM_AntarcticaGEOS.errlog
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Antarctica_input_control.py)
${ISSM_DIR}/bin/issm.exe StressbalanceSolution $(pwd) ISSM_AntarcticaGEOS 2>> ISSM_AntarcticaGEOS.errlog
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Antarctica_input_transient.py)
