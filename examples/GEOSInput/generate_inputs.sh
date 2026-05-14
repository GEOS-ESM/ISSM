./generate_greenland_input.sh
./generate_antarctica_input.sh
source issm_env
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python global_mesh.py)
