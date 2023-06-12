# Files

`01_init_data.py` : create a data json file with the basic data
`02_update_data.py` : populate the json file with missing parameters that are computed automatically
`03_makexml.py` : build the xml file for geos following the json data file
`04_lock_files.sh` : [optional] set .xml and .json file to readonly to avoir possible later confusion
`05_cutHdf5.py`: actually not working. Cut the hdf5 into numpy array
`06_plotres.py`: plot the numpy arrays provided by `05_cutHdf5.py`
`notebook_viewing_wavefield.ipynb`: Jupyter notebook for post processing

# How to use ?

1. `python3 01_init_data.py` to init the `data.json` file.
2. Modify the json file accordingly
3. `python3 02_update_data.py` to init the `data.json` file.
4. Check the json file and modify what you would want. Be careful though what you change...
5. `python3 03_makexml.py`: build the `vti.xml` file for `GEOS`
6. `. ./lock_files.sh`: avoid further confusion by locking the files (read only)
7. `geos vti.xml`: run the simulation
8. `jupyter notebook --no-browser --port=8891 --ip=0.0.0.0`: launch the notebook
