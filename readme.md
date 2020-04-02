Fluctuations in Undulator Radiation

Add project directory to PYTHONPATH. For example, since I have the project folder `fur` in the home directory of my user `ilobach`, I added 
`export PYTHONPATH="${PYTHONPATH}:/home/ilobach/fur"`
to the end of my `.bashrc` file.

The python version has to be 3.6.5 for SRW to work (`pip install vinyl-srw`). The rest of the packages is given in `requirements.txt`.

To recreate the virtual environment with `virtualenv` do the following in the `\fur` directory:
`
$ python3 -m virtualenv env
$ source env/bin/activate
$ pip install -r requirements.txt
`

![Demo](demo.gif)
