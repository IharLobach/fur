Fluctuations in Undulator Radiation

See [https://github.com/IharLobach/ursse](https://github.com/IharLobach/ursse) for experiments with a single electron.


Add project directory to PYTHONPATH. For example, since I have the project folder `fur` in the home directory of my user `ilobach`, I added 
`export PYTHONPATH="${PYTHONPATH}:/home/ilobach/fur"`
to the end of my `.bashrc` file.

The python version has to be 3.6.5 for SRW to work (`pip install vinyl-srw`). The rest of the packages is given in `requirements.txt`.

To recreate the virtual environment with `virtualenv` do the following in the `\fur` directory:
```
$ python3 -m virtualenv env --python=python3.6.5
$ source env/bin/activate
$ pip install -r requirements.txt
```

Some plots use latex for text rendering. Use

```
sudo apt install texlive-full
```

to install latex on your machine, or set option `usetex=False` in the plots.

Some of the figures require `arial.ttf` font to be installed on the ubuntu machine. It can be achieved by running
```shell
sudo apt-get install --reinstall ttf-mscorefonts-installer
```

![Demo](demo.gif)
