# Parcel Simulator 1D
### The numerical model of convective air parcel ascent.
![Github Actions](https://github.com/Quba1/ParcelSimulator1D/workflows/g++%20build/badge.svg)

The easiest way to compile the model is to use a Linux machine and `g++`.
To do this run those commands in a suitable directory:
```bash
git clone https://github.com/Quba1/ParcelSimulator1D.git
cd ParcelSimulator1D
make all
```
Now edit `model.conf` and `parcel.conf` in `config` directory to setup your simulation, and run the model:
```bash
./simulator.exe
```
The model output is located in `./output` directory.

You can also use your own input file. Simply copy sample profile in `input` directory and modify it with your own values.

To remove all created executables run:
```bash
make clean
```
