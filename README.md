# hwo_pol
Polarization aberration simulation in support of NASA's Habitable Worlds Observatory

Tentatively referred to as `hwop` for short, on behalf on Max Millar-Blanchaer

## Installation
To install `hwo_pol` we assume you have a working Python 3.8+ environment,
with anaconda and pip installed.


0. Install Poke
If you want to create new Jones pupils without using the supplied Rayfronts,
you will need to install Poke. Follow the installation instructions on the 
[Poke repository](https://github.com/Jashcraf/poke).

1. Clone the repository
```bash
git clone https://github.com/Jashcraf/hwo_pol.git
```

2. Install the package
This will install the package in development mode, so any changes you make to
the code will be reflected in the installed package. This is important for 
including the coating data.

```bash
cd hwo_pol
pip install -e .
```

3. Populate the coating folder
The coating recipes are stored in the `coating_recipes` folder at the top level
of the repository. This data is in part protected or otherwise not mine to 
distribute, but it can be provided on request to Jaren Ashcraft. Please open
an issue on the repository if you need this data.

Once given the data, copy the files from the folder into the `coating_recipes`
folder in this directory.

4. Test the installation
For a simple Jones pupil computation, run the `test_jones_pupil.py` script in
the `scripts` folder.
