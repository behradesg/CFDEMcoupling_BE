# CFDEMcoupling

CFDEM®coupling stands for Computational Fluid Dynamics (CFD) - Discrete Element Method (DEM) coupling. It combines the open source packages OpenFOAM® (CFD) and LIGGGHTS® (DEM) to simulate particle-laden flows. CFDEM®coupling is part of the [CFDEM®project](https://www.cfdem.com).

[![CircleCI](https://circleci.com/gh/ParticulateFlow/CFDEMcoupling.svg?style=shield&circle-token=7e8118524babddbefccf4e3608a7545d405acbb4)](https://circleci.com/gh/ParticulateFlow/CFDEMcoupling)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

## Disclaimer

> This is an academic adaptation of the CFDEM®coupling software package, released by the
[Department of Particulate Flow Modelling at Johannes Kepler University in Linz, Austria.](https://www.jku.at/pfm)
> LIGGGHTS® and CFDEM® are registered trademarks, and this offering is not approved or
endorsed by DCS Computing GmbH, the official producer of the LIGGGHTS® and CFDEM®coupling software.
> This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.

## installation
Make sure OpenFOAM 6.0 is set up correctly and LIGGGHTS is installed as well. Clone the 
CFDEMcoupling source from the repository:

```bash
cd $HOME
mkdir CFDEM
cd CFDEM
git clone git@github.com:behradesg/CFDEMcoupling_BE.git
```

Open the bashrc file of CFDEMcoupling

```bash
gedit ~/CFDEM/CFDEMcoupling_BE/etc/bashrc &
```

Edit the lines marked as `USER EDITABLE PART` to reflect your installation paths correctly. Save the bashrc file and reload it:

```bash
source ~/CFDEM/CFDEMcoupling_BE/etc/bashrc
```

Entering $CFDEM_PROJECT_DIR in a the terminal should now give "... is a directory"

Check if everything is set up correctly:

```bash
cfdemSysTest
```

In case LIGGGHTS has already been compiled via cmake, it is possible to just compile LIGGGHTS-related sub-libraries using:
```bash
cfdemCompLIGlib
```

If the compilation fails with a message like

```bash
No rule to make target `/usr/lib/libpython2.7.so'
```

you probably need to create a symbolic link to the library in question.

Compile CFDEMcoupling (library, solvers and utilities) in one go

```bash
cfdemCompCFDEM
```

or alternatively step by step

```bash
cfdemCompCFDEMsrc
cfdemCompCFDEMsol
cfdemCompCFDEMuti
```

Find the log files of the compile process

```bash
cd ~/CFDEM/CFDEMcoupling_BE/etc/log
ls
```

If the file *log_compile_results_success* is present, compilation was successful.


## License

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

- This software is distributed under the [GNU General Public License](https://opensource.org/licenses/GPL-3.0).
- Copyright © 2009-     JKU Linz
- Copyright © 2012-2015 DCS Computing GmbH, Linz
- Some parts of CFDEM®coupling are based on OpenFOAM® and Copyright on these
  parts is held by the OpenFOAM® Foundation (www.openfoam.org)
  and potentially other parties.
- Some parts of CFDEM®coupling are contributed by other parties, which are
  holding the Copyright. This is listed in each file of the distribution.
