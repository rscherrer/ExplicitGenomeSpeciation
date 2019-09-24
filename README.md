# ExplicitGenomeSpeciation

Branch|[![Travis CI logo](ci_setup/pics/TravisCI.png)](https://travis-ci.org)|[![Codecov logo](ci_setup/pics/Codecov.png)](https://www.codecov.io)
---|---|---
master|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=master)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=master)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/master)
develop|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=develop)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=develop)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/develop)

An individual based simulation of adaptive speciation with explicit genome, additive and non-additive genetic effects.

# Compile

module load Qt5
qmake EGS.pro
make --silent release

# Run

```bash
./EGS
```
or
```bash
./EGS parameters.txt
```

# Submit to Peregrine

```bash
./launcher.sh

# Install Python dependencies

```bash
./pystaller.sh # if needed, to install Python dependencies
```

# Read single files

```bash
./reader.py <variable_name>.dat
```

or

```bash
python reader.py <variable_name>.dat
```

# Plot single files

```bash
./plotter.py <variable_name>.dat
```

or

```bash
python plotter.py <variable_name>.dat
```


