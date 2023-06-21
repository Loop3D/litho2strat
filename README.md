# litho2strat

A statistical framework for estimation of local stratigraphy with uncertainty from drillhole lithology data.

Dependencies: numpy, matplotlib, networkx

To install:
```shell
pip install .
```

Usage:
```shell
python3 litho2strat.py -p <Parfile>
```

Example:
```shell
python3 litho2strat.py -p ./parfiles/Parfile_SA.txt
```

The module dependancies:  
(The graph is obtained by running: pydeps litho2strat.py --max-bacon 3 --exclude matplotlib numpy networkx)


![litho2strat](https://user-images.githubusercontent.com/33440028/214574327-12101efc-a4c2-49f0-a251-86cc60b9d5a2.svg)
