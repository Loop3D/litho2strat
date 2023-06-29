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

![litho2strat](https://github.com/Loop3D/litho2strat/assets/33440028/1327f3a2-ba94-4c19-b8ac-ba62832b5f26)
