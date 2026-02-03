# litho2strat

A statistical framework for automated stratigraphic interpretation from drillhole lithology data with uncertainty quantification.

## Overview

Many legacy drillhole datasets contain only lithological descriptions without stratigraphic information, limiting their utility for regional 3D geological modeling. `litho2strat` addresses this challenge by:

- Systematically exploring all geologically plausible stratigraphic orderings for individual drillholes
- Correlating solutions across multiple drillholes to identify geologically consistent interpretations
- Combining constraints from lithological descriptions with stratigraphic relationships from regional maps
- Quantifying uncertainty by generating multiple plausible stratigraphic interpretations

This approach enables the transformation of lithology-only drillhole logs into stratigraphic sequences with quantified uncertainty, providing valuable constraints for subsurface modeling, resource estimation, and data acquisition planning.

## Reference

If you use this code, please cite:

Ogarko, V. and Jessell, M.: Automated stratigraphic interpretation from drillhole lithological descriptions with uncertainty quantification: litho2strat 1.0, Geosci. Model Dev., 19, 1007–1025, https://doi.org/10.5194/gmd-19-1007-2026, 2026.

This paper provides a detailed explanation of the methodology implemented in this code.

## Dependencies

numpy, matplotlib, networkx

## Installation
```shell
pip install .
```

## Usage
```shell
python3 litho2strat.py -p <Parfile>
```

### Example
```shell
python3 litho2strat.py -p ./parfiles/Parfile_SA.txt
```

## Module Dependencies

The graph is obtained by running: `pydeps litho2strat.py --max-bacon 3 --exclude matplotlib numpy networkx`

![litho2strat](https://github.com/Loop3D/litho2strat/assets/33440028/1327f3a2-ba94-4c19-b8ac-ba62832b5f26)
