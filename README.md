multi_cbs
=========

multiple peaks fit using root and  RooFit libraries

Usage:
---------
Two peaks example with different backgrounds and automatic peak start
```
.L ro_cbs.C+

ro_cbs( hh, 2, NULL, "pn")
ro_cbs( hh, 2, NULL, "p0")
ro_cbs( hh, 2, NULL, "p1")
ro_cbs( hh, 2, NULL, "p2")
```

--------
Polynomial background is Chebyshev for better convergation.

initial simnple guess based on equal distances or double* array of peak positions

Compile, else there is some error of FILE* log.

**Background**

pn ... no backgroung

p0 ... constant background

p1 ... linear bg

p4 ... maximum as now