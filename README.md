multi_cbs
=========

multiple peaks fit using root and  RooFit libraries

Usage:
---------
Two peaks example with different backgrounds and automatic peak start
```
ro_cbs( hh, 2, NULL, "pn")
ro_cbs( hh, 2, NULL, "p0")
ro_cbs( hh, 2, NULL, "p1")
ro_cbs( hh, 2, NULL, "p2")
```