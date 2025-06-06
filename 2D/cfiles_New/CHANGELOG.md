# Changelog

## 2025-06-06
- Renamed `NURBSbasis.c` to `NURBSbasisOld.c` for archival purposes.
- Added new file `NURBSbasis.c`:
  - Supports multiple basis function indices (`i`) as input.
  - Returns arrays of NURBS basis values and their derivatives at a given parametric point.
- Added new file `NURBSinterpolation3d.c`.
