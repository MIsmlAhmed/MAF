# summaflow: create summa input files from the agnostic component.

prerequisites:

1. Model agnositc ouptus
2. MESH forcing nc (can be obtained by running the MESH flow notebook)

Assumptions:

1. This notebook assumes one HRU per GRU (HRU = GRU = subbasin). This is commonly used in SUMMA for large-scale simulations.
2. This notebook is based on using NALCMS and SoilGrid inputs so that it converts their indexes to equivalent ones by USGS landuse and ROSETTA soil classes.