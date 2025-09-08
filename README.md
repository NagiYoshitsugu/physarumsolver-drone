# Physarum Solver for drone flight paths
## Overview
This project implements a method that uses Physarum Solver to optimize drone flight paths based on radio signal strength.

## Input file
edges.csv  
Required columns:  
|Column	|Type	|Description|
|----|----|----|
|source	|int	|Start node index (0–63)|
|target	|int	|End node index (0–63)|
|length	|float	|Edge length (positive number)|
|D_max	|float	|Maximum allowable conductivity for the edge|

Notes:  
- Node indices must match a 4×4×4 cubic grid layout.  
- Lengths may be min–max normalized if desired for scaling purposes.  

## Output files
#### 1. physarum_results.csv  
Final conductivities of all edges after simulation.  
Columns:  
|Column |Description|
|----|----|
|source |Start node index|
|target |End node index|
|D_final |Final conductivity after iterations|
|D_max |Maximum allowed conductivity|

#### 2. physarum_path.csv
The path from start to end that carries the maximum flow.  
Columns:  
|Column|	Description|
|----|----|
|source|	Start node index of edge in path|
|target|	End node index of edge in path|
|Q|	Flow value through the edge|


#### 3. 3D Visualization PNGs
- physarum_network.png  
Full network, edges colored by final conductivity D_final.  
- physarum_path.png  
Maximum flow path in red; other edges are light gray.  
- physarum_dmax.png  
Network edges colored by D_max.  

## Example
|signal strength|flight path|
|----|----|
|<img width="1200" height="1678" alt="signal strength" title="signal strength" src="https://github.com/user-attachments/assets/e4d6bb06-b108-4b3c-a78a-9846e0912e7f" />|<img width="1000" height="1907" alt="physarum_path" title="flight path" src="https://github.com/user-attachments/assets/57e91bee-445e-4df5-a016-8328eeca9ef4" />|



## References
Tero, A., Kobayashi, R., & Nakagaki, T. (2006). A mathematical model for adaptive transport network in path finding by true slime mold. Journal of Theoretical Biology, 244(4), 553–564.  
