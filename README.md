# 2-BESS
Code operation starts from the main function and is divided into five tasks.  
## Task 1 modulation of signal
Generate voltage levels and corresponding time intervals based on the number of batteries and battery voltage.
## Task 2 user input
Set system and environment parameters.
## Task 3 connection design(layer 1 converter & battery use)
Using distribution flattening method to get Expected Set of Batteries(the battery power（current）capability)，then sovle MILP problem to select the best layer 1 converter connections and the best battery connections.   
contain function: func_connection_design , func_energyflow_connection
## Task 4 MC simulation
Fix batteries and layer 1 converter connections, MC simulations are performed to obtain the average Up and system efficiency  
contain function: func_MC , func_energyflow_MC
## Task 5 Output PQ
Calculate the system output when the load is not purely resistive.  
contain function: func_PQ , func_energyflow_PQ
