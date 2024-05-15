@echo off
for /L %%i in (5,1,8) do (
    python ..\motif_search_main.py ^
    --input_type worm_wiring_xlsx ^
    --input_network_file "../networks/data/Cook_2019/SI 2 Synapse adjacency matrices.xlsx" ^
    --sheet_name "herm chem synapse adjacency" ^
    --synapse_threshold %%i ^
    --bin_file "../results/synapse_threshold_analysis/example_k3_m%%i.bin "^
    --sub_graph_algorithm triadic_census ^
    --run_motif_criteria
)
