@echo off
python motif_search_main.py ^
--input_type durbin_txt ^
--input_network_file "networks/data/Durbin_1986/neurodata.txt" ^
--synapse_threshold 25 ^
--filter_sex_type herm ^
--filter_syn_type chem ^
--sub_graph_algorithm mfinder_i ^
--k 3 ^
--sim 1 ^
--bin_file "durbin_herm_chem_k4_m5.bin" ^
--run_motif_criteria

