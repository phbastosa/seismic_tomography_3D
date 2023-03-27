#!/bin/bash

# python3 SimBGf_2022_models.py

# cd ../run/

# ./run -build_all

# cd ../SimBGf_2022

# # Podvin & Lecomte (1991)
# ./../bin/eikonal.exe parameters/central_100m_classic.txt
# ./../bin/eikonal.exe parameters/externs_100m_classic.txt
# ./../bin/eikonal.exe parameters/central_50m_classic.txt
# ./../bin/eikonal.exe parameters/externs_50m_classic.txt
# ./../bin/eikonal.exe parameters/central_25m_classic.txt
# ./../bin/eikonal.exe parameters/externs_25m_classic.txt

# # Jeong & Whitaker (2008) 
# ./../bin/eikonal.exe parameters/central_100m_fim.txt
# ./../bin/eikonal.exe parameters/externs_100m_fim.txt
# ./../bin/eikonal.exe parameters/central_50m_fim.txt
# ./../bin/eikonal.exe parameters/externs_50m_fim.txt
# ./../bin/eikonal.exe parameters/central_25m_fim.txt
# ./../bin/eikonal.exe parameters/externs_25m_fim.txt

# #  Noble, Gesret and Belayouni (2014) 
# ./../bin/eikonal.exe parameters/central_100m_fsm.txt
# ./../bin/eikonal.exe parameters/externs_100m_fsm.txt
# ./../bin/eikonal.exe parameters/central_50m_fsm.txt
# ./../bin/eikonal.exe parameters/externs_50m_fsm.txt
# ./../bin/eikonal.exe parameters/central_25m_fsm.txt
# ./../bin/eikonal.exe parameters/externs_25m_fsm.txt

python3 SimBGf_2022_figures.py
