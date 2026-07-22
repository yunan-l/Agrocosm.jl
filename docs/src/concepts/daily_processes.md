# Daily process order

Process order is part of the model definition because same-day water, carbon,
nitrogen, and heat exchanges are coupled.

The current C3/C4 drivers follow this high-level sequence:

1. Read daily climate and update rolling climate history.
2. Apply cultivation and scheduled fertilizer/manure events.
3. Apply tillage and bioturbation.
4. Compute surface albedo and potential evaporation/radiation variables.
5. Update snow, hydraulic properties, litter properties, and soil temperature.
6. Decompose existing litter/SOM and update mineral nitrogen.
7. Advance crop phenology and harvest; route new residues for future days.
8. Apply canopy and litter interception, infiltration, and percolation.
9. Compute absorbed radiation, temperature stress, photosynthesis, and
   water-limited conductance.
10. Update crop carbon, nitrogen, soil evaporation, and transpiration removal.
11. Apply denitrification, volatilization, and leaching losses.
12. Record conservation diagnostics and user outputs.

New harvest residues enter after the day's decomposition and become eligible
on the following day. Mineralization occurs before crop uptake, so newly
mineralized nitrogen is available on the same day.

Changing this order can change the trajectory even when individual equations
are unchanged. C3/C4 source-order tests protect the audited contract.
