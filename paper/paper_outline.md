# Paper Outline 

## Introduction: 
### Motivation: 
Ocean carbon and heat uptake and storage are very important processes in regulating the s climate. As shown in Frolicher et al. (2015) there is considerable uncertainty in carbon and heat uptake and storage as represented by CMIP5 models, especially in the Southern Ocean. A further complicating factor is open ocean deep convection in the Weddell Sea. Often referred to as the Weddell polynya, these deep convective events are important for the formation of Antarctic Bottom water and the global heat and carbon budget (Killworth, 1983). Climate models simulate deep convection oscillations to varying degrees and frequencies (de Lavergne et al., 2014); however, with few observational constraints it is difficult to determine how well modeling groups represent Southern Ocean dynamics.

### Previous Papers using Aredi Simulation Suite: 
***Gnanadesikan et al, 2015: Isopycnal mixing by mesoscale eddies significantly impacts oceanic anthropogenic carbon uptake.***

Results:

Results from historical-ish simulations: 
Aredi 400 mixing simulation results in the least amount of anthropogenic carbon flux by the ocean (2.4 Gt/yr in final year → year 2000?). Aredi 2400 mixing simulation results in the most anthropogenic carbon flux (2.8 Gt/yr). Larger difference seen in Southern Ocean (Figures 2a and b). Aberzonal and Aber2D lateral mixing estimates lie between these two extreme simulations. Low lateral mixing coefficient in Southern Ocean is more important than the high values in the gyre interior with respect to carbon flux. 

Consistent with the carbon flux, the total inventory of anthropogenic carbon (at year 1995) is largest for the Aredi 2400 simulation (118.7 GtC) and smallest for the Aredi 400 simulation (102 GtC) with the Abernathy estimates falling between the two. This range is about half of the variability across the CMIP5 models, with these estimates being at the higher end.  

Doubling CO2 simulations:
Including the circulation changes due to changes is radiative forcing in the instantaneous doubling CO2 simulation, the spatial pattern of carbon storage remains the same as the biological-only historical simulation. In the Aredi 400 mixing simulation, the primary region of anthropogenic carbon storage is in the N. Atlantic. As the mixing parameter increases, more anthropogenic carbon is stored in the N. Pacific. This increased storage is due to the higher mixing destabilizing the halocline and allows old waters to come up to the surface and take up anthropogenic carbon (NOT SHOWN).

### Aim:


## Methods: 
### Aredi Simulations:  
We use the GFDL ESM2Mc, a coarse resolution configuration of the GFDL ESM2Mc. We use different parameterizations for the lateral mixing coefficient, Aredi (Redi, 1982), and mesoscale eddy advenction, Gent-McWilliams (Gent & McWilliams, 1990), to asses the impact on the heat and carbon content. Here we analyze 4 pre-industrial control simulations: 
* Three simulations where Aredi is constant at 400, 800, and 2400 m2  
* GM minimum value is increased from 200 to 600 m2s-1 (Aredi = 800 m2s-1)


## Current Results: 
1. Simply changing the along-isopycnal mixing parameterization significantly changes the convective variability in the GCM. Additionally, by changing the GM miximum value, convection ceases in the Weddell Sea. 

<img align="right" img src="paper_outline_figures/aredi_heat_depth.png" alt="alt text" width="400">
Annually-averaged Weddell Sea subsurface temperature for different mixing simulations. Black line indicates average MLD. Low-mixing simulation (Aredi = 400 m2s-1) show periodic deep convection events in the Weddell Sea with periods of no convection where subsurface heat build up. High mixing simulation (Aredi = 2400 m2s-1) shows constant convection with no subsurface heat build up.



2. While raising the GM minimum value stops Weddell Sea convection, there is not a change in the average ocean circulation over the period of the model run. Therefore changing the parameterization changes convection, but not the overall dynamics. 

## To Do List: 
