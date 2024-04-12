# Closing the Service: Contrasting activity-based and time-based systematic closure policies
## Antonio Castellanos, Andrew Daw, Amy Ward, Galit B. Yom-Tov

This repository contains the source code and supplementary materials for the research paper titled "Closing the Service: Contrasting activity-based and time-based systematic closure policies,".

### Authors and Affiliations
- **Antonio Castellanos**
  - Booth School of Business, The University of Chicago, Chicago, IL, USA
- **Andrew Daw**
  - Marshall School of Business, University of Southern California, Los Angeles, CA, USA
- **Amy Ward**
  - Booth School of Business, The University of Chicago, Chicago, IL, USA
- **Galit B. Yom-Tov**
  - Faculty of Data and Decision Sciences, Technion—Israel Institute of Technology, Haifa, Israel

### Contact
For any additional questions or feedback, please contact:
- Antonio Castellanos: [Antonio.Castellanos@chicagobooth.edu](mailto:Antono.Castellanos@chicagobooth.edu)
- Andrew Daw: [dawandre@usc.edu](mailto:dawandre@usc.edu)
- Amy Ward: [Amy.Ward@chicagobooth.edu](mailto:Amy.Ward@chicagobooth.edu)
- Galit B. Yom-Tov: [gality@technion.ac.il](mailto:gality@technion.ac.il)



### Abstract
We examine different policies for systematic service closure in messaging service systems. The system is modeled as an $M/UHP/1$ queue, where service times follow a history-based Hawkes cluster process. We propose and examine stopping-time rules that balance between queue length and the probability of prematurely closing conversations. In a  simulation study, we compare two families of systematic closure policies: the first relies on predictive information regarding service progress, i.e., the conversation's activity levels, while the second relies on elapsed time without activity. When restricted to static threshold policies, both families provide similar performance. However, when allowing the threshold to vary with the system state, activity-level policies outperform the inactive-time policies. 
Moreover, a large difference is observed between static and dynamic threshold policies. We therefore conclude that state-dependent (i.e., dynamic) activity-based policy is the most promising candidate to achieve optimal closure rules.

### Repository Contents
- `Figure1_Hawkes_Intensity_Maker.m`: UHP simulation for Figure 1, and figure maker. 
- `Runner_Static_Policies_and_Figure2_Maker.m`: Static policies runner and Figure 2 maker. Requires: `tauFuncSim.m`, `uhpAbsorbCheckRep.m`, `uhpAbsorbCheckSim.m`, `randBorel.m`, `randPF3.m`.
- `Runner_Dynamic_Policies_and_Figure3_Maker.m`: Dynamic policies runner and Figure 3 maker. Requires: `mHawkes1Rep_Online.m` and `mHawkes1Sim_Online.m`.

### Usage
Instructions on how to setup and run the code: run the 3 previous code scripts, they simulate, create figures and save them. 


