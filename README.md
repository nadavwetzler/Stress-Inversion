# Stress-Inversion
A Matlab code to calculate stress field from FMS data over a range of friction coefficients 

# Features
The code calculates the principal axis of the stress state from focal mechannism information
Use FMS data at /data folder as a tamplate for your one dataset

To identify the fault plane from the two nodal planes, determination of the faults frictional strength and calculate the local stress field we used stress-inversion approach developed to calculate the stress state associated with a set of focal-mechanisms (Reches et al., 1992; Busetti and Reches, 2014). This approach was customized into MATLAB environment to calculate the orientations and relative magnitudes of the three principal stress axes (3D stress tensor) for a group of FMSs under three assumptions: 
1) All the earthquakes in the group occurred under the same stress state. 
2) Slip along a fault occurs in the direction of maximum resolved shear stress 
3) The shear and normal stress on the faults satisfy the Coulomb failure criterion.

Each pair of the nodal planes is tested with respect to the Principle Axes Misfit Angle (PAM), which is the angle between the ideal stress axes of each nodal plane and general stress axes of the entire group according to the optimal mechanical condition for faulting. The quality of the calculated stress tensor is represented by the confidence levels, calculated by bootstrapping method, for 500 random samples of the original FMS group. 
Stress is inverted for each friction coefficient that ranges 0.1 to 0.8, weighted by the earthquake magnitude, and fault planes are selected according to smallest PAM between the two FMS nodal planes. To ensure coherent selection, two more criteria are set: 
1) PAM is smaller than 30 deg. 
2) The aperture between the PAM of the two nodal planes is larger than 10%.

Busetti, S., Jiao, W., Reches, Z., 2014. Geomechanics of hydraulic fracturing microseismicity: Part 1. Shear, hybrid, and tensile events. Am. Assoc. Pet. Geol. Bull. 98, 2439–2457. https://doi.org/10.1306/05141413123

Reches, Z., Baer, G., Hatzor, Y., 1992. Constraints on the strength of the upper crust from stress inversion of fault slip data. J. Geophys. Res. 97, 12481. https://doi.org/10.1029/90JB02258

# Run
Download the files

Add the contant of the folder to your Matalb path

Run run_Stress_inv.m file and you should soon see the results

Depending on the resolution of your screen, the lagend may overlap the bottom of subplot. you can shift it down
