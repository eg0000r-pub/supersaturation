# supersaturation
This MATLAB program models vapor supersaturation in a laminar flow. It then uses the the solution to calculate the amount of vapor condensed on aerosol particles in the flow. Tube wall temperature profile is measured experimentally and provided in a csv file. The program then solves the heat PDE to find the temperature at all points in the flow:
<p align="center">
    <img src="images/plt1.png"/>
</p>
The program assumes that vapor at the inlet and hear the tube wall is saturated. Based on that, it solves the diffusion PDE to find vapor concentration at every point in the flow:
<p align="center">
    <img src="images/plt2.png"/>
</p>
The program then calculates vapor supersaturation, which is a function of temperature and vapor concentration:
<p align="center">
    <img src="images/plt3.png"/>
</p>
<p align="center">
    <img src="images/plt4.png"/>
</p>
Finally, the program solves an ODE to model vapor condensation on aerosol particles. The equation is solved at different radial positions and the solutions are averaged:
<p align="center">
    <img src="images/plt5.png"/>
</p>
