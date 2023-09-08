# Badaletal2023
MATLAB 2022b code:
1) parMacrosimulation.m (main file)
2) macro_function.m (MATLAB function)
3) imageOut.m (MATLAB function)

This is a 2D model created to study the coexistence of yeast and bacteria on a plate-based assay. These MATLAB scripts of the model were submitted along with the paper titled "A dynamic fluid landscape mediates the spread of bacteria", written by Divakar Badal(1), Aloke Kumar(2), Varsha Singh(1,3), and Danny Raj M(4);(1) Department of Bioengineering, Indian Institute of Science, Bengaluru, India; (2)Department of Mechanical Engineering, Indian Institute of Science, Bengaluru, India; (3)Department of Developmental Biology and Genetics, Indian Institute of Science, Bengaluru, India; (4)Department of Chemical Engineering, Indian Institute of Science, Bengaluru, India. 

We created a grid of 100 x 100, representing a 10 x 10 cm area. We randomly seeded yeast agents (area of each agent = 6.3617e-07 cm2) on the grid and placed bacteria agents (area of each agent = 5.8905e-09 cm2) on the central pixel. We inoculated yeast agents 9.4 times more than the number of bacteria agents. The carrying capacity of a pixel on the grid is a 0.01 cm^2 area. Therefore, each pixel can hold several yeast and/or bacteria agents. We allowed the agents to grow at a certain growth rate and consume nutrients. Thus, in a pixel, the number of agents can easily cross carrying capacity, which leads to overflow in the immediate neighbor such that the nutrient-rich neighbor will have a higher probability of getting these cells. 


Parameters: 

Yeast division time = 4500 sec

Bacteria division time = 2100 sec

Fluid effect constant = 0.9e-5

dt = 1/3600 hours

Fluid pumping rate = 0.3e-9 cm3/sec/cell

Diffusivity of nutrient = 1.24e-5 cm2/sec

Diffusivity of fluid = 2.3e-5 cm2/sec

Bacteria nutrient consumption rate = 1.94e-12 m/cell/sec

Yeast nutrient consumption rate = 4.4e-12 m/cell/sec

area of bacteria agent = 5.8905e-09 cm2

area of yeast agent = 6.3617e-07 cm2

Simulation time = 6 hours


Input parameters:

Number of realizations (1-1000)

Fluid Pumping Rate Factor (0.5-1.5)

Agent count Multiplier (0.5-1.5)

Growth rate Multiplier (0.5-1.5)

imFlag (1 or 0)


Exporting each frame

if ‘imFlag’ is 1, then script will create folder named 'Frames' at current directory. Inside 'Frames' folder, there will be three folders namely 'FluidLayer', 'CoExit', and 'YeastLayer'. Further, inside each of these folders there will be three folders representing without yeast, with yeast with fluid pool, and with yeast and without fluid pool named as '0', '1', and '2' folder names respectively. 
