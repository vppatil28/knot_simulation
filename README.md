# knot_simulation

Code for simulating bend knots and 1-tangles (knots on one rope) using the Kirchhoff model for elastic fibers.

In this repository:
  
  one_tangles.m - code for simulating 1-tangles. Requires fiber.m, input_one_tangles.txt, the initial_conditions folder and a        results folder to be in the same directory in order to run.
  
  bend_knots.m - code for simulating bend_knots. Requires fiber.m, input_bend_knots.txt, the initial_conditions folder and a          results folder to be in the same directory in order to run.
  
  fiber.m - contains colormap used in one_tangles.m and bend_knots.m
  
  input_one_tangles.txt - plaintext file containing input parameters for one_tangles.m and instructions on what they mean.
  
  input_bend_knots.txt - plaintext file containing input parameters for bend_knots.m and instructions on what they mean.
  
  initial_conditions - folder containing .mat files which store the initial data for both bend_knots.m and one_tangles.m
  
  results - folder which will contain output files from running bend_knots.m and one_tangles.m


Instructions for using bend_knots.m

  1) Some choices of initial configurations are given in .mat files in the initial_conditions folder. Each .mat file containing a bend knot has two matrices, A (size 3 x n1) and B (size 3 x n2) representing the spatial positions of the two ropes in the bend knot. The files are named after the bend knot they represent. In most cases, two numbers are also included after the knot name, referring to which ends of the bend knot are pulled.
  
  2) A choice of bend knot, material parameters and pulling force can be specified by editing the top row of input_bend_knots.txt (more instructions are also in this text file). The bend knot can also be specified by directly editing bend_knots.m
  
  3) Running bend_knots.m will simulate the chosen bend, with the given material parameters and pulling force. The simulation has 10000 timesteps, corresponding to 1 second of real time. These parameters can also be edited within bend_knots.m if desired.
  
  4) Output from bend_knots.m will appear in the form of a video and a .mat file in the results folder. Output data consist of a 3 x n x 100 array for each rope, representing the spatial position the rope at 100 evenly spaced time intervals across the simulation. There is also an m x 100 array for each rope in the .mat file, representing the twist along it evolving in time. For both position and twist arrays, the last dimension is the time dimension. The twist array stores twist density multiplied by ds at each timepoint (ds is an arclength element of the rope). 
  
  5) The output video is colored by strain according to the colormap in fiber.m. (This colormap represents the strain - color relationship of a particular optomechanical fiber.)


Instructions for using one_tangles.m are essentially the same, except that there is only one rope being simulated.
