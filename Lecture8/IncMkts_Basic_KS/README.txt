To run the Matlab codes you will first need to download COMPECON toolbox from Paul Fackler's website:

http://www4.ncsu.edu/~pfackler/compecon/toolbox.html

The files called in the code can either be added to the same directory as the code, or you can add the COMPECON directory to the 
Matlab path. Instructions are on their website.

The main file to run is start.m, there are numerous options that can be played around with.

The model is described in Examplemodel.pdf

One notable thing to check out is how much better the model does under an approximation of the continuous AR(1) process rather than the 
Rouwenhurst discretization. I was a little hazy on this during class since I hadn't really tried it out properly myself. 
Turns out its excellent. You can go bac and forth between AR1 and Rouwenhurst using options.AR1. 
With options.AR1='Y' you can change the order of approximation of the value function in the z-dimension to any order. 
With options.AR1='N', it must be of order 1.
The continuous approximation allows for almost exact simulation and the computation of impulse reponse functions which I've also included.

The current [default] options are currently set:
options.AC 	= 'Y'
options.ACdown 	= 'Y' 		% Removing this means you'll get kinks in the value function as firms will downwardsly adjust capital in 
				  jumps if capital is larger than optimal for that level of productivity
options.AR1 	= 'Y' 
spliorder 	= [3,3]; 	% Cubic splines in both dimensions
n 		= [10,5]; 	% 10 nodes for capital, 5 for productivity in the approximation of the value function

Obviously all of these can be changed.

Feel free to get in touch if you don't understand any aspect of the code / the notes.

- Simon Mongey