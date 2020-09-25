%problems to be solved:

%1. hydraulic concuctivity, use different formulations fro Richards Equation and bucket W
%2. allow seepage flow downwards for uppermost unconfined aquifer downwards,
%even if aquifers don't overlap - add this to the overlap-function!
%3. improve routing of water in the push-step, so that advection works properly - difficult one, probably needs to make the nomal routing (switching flow between cells) for water flowing down!  
%4. use the correct diatncealso for unconfined aquifers
%5. add other classes like snow and lakes
%6. Difficult: check if a class with index 1 is filled up and then recalculate the heads, partition timestep for different phases
%7. Implement for excess ice
%8. Implement for Richards Eq
%9. Implement a fucntionality that detects the ground surface and matches
%the elevations at the ground surface