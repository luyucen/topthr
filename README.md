# topthr
A thresholding method for the minimal compliance problem which uses characteristic functions as the design variable. The method is unconditionally energy stable with proven energy decaying property. 
The finite element implementation is based on the 88 lines of MATLAB code written by E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov and O. Sigmund [1].
This folder contains tests on three examples with different boundary conditions and loads:\
Example 1. Traction force applies on the right bottom corner of a cantilever beam with its left boundary fixed;\
Example 2. Point force applies on the middle of the right boundary of a cantilever beam with its left boundary fixed;\
Example 3. Half of an MBB beam.\
Two continuation schemes have been used:
1. Continuation based on decreasing Emin. For example, run `demo_example1_200X100` in MATLAB, to apply this continuation method on Example 1 with mesh resolution 200X100;
2. Continuation based on decreasing sensitivity filter. For example, run `demo_sensitivity_continuation_example1` to apply this continuation method on Example 1 with mesh resolution 100X50, 200X100, 300X150 and 400X200.

# Reference
[1] Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B. S., & Sigmund, O. (2011). Efficient topology optimization in MATLAB using 88 lines of code. Structural and Multidisciplinary Optimization, 43, 1-16.
