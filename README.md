The project contains an implementation of the robust snapshot radio SLAM algorithm presented in [1].

1. Download the project to your computer and run the main file which should produce the following result:

   LOS, estimates: 32/45, rmse: 0.2882 [m] 1.9456 [deg] 1.0554 [ns], time: 0.0103/0.0853 [ms] <br />
   NLOS, estimates: 13/45, rmse: 0.4886 [m] 2.2702 [deg] 2.1263 [ns], time: 5.9994/0.1686 [ms] <br />
   ALL, estimates: 45/45, rmse: 0.3578 [m] 2.0447 [deg] 1.4485 [ns], time: 1.7405/0.1093 [ms] <br />

3. If Matlab throws an error from the mex-files, you need to either compile the mex-files on your computer or set "params.MEX = false" to use the m-file implementation of the algorithm. The mex-files can be compiled using CompileCLibraries.m found in the folder ".\mex files". Further instructions can be found in CompileCLibraries.m.     

References:

[1] Ossi Kaltiokallio, Elizaveta Rastorgueva-Foi, Jukka Talvitie, Yu Ge, Henk Wymeersch and Mikko Valkama, "Robust snapshot radio SLAM," accepted to IEEE Transactions on Vehicular Technology. [Online] Available: https://arxiv.org/abs/2404.10291
