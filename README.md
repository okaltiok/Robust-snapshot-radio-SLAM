The projects contains an implementation of the robust snapshot radio SLAM algorithm presented in [1].

1. Download the project to your computer and run the main file which should produce the following result:

   LOS, estimates:  32/45, rmse: 0.2882 [m] 1.9456 [deg] 1.0554 [ns], time: 0.0307/0.3894 [ms] \n
   
   NLOS, estimates: 13/45, rmse: 0.4886 [m] 2.2702 [deg] 2.1263 [ns], time: 8.7723/0.3391 [ms]

   ALL, estimates:  45/45, rmse: 0.3578 [m] 2.0447 [deg] 1.4485 [ns], time: 2.5560/0.3749 [ms]

3. If Matlab throws an error from the mex-files, you need to either compile the mex-files on your computer or set "params.MEX = false". The mex-files can be compiled using CompileCLibraries.m found in the folder ".\mex files". Further instructions can be found in CompileCLibraries.m.     

References:

[1] Ossi Kaltiokallio, Elizaveta Rastorgueva-Foi, Jukka Talvitie, Yu Ge, Henk Wymeersch and Mikko Valkama, "Robust snapshot radio SLAM," accepted to IEEE Transactions on Vehicular Technology. [Online] Available: https://arxiv.org/abs/2404.10291
