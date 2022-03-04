# pointPositioning

The pointPositioning project aims to provide the user with modules containing python functions for processing GPS and GALILEO data.

Goals:

1.  Raw data reader
    Allow the user to read data from NAV and OBS RINEX, both 2.x and 3.x and including Hatanaka compressed OBS and SP3
    and store them in DataSet structures suitable for analysis and plotting.
    
2.  Computation of satellites’ Earth-fixed position, velocity and clock offset

3.  Definitions of ad hoc functions for the estimation of atmospheric disturbances, in particular:
    -   Tropospheric delay is computed with Saastamoinen model.
    -   Ionospheric delay is computed with Klobuchar mdel for GPS observations and NeQuick G algorithm for GALILEO observations

4.  Point positioning in single epoch with Least Squares on code observations
    Estimate of the position [X, Y, Z] of a receiver and its clock offsets with respect to GPS and GALILEO System Time
    
5.  Transformation of coordinates between different reference systems

The functions are organized in modules and are contained in the codepos folder, except for the implementation of the NequikG algorithm which is contained in NQ folder.

An example of use of the software is provided in TUTORIAL-NOTEBOOK jupyter notebook.

Required libraries:

-   pandas
-   numpy
-   math
-   matplot
-   geoRinex
