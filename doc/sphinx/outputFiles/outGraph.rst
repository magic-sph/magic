.. _secGraphFile:

Graphic files ``G_#.TAG``
=========================

These are fortran unformatted files containing 3D data (in the form vector_array(phi, theta, r) ) which can be used to visualize the solution. They are written in chunks of latitude for one radial level at a time by the subroutine :f:subr:`graphOut <graphout_mod/graphout>` or by :f:subr:`graphOut_mpi <graphout_mod/graphout_mpi>` depending on whether ``USE_MPI`` is set to ``Yes`` or ``No`` in the Makefile. The structure of the file looks like below:

  .. code-block:: fortran

      !-------------
      ! Header
      !-------------

      version

      runid

      time, n_r_max, n_theta_max, n_phi_tot, n_r_ic_max-1, minc, nThetasBs, ra, ek, pr, prmag, radratio, sigma_ratio

      theta(1:n_theta_max)

      !-------------------------------------------------------

      !-------------------
      !Graphout_version_9
      !-------------------
      
      !This version is written when the code uses MPI (USE_MPI=yes). Parallel chunks of fields are written for different radial levels. Chunks in thete are written in parallel using OpenMP

      !---------------
      ! Data
      !---------------
      
      !-------------
      ! Block N
      !-------------

      n_r-1, r(n_r)/r(1), n_theta_start, n_theta_stop           !Radial index, radius in terms of r_cmb, start and stop of the theta block

      sr(1:n_phi_tot, n_theta_start:n_theta_stop, n_r)          !Entropy

      vr(1:n_phi_tot, n_theta_start:n_theta_stop, n_r)          !Radial velocity

      vt(1:n_phi_tot, n_theta_start:n_theta_stop, n_r)          !Theta component of velocity

      vp(1:n_phi_tot, n_theta_start:n_theta_stop, n_r)          !Zonal (phi component) of velocity

      if(l_mag):						!For a magnetic run

         br(1:n_phi_tot, n_theta_start:n_theta_stop, n_r)       !Radial magnetic field

	 bt(1:n_phi_tot, n_theta_start:n_theta_stop, n_r)       !Theta component of magnetic field

	 bp(1:n_phi_tot, n_theta_start:n_theta_stop, n_r)       !Zonal (phi component) of magnetic field
      

      !-------------------
      !Graphout_version_7
      !-------------------

      !This version is written when the code does not use MPI (USE_MPI=no). Chunks in theta are written in parallel with OpenMP.

      !---------------
      ! Data
      !---------------
      
      !-------------
      ! Block N
      !-------------

      n_r-1, r(n_r)/r(1), n_theta_start, n_theta_stop

      !-----------------
      ! Entropy
      !-----------------

      sr(1,n_theta_start,n_r)              !n_phi = 1, n_theta = n_theta_start, n_r
      sr(2,n_theta_start,n_r)	           !n_phi = 2, n_theta = n_theta_start, n_r
      ...
      sr(n_phi_tot,n_theta_start,n_r)      !n_phi = n_phi_tot, n_theta = n_theta_start, n_r
      sr(1,n_theta_start+1,n_r)            !n_phi = 1, n_theta = n_theta_start+1, n_r
      ...
      sr(n_phi_tot,n_theta_start+1,n_r)
      ...
      sr(1,n_theta_stop,n_r)               !n_phi = 1, n_theta = n_theta_stop, n_r
      sr(2,n_theta_stop,n_r)               !n_phi = 2, n_theta = n_theta_stop, n_r
      ...
      sr(n_phi_tot,n_theta_stop,n_r)       !n_phi = n_phi_tot, n_theta = n_theta_stop, n_r 
      
      !-----------------
      ! Radial velocity
      !-----------------

      vr(1,n_theta_start,n_r)              !n_phi = 1, n_theta = n_theta_start, n_r
      vr(2,n_theta_start,n_r)	           !n_phi = 2, n_theta = n_theta_start, n_r
      ...
      vr(n_phi_tot,n_theta_start,n_r)      !n_phi = n_phi_tot, n_theta = n_theta_start, n_r
      vr(1,n_theta_start+1,n_r)            !n_phi = 1, n_theta = n_theta_start+1, n_r
      ...
      vr(n_phi_tot,n_theta_start+1,n_r)
      ...
      vr(1,n_theta_stop,n_r)               !n_phi = 1, n_theta = n_theta_stop, n_r
      vr(2,n_theta_stop,n_r)               !n_phi = 2, n_theta = n_theta_stop, n_r
      ...
      vr(n_phi_tot,n_theta_stop,n_r)       !n_phi = n_phi_tot, n_theta = n_theta_stop, n_r 

      !-----------------------------
      ! Theta component of velocity
      !-----------------------------

      vt(1,n_theta_start,n_r)              !n_phi = 1, n_theta = n_theta_start, n_r
      vt(2,n_theta_start,n_r)	           !n_phi = 2, n_theta = n_theta_start, n_r
      ...
      vt(n_phi_tot,n_theta_start,n_r)      !n_phi = n_phi_tot, n_theta = n_theta_start, n_r
      vt(1,n_theta_start+1,n_r)            !n_phi = 1, n_theta = n_theta_start+1, n_r
      ...
      vt(n_phi_tot,n_theta_start+1,n_r)    
      ...
      vt(1,n_theta_stop,n_r)               !n_phi = 1, n_theta = n_theta_stop, n_r
      vt(2,n_theta_stop,n_r)               !n_phi = 2, n_theta = n_theta_stop, n_r
      ...
      vt(n_phi_tot,n_theta_stop,n_r)       !n_phi = n_phi_tot, n_theta = n_theta_stop, n_r 


      !----------------------------------
      ! Zonal (phi component) of velocity
      !----------------------------------

      vp(1,n_theta_start,n_r)              !n_phi = 1, n_theta = n_theta_start, n_r
      vp(2,n_theta_start,n_r)	           !n_phi = 2, n_theta = n_theta_start, n_r
      ...
      vp(n_phi_tot,n_theta_start,n_r)      !n_phi = n_phi_tot, n_theta = n_theta_start, n_r
      vp(1,n_theta_start+1,n_r)            !n_phi = 1, n_theta = n_theta_start+1, n_r
      ...
      vp(n_phi_tot,n_theta_start+1,n_r)
      ...
      vp(1,n_theta_stop,n_r)               !n_phi = 1, n_theta = n_theta_stop, n_r
      vp(2,n_theta_stop,n_r)               !n_phi = 2, n_theta = n_theta_stop, n_r
      ...
      vp(n_phi_tot,n_theta_stop,n_r)       !n_phi = n_phi_tot, n_theta = n_theta_stop, n_r 


      if(l_mag):                           !Only if it is a magnetic case

      !----------------------
      ! Radial magnetic field
      !----------------------

      br(1,n_theta_start,n_r)              !n_phi = 1, n_theta = n_theta_start, n_r
      br(2,n_theta_start,n_r)	           !n_phi = 2, n_theta = n_theta_start, n_r
      ...
      br(n_phi_tot,n_theta_start,n_r)      !n_phi = n_phi_tot, n_theta = n_theta_start, n_r
      br(1,n_theta_start+1,n_r)            !n_phi = 1, n_theta = n_theta_start+1, n_r
      ...
      br(n_phi_tot,n_theta_start+1,n_r)    
      ...
      br(1,n_theta_stop,n_r)               !n_phi = 1, n_theta = n_theta_stop, n_r 
      br(2,n_theta_stop,n_r)               !n_phi = 2, n_theta = n_theta_stop, n_r
      ...
      br(n_phi_tot,n_theta_stop,n_r)       !n_phi = n_phi_tot, n_theta = n_theta_stop, n_r 

      
      !----------------------------------
      ! Theta component of magnetic field
      !----------------------------------

      bt(1,n_theta_start,n_r)              !n_phi = 1, n_theta = n_theta_start, n_r
      bt(2,n_theta_start,n_r)	           !n_phi = 2, n_theta = n_theta_start, n_r
      ...
      bt(n_phi_tot,n_theta_start,n_r)      !n_phi = n_phi_tot, n_theta = n_theta_start, n_r
      bt(1,n_theta_start+1,n_r)            !n_phi = 1, n_theta = n_theta_start+1, n_r
      ...
      bt(n_phi_tot,n_theta_start+1,n_r)
      ...
      bt(1,n_theta_stop,n_r)               !n_phi = 1, n_theta = n_theta_stop, n_r
      bt(2,n_theta_stop,n_r)               !n_phi = 2, n_theta = n_theta_stop, n_r
      ...
      bt(n_phi_tot,n_theta_stop,n_r)       !n_phi = n_phi_tot, n_theta = n_theta_stop, n_r 

     
      !----------------------------------------
      ! Zonal (phi component) of magnetic field
      !----------------------------------------

      bp(1,n_theta_start,n_r)              !n_phi = 1, n_theta = n_theta_start, n_r
      bp(2,n_theta_start,n_r)	           !n_phi = 2, n_theta = n_theta_start, n_r
      ...
      bp(n_phi_tot,n_theta_start,n_r)      !n_phi = n_phi_tot, n_theta = n_theta_start, n_r
      bp(1,n_theta_start+1,n_r)            !n_phi = 1, n_theta = n_theta_start+1, n_r
      ...
      bp(n_phi_tot,n_theta_start+1,n_r)
      ...
      bp(1,n_theta_stop,n_r)               !n_phi = 1, n_theta = n_theta_stop, n_r
      bp(2,n_theta_stop,n_r)               !n_phi = 2, n_theta = n_theta_stop, n_r
      ...
      bp(n_phi_tot,n_theta_stop,n_r)       !n_phi = n_phi_tot, n_theta = n_theta_stop, n_r 

      !-----------------
      !Subsequent blocks
      !-----------------
      
      !Block N+1 in both cases have data at the same radial level but the next theta chunk (n_theta_start + nThetaB, n_theta_stop + n_thetaB)
      
      !After data for all the theta blocks have been written for one radial level, everything above is repeated for the next radial level

The graphic files can be read and visualized using the :py:class:`Surf <magic.Surf>` class as follows:

     >>> S = Surf(tag='TAG')
     >>> S.surf(field = 'vr', r = 0.5, cmap = 'jet', levels = 50)     #Surface map of radial velocity
     >>> S.slice(field = 'br', lon_0 = [0])                           #Longitudinal Slice of radial magnetic field
     >>> S.equat(field = 'entropy')                                   #Equatorial slice of entropy 
