module rieutord
   
   !Module containing variables for inertial mode forcing like Rieutord.
   !Forcing is implemented in updateZ.
   !
   !A toroidal boundary force of symmetry (m,m) is used.
   !Can be applied to inner or outer boundary.

   use precision_mod
  
   logical :: l_Ri                  !Decide whether to use forcing at all
   logical :: l_RiIc,l_RiMa         !Switches to decide which boundary
                                    !should be forced           
   integer :: m_RiIc,m_RiMa         !Order of forcing at boundaries
   real(cp) :: amp_RiIc,omega_RiIc  !Inner boundary
   real(cp) :: amp_RiMa,omega_RiMa  !Outer boundary

end module rieutord
