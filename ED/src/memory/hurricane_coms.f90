
module hurricane_coms !!!Jiaying
   use ed_max_dims, only : ndbhmax          &
                         , n_pft

   implicit none

   integer                            :: hurricane_ndbhclass      
   real                               :: hurricane_wind_threshold
   real                               :: hurricane_distrate_a
   real                               :: hurricane_distrate_b
   real , dimension(ndbhmax*n_pft)    :: hurricane_s_a
   real , dimension(ndbhmax*n_pft)    :: hurricane_s_b
   real , dimension(ndbhmax*n_pft)    :: hurricane_s
   real , dimension(ndbhmax+1)        :: hurricane_dbh_class
   real                               :: hurricane_defoliaterate


end module hurricane_coms
