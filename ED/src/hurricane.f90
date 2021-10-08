!==========================================================================================!
!    This subroutine will calculate hurricane disturbances rate               !!!Jiaying   !
!      accumulate daily disturbance rate to monthly and patch level to polygon level       !
!==========================================================================================!
subroutine hurricane_distrate(cgrid)
   use ed_state_vars , only : edtype                 & ! structure
                            , polygontype            & ! structure
                            , sitetype               & ! structure
                            , patchtype              ! ! structure
   use ed_misc_coms  , only : current_time           ! ! intent(in)
   use hurricane_coms, only : hurricane_distrate_a   & ! intent(in)
                            , hurricane_distrate_b   & ! intent(in)
                            , hurricane_s_a          & ! intent(in)
                            , hurricane_s_b          & ! intent(in)
                            , hurricane_s            & ! intent(inout)
                            , hurricane_wind_threshold  ! ! intent(in)
   use consts_coms   , only : lnexp_min              & ! intent(in)
                            , lnexp_max              ! ! intent(in)

   implicit none
   !----- Arguments --------------------------------------------------------------------!
   type(edtype)      , target     :: cgrid
   type(polygontype) , pointer    :: cpoly
   type(sitetype)    , pointer    :: csite
   type(patchtype)   , pointer    :: cpatch
   integer                        :: ipy
   integer                        :: isi
   integer                        :: ipa
   integer                        :: ico
   integer                        :: imon
   real                           :: n_large_trees     ! # of large trees in each patch 
   real                           :: n_small_trees     ! # of small trees in each patch
   real                           :: prop_nl           ! proportion of large trees in each patch 
   real                           :: distrate          ! disturbance rate 0-1
   real                           :: hurricane_m_patch ! area mortality (%) of each site
						       ! weighted by the area of each patch 

   !----- Current month. ------------------------------------------------------------------!
   imon = current_time%month
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Loop over polygons                                                               !
   !---------------------------------------------------------------------------------------!
   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      !---- initialize the hurricane disturbance and survival rates after each new_year  --!
      if ( current_time%month == 1 .and. current_time%date == 2 ) then 
          hurricane_s = 1.0
          cpoly%lambda_hurricane(:,:) = 0.0
      end if 

      !------------------------------------------------------------------------------------!
      !     Loop over all sites.                                                           !
      !------------------------------------------------------------------------------------!
      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         if ( cpoly%met(isi)%vels >= hurricane_wind_threshold ) then
        !----------------------------------------------------------------------------------!
        !     Loop over patches.                                                           !
        !----------------------------------------------------------------------------------!
         n_large_trees = 0
         n_small_trees = 0
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
	    ! --------------------------------------------------------------- !
	    !  Calculate Structure
	    ! --------------------------------------------------------------- !
            do ico=1,cpatch%ncohorts
               if ( cpatch%dbh(ico)>=10 )  then
                   n_large_trees = n_large_trees + cpatch%nplant(ico) * csite%area(ipa)
               else if ( cpatch%dbh(ico)<10 .and. cpatch%dbh(ico)>= 2.5 ) then
                   n_small_trees = n_small_trees + cpatch%nplant(ico) * csite%area(ipa)
               end if
            end do
         end do patchloop
	 ! --------------------------------------------------------------- !
         ! assign mortality based on structure
	 ! --------------------------------------------------------------- !
         prop_nl = n_large_trees / (n_large_trees + n_small_trees)
         distrate = 1. / ( 1. + exp( - hurricane_distrate_a * (prop_nl - hurricane_distrate_b)) ) 
         cpoly%lambda_hurricane(imon,isi) = cpoly%lambda_hurricane(imon,isi)            & 
                                             + min( lnexp_max, -log(1-distrate) )
         hurricane_s = hurricane_s                                                      &
                         * 1. / ( 1. + exp( - hurricane_s_a * (prop_nl - hurricane_s_b)) )
         write(*,'(a,f5.4)') 'prop_nl', prop_nl
         !---------------------------------------------------------------------------------!
         end if
      end do siteloop
      !------------------------------------------------------------------------------------!
   end do polyloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine hurricane_distrate
!==========================================================================================!
!==========================================================================================!


