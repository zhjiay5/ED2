!==========================================================================================!
!==========================================================================================!
!    This module will calculate hurricane disturbances rate                   !!!Jiaying   !
! to cohorts. This is done once a day.                                                     !
!------------------------------------------------------------------------------------------!

module hurricane
  implicit none
  contains

  !==========================================================================================!
  !==========================================================================================!
  subroutine hurricane_mortrate(csite,cpatch,ico,ipft,vels)
      use ed_state_vars , only : sitetype                   & ! intent(in)
                               , patchtype                  ! ! Structure
      use hurricane_coms, only : treefallrate_wind_S        & ! intent(in)
                               , treefallrate_wind_L        & ! intent(in)
                               , treefall_wind_threshold    & ! intent(in)
                               , treefall_ndbhclass         & ! intent(in)
                               , treefall_dbh_class         ! ! intent(in)
      use ed_max_dims, only: n_pft

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype),  target     :: csite           ! Current site
      type(patchtype), target     :: cpatch          ! Current patch
      integer        , intent(in) :: ico             ! Current cohort ID
      integer                     :: ipft	     ! pft ID
      real           , intent(in) :: vels            ! velocity of current polygon
      !----- Local variables --------------------------------------------------------------!
      integer                     :: idbhclss        ! DBH class for windthrow mortality
      integer                     :: jpa	     ! patch ID
      integer                     :: jco	     ! cohort ID
      real                        :: n_large_trees   ! 
      real                        :: n_small_trees   ! 
      real                        :: windthrow_mort 

      cpatch%windthrow_mort(ico) = 0.0
      if ( vels >= treefall_wind_threshold ) then
         csite%age(:) = 0.0
	 ! --------------------------------------------------------------- !
	 !  Calculate Structure
	 ! --------------------------------------------------------------- !
          n_large_trees = 0
          n_small_trees = 0
          do jpa=1,csite%npatches
          do jco=1,csite%patch(jpa)%ncohorts
              if ( csite%patch(jpa)%dbh(jco)>=10 )  then
                    n_large_trees = n_large_trees + csite%patch(jpa)%nplant(jco)*csite%area(jpa)
              else if ( csite%patch(jpa)%dbh(jco)<10 .and. csite%patch(jpa)%dbh(jco)>= 2.5 ) then
                    n_small_trees = n_small_trees + csite%patch(jpa)%nplant(jco)*csite%area(jpa)
              end if
          end do
          end do
	 !  End Calculating Structure
	 ! --------------------------------------------------------------- !
	 ! assign mortality based on structure
	 ! --------------------------------------------------------------- !
          do idbhclss=1, treefall_ndbhclass
          if ( cpatch%dbh(ico)>treefall_dbh_class(idbhclss) .and. cpatch%dbh(ico)<=treefall_dbh_class(idbhclss+1) ) then
                if ( n_small_trees > n_large_trees ) then
                    cpatch%windthrow_mort(ico) = treefallrate_wind_S((idbhclss-1)*n_pft+ipft)
                else
                    cpatch%windthrow_mort(ico) = treefallrate_wind_L((idbhclss-1)*n_pft+ipft)
                end if

            write(*,'(3(a,f4.1),(a,i2),(a,f5.3))') 'vels=', vels, '; ws0=',treefall_wind_threshold, &
                         '; dbh=', cpatch%dbh(ico), &
                         '; pft=', ipft, &
                         '; mort=',cpatch%windthrow_mort(ico)

                goto 101
          end if
          end do
	 ! --------------------------------------------------------------- !
      end if

      101 continue

      return
  end subroutine hurricane_mortrate
  !==========================================================================================!
  !==========================================================================================!



  !==========================================================================================!
  !==========================================================================================!
  subroutine hurricane_windthrow(csite,ipa)
      use ed_state_vars, only : patchtype & ! structure
                              , sitetype  ! !
      use pft_coms     , only : c2n_leaf  & ! intent(in)
                              , c2n_stem  & ! intent(in)
                              , l2n_stem  ! ! intent(in)
      use decomp_coms  , only : f_labile  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: ipft
      real                         :: plant_litter
      real                         :: plant_litter_f
      real                         :: plant_litter_s
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !      Add fine root and leaf turnover to the litter.                                !
      !------------------------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         plant_litter   = ( cpatch%balive(ico) + cpatch%bdead(ico) + cpatch%bstorage(ico) )   &
                          * cpatch%nplant(ico) *( 1 - exp(-cpatch%windthrow_mort(ico)) )
         plant_litter_f = plant_litter * f_labile(ipft)
         plant_litter_s = plant_litter - plant_litter_f

         csite%fsc_in(ipa) = csite%fsc_in(ipa) + plant_litter_f
         csite%fsn_in(ipa) = csite%fsn_in(ipa) + plant_litter_f / c2n_leaf(ipft)
         csite%ssc_in(ipa) = csite%ssc_in(ipa) + plant_litter_s
         csite%ssl_in(ipa) = csite%ssl_in(ipa) + plant_litter_s * l2n_stem / c2n_stem(ipft)

         cpatch%nplant(ico) = cpatch%nplant(ico) * exp(-cpatch%windthrow_mort(ico))

         !---- it is a new day, reset windthrow mortality---------------------------------!
         cpatch%windthrow_mort(ico)=0.0
      end do
      return
  end subroutine hurricane_windthrow
  !==========================================================================================!
  !==========================================================================================!



  !==========================================================================================!
  !==========================================================================================!
  subroutine hurricane_defoliate(cpatch,ico,vels, cb_decrement)
      use ed_state_vars, only : patchtype & ! structure
                              , sitetype  ! !
      use hurricane_coms,only : defoliaterate_wind         & ! intent(in)
                              , treefall_wind_threshold    ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target        :: cpatch          ! Current patch
      integer        , intent(in)    :: ico             ! Current cohort ID
      real           , intent(in)    :: vels            ! velocity of current polygon
      real           , intent(inout) :: cb_decrement    ! decrement of carbon balance
      !----- Local variable----------------------------------------------------------------!
      real                           :: bleaf_loss      ! leaf loss due to hurricane
      !------------------------------------------------------------------------------------!

      !--------------------------------------------------! 
      !               Initialize Leafloss                !
      !--------------------------------------------------! 
      bleaf_loss = 0.0

      !--------------------------------------------------! 
      !               Update leafloss                    !
      !--------------------------------------------------! 
      if ( vels >= treefall_wind_threshold ) then
        bleaf_loss = cpatch%bleaf(ico) * defoliaterate_wind 
      end if

      !--------------------------------------------------! 
      !               Update storages                    !
      !--------------------------------------------------! 
      cpatch%bleaf(ico) = cpatch%bleaf(ico) - bleaf_loss
      cpatch%balive(ico) = cpatch%balive(ico) - bleaf_loss

      !--------------------------------------------------! 
      !               Update decrement                   !
      !     will be used to update carbon balance        !
      !--------------------------------------------------! 
      cb_decrement  = cb_decrement + bleaf_loss

      !--------------------------------------------------! 
      !       Accumulate to leaf maintenance             ! 
      !      will be used to calculate litter            !
      !--------------------------------------------------!
      cpatch%leaf_maintenance(ico) = cpatch%leaf_maintenance(ico) + bleaf_loss

      return
  end subroutine hurricane_defoliate
  !==========================================================================================!
  !==========================================================================================!



 ! subroutine hurricane_landslide(cgrid)

 ! end subroutine hurricane_landslide


end module hurricane

