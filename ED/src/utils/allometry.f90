!==========================================================================================!
!==========================================================================================!
! MODULE: ALLOMETRY  
!> \brief This module is a library with several allometric relationships.
!------------------------------------------------------------------------------------------!
module allometry

implicit none

contains
   !=======================================================================================!
   !=======================================================================================!
   real function h2dbh(h,ipft)

      use pft_coms    , only : is_tropical & ! intent(in)
         , b1Ht        & ! intent(in), lookup table
         , b2Ht        & ! intent(in), lookup table
         , hgt_ref     ! ! intent(in)
      use ed_misc_coms, only : iallom      ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: h
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Select which type of model we are running.                                    !
      !------------------------------------------------------------------------------------!
      !----- Size- and age-structure (typical ED model). ----------------------------------!
      if (is_tropical(ipft)) then
         select case (iallom)
            case (0,1)
               !----- Default ED-2.1 allometry. -------------------------------------------!
               h2dbh = exp((log(h)-b1Ht(ipft))/b2Ht(ipft))
            case default
               !----- Poorter et al. (2006) allometry. ------------------------------------!
               h2dbh =  ( log(hgt_ref(ipft) / ( hgt_ref(ipft) - h ) ) / b1Ht(ipft) )       &
                  ** ( 1.0 / b2Ht(ipft) )
         end select
      else ! Temperate
         h2dbh = log(1.0-(h-hgt_ref(ipft))/b1Ht(ipft))/b2Ht(ipft)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function h2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function dbh2h(ipft, dbh)
      use pft_coms    , only : is_tropical & ! intent(in)
         , dbh_crit    & ! intent(in)
         , b1Ht        & ! intent(in)
         , b2Ht        & ! intent(in)
         , hgt_ref     & ! intent(in)
         , hgt_max     ! ! intent(in)
      use ed_misc_coms, only : iallom      & ! intent(in)
         , ibigleaf    ! ! intent(in)
      use ed_state_vars, only:  patchtype  ! ! structure

      !----- Arguments --------------------------------------------------------------------!
      integer       , intent(in) :: ipft
      real          , intent(in) :: dbh
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Select which type of model we are running.                                    !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
         case (0)
            !----- Size- and age-structure (typical ED model). ----------------------------!
            if (is_tropical(ipft)) then
               mdbh = min(dbh,dbh_crit(ipft))
               select case (iallom)
                  case (0,1)
                     !----- Default ED-2.1 allometry. -------------------------------------!
                     dbh2h = exp (b1Ht(ipft) + b2Ht(ipft) * log(mdbh) )
                  case default
                     !----- Poorter et al. (2006) allometry. ------------------------------!
                     dbh2h = hgt_ref(ipft) * (1. - exp(-b1Ht(ipft) * mdbh ** b2Ht(ipft)))
               end select
            else !----- Temperate PFT allometry. ------------------------------------------!
               dbh2h = hgt_ref(ipft) + b1Ht(ipft) * (1.0 - exp(b2Ht(ipft) * dbh))
            end if

         case (1)
            !------------------------------------------------------------------------------!
            !     Big-leaf version of ED. DBH is not really meaningful, but in the big-leaf!
            ! model the typical allometry doesn't really make sense so we impose maximum   !
            ! height.                                                                      !
            !------------------------------------------------------------------------------!
            dbh2h = hgt_max(ipft)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end function dbh2h
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   real function dbh2bd(dbh,ipft)

      use pft_coms    , only : C2B         & ! intent(in)
         , b1Bs_small  & ! intent(in), lookup table
         , b2Bs_small  & ! intent(in), lookup table
         , b1Bs_large  & ! intent(in), lookup table
         , b2Bs_large  & ! intent(in), lookup table
         , is_grass    & ! intent(in)
         , dbh_crit    ! ! intent(in), lookup table
      use ed_misc_coms, only : igrass      ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!

         if (is_grass(ipft) .and. igrass==1) then
            dbh2bd = 0.0
         else if (dbh <= dbh_crit(ipft)) then
            dbh2bd = b1Bs_small(ipft) / C2B * dbh ** b2Bs_small(ipft)
         else
            dbh2bd = b1Bs_large(ipft) / C2B * dbh ** b2Bs_large(ipft)
      end if
      return
   end function dbh2bd
   !=======================================================================================!
   !=======================================================================================!

   !=======================================================================================!
   !=======================================================================================!
   ! FUNCTION DBH2SF             
   !< \brief Calculate sapwood fraction from DBH
   !< \author Xiangtao Xu, 31 Jan. 2018
   !---------------------------------------------------------------------------------------!
   real function dbh2sf(dbh,ipft)

      use consts_coms   , only : pio4               ! ! intent(in)
      use pft_coms      , only : is_grass           & ! intent(in)
                               , b1SA               & ! intent(in)
                               , b2SA               ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!

      if (is_grass(ipft)) then
          ! All of the 'stem' of grass is sapwood
          ! But we actually not use the sapwood area for grass
          ! This is just to keep consistency across species
          dbh2sf = 1.
      else
          ! For trees
          dbh2sf = min(1.0,                             &
                       (b1SA(ipft) * dbh ** b2SA(ipft)) & ! sapwood area cm2
                    /  (pio4 * dbh ** 2 )               & ! total basal area
                      )
      endif

      return 
   end function dbh2sf
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the diameter at the breast height given the structural     !
   ! (dead) biomass and the PFT type.  In addition to these inputs, it is also required to !
   ! give a first guess for DBH, which can be the previous value, for example.   This is   !
   ! only used for the new allometry, because we must solve an iterative root-finding      !
   ! method.                                                                               !
   !---------------------------------------------------------------------------------------!
   real function bd2dbh(ipft, bdead)
      use pft_coms    , only : b1Bs_small  & ! intent(in), lookup table
         , b2Bs_small  & ! intent(in), lookup table
         , b1Bs_large  & ! intent(in), lookup table
         , b2Bs_large  & ! intent(in), lookup table
         , bdead_crit  & ! intent(in), lookup table
         , C2B         ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in) :: ipft      ! PFT type                            [         ---]
      real   , intent(in) :: bdead     ! Structural (dead) biomass           [   kgC/plant]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Decide which coefficients to use based on the critical bdead.                   !
      !------------------------------------------------------------------------------------!
         if (bdead <= bdead_crit(ipft)) then
            bd2dbh = (bdead / b1Bs_small(ipft) * C2B)**(1.0/b2Bs_small(ipft))
         else
            bd2dbh = (bdead / b1Bs_large(ipft) * C2B)**(1.0/b2Bs_large(ipft))
      end if
      !------------------------------------------------------------------------------------!

      return
   end function bd2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function computes the maximum leaf biomass [kgC/m2] given the size (either  !
   ! height or DBH).  Trees would always use DBH as the starting point, but grasses may    !
   ! use DBH (old style) or height (new style).  This replaces dbh2bl and h2bl with a      !
   ! single generic function that should be used by all plants.                            !
   !---------------------------------------------------------------------------------------!
   real function size2bl(dbh,hite,ipft)
      use pft_coms    , only : dbh_adult   & ! intent(in), lookup table
         , dbh_crit    & ! intent(in), lookup table
         , C2B         & ! intent(in), lookup table
         , b1Bl_small  & ! intent(in), lookup table
         , b2Bl_small  & ! intent(in), lookup table
         , b1Bl_large  & ! intent(in), lookup table
         , b2Bl_large  & ! intent(in), lookup table
         , is_liana    & ! intent(in)
         , is_grass    ! ! intent(in)
      use ed_misc_coms, only : igrass      ! ! intent(in)
      use ed_state_vars, only: patchtype  ! ! structure

      !----- Arguments --------------------------------------------------------------------!
      real          , intent(in) :: dbh
      real          , intent(in) :: hite
      integer       , intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      real                :: ldbh
      real                :: liana_dbh_crit = 26.0 !< liana specific critical dbh
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! lianas dbh_crit would be too small for the leaf biomass. We choose ldbh to be the  !
      ! DBH that give lianas the same (~) Bl_max as late successionals                     !
      !------------------------------------------------------------------------------------!
      if (is_liana(ipft)) then
         ldbh = liana_dbh_crit
      else
         ldbh = dbh_crit(ipft)
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Get the DBH, or potential DBH in case of grasses.                              !
      !------------------------------------------------------------------------------------!
      if (is_grass(ipft) .and. igrass == 1) then 
         !---- Use height for new grasses. ------------------------------------------------!
         mdbh   = min(h2dbh(hite,ipft),ldbh)
      else
         !---- Use dbh for trees. ---------------------------------------------------------!
         mdbh   = min(dbh,ldbh)
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find leaf biomass depending on the tree size.                                  !
      !------------------------------------------------------------------------------------!
         if (mdbh < dbh_adult(ipft)) then
            size2bl = b1Bl_small(ipft) / C2B * mdbh ** b2Bl_small(ipft)
         else
            size2bl = b1Bl_large(ipft) / C2B * mdbh ** b2Bl_large(ipft)
      end if
      !------------------------------------------------------------------------------------!


      return
   end function size2bl
      !====================================================================================!
      !====================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !      This function determines an effective DBH for grasses given their leaf biomass.  !
   ! DBH has no real meaning for grasses with the new allometry.                           !
   !---------------------------------------------------------------------------------------!
   real function bl2dbh(bleaf,ipft)
      use pft_coms,     only : dbh_crit    & ! intent(in), lookup table
         , hgt_max     & ! intent(in), lookup table
         , is_grass    & ! intent(in)
         , C2B         & ! intent(in)
         , b1Bl_small  & ! intent(in), lookup table
         , b2Bl_small  & ! intent(in), lookup table
         , b1Bl_large  & ! intent(in), lookup table
         , b2Bl_large  & ! intent(in), lookup table
         , bleaf_adult ! ! intent(in), lookup table
      use ed_misc_coms, only : igrass      ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)      :: bleaf
      integer, intent(in)      :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      !------------------------------------------------------------------------------------!



         !----- Find out whether this is an adult tree or a sapling/grass. ----------------!
         if (bleaf < bleaf_adult(ipft)) then
            mdbh = (bleaf * C2B / b1Bl_small(ipft) ) ** (1./b2Bl_small(ipft))
         else
            mdbh = (bleaf * C2B / b1Bl_large(ipft) ) ** (1./b2Bl_large(ipft))
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     For grasses, limit maximum effective dbh by maximum height.                 !
         !---------------------------------------------------------------------------------!
         if (is_grass(ipft) .and. igrass == 1) then
            bl2dbh = min(mdbh, h2dbh(hgt_max(ipft),ipft))
         else
            bl2dbh = min(mdbh, dbh_crit(ipft))
      end if
      !------------------------------------------------------------------------------------!


      return
   end function bl2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function determines the height for grasses given their leaf biomass.        !
   !---------------------------------------------------------------------------------------!
   real function bl2h(bleaf,ipft)
      use pft_coms,      only:  hgt_max    ! ! intent(in), lookup table
      
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)      :: bleaf
      integer, intent(in)      :: ipft
      !------------------------------------------------------------------------------------!

      !----- Use existing allometric equations to convert leaves to height. ---------------!
      bl2h = min(hgt_max(ipft),dbh2h(ipft,bl2dbh(bleaf,ipft)))

      return
   end function bl2h
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !    Canopy Area allometry from Dietze and Clark (2008).                                !
   !---------------------------------------------------------------------------------------!
   real function dbh2ca(dbh,hite,sla,ipft)
      use ed_misc_coms, only : iallom      ! ! intent(in)
      use pft_coms    , only : dbh_crit    & ! intent(in)
         , hgt_max     & ! intent(in)
         , is_grass    & ! intent(in)
         , is_liana    & ! intent(in)
         , b1Ca        & ! intent(in)
         , b2Ca        ! ! intent(in)
      use ed_misc_coms, only : igrass      ! ! intent(in)
      use ed_state_vars, only: patchtype  ! ! structure

      !----- Arguments --------------------------------------------------------------------!
      real          , intent(in) :: dbh
      real          , intent(in) :: hite
      real          , intent(in) :: sla
      integer       , intent(in) :: ipft
      !----- Internal variables -----------------------------------------------------------!
      real                :: loclai ! The maximum local LAI for a given DBH
      real                :: ldbh
      real                :: liana_dbh_crit = 26.0 !< liana specific critical dbh
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! lianas dbh_crit would be too small for the leaf biomass. We choose ldbh to be the  !
      ! DBH that give lianas the same (~) Bl_max as late successionals                     !
      !------------------------------------------------------------------------------------!
      if (is_liana(ipft)) then
         ldbh = liana_dbh_crit
      else
         ldbh = dbh_crit(ipft)
      end if
      !------------------------------------------------------------------------------------!

      if (dbh < tiny(1.0)) then
         loclai = 0.0
         dbh2ca = 0.0
      else

         !----- make this function generic to size, not just dbh. -------------------------!

         loclai = sla * size2bl(dbh,hite,ipft)

         select case (iallom)
            case (0)
               !----- No upper bound in the allometry. ------------------------------------!
               dbh2ca = b1Ca(ipft) * dbh ** b2Ca(ipft)

            case default
               !----- Impose a maximum crown area. ----------------------------------------!
               if (is_grass(ipft) .and. igrass==1) then
                  dbh2ca = b1Ca(ipft) * min(dbh,h2dbh(hgt_max(ipft),ipft) ) ** b2Ca(ipft)
               else
                  dbh2ca = b1Ca(ipft) * min(dbh, ldbh                     ) ** b2Ca(ipft)
               end if
         end select
      end if

      !----- Local LAI / Crown area should never be less than one. ------------------------!
      dbh2ca = min (loclai, dbh2ca)

      return
   end function dbh2ca
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the standing volume of a tree.                              !
   !---------------------------------------------------------------------------------------!
   real function dbh2vol(hgt,dbh,ipft)
      use pft_coms    , only : b1Vol    & ! intent(in)
         , b2Vol    ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: hgt
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!


      dbh2vol = b1Vol(ipft) * hgt * dbh ** b2Vol(ipft)

      return
   end function dbh2vol
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   integer function dbh2krdepth(hite,dbh,ipft,lsl)
      use ed_misc_coms, only : iallom   ! ! intent(in)
      use grid_coms   , only : nzg      ! ! intent(in)
      use soil_coms   , only : slz      ! ! intent(in)
      use pft_coms    , only : b1Rd     & ! intent(in)
         , b2Rd     ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: hite
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      integer, intent(in) :: lsl
      !----- Local variables --------------------------------------------------------------!
      real                :: volume
      real                :: root_depth
      integer             :: k
      !------------------------------------------------------------------------------------!


      !----- Grasses get a fixed rooting depth of 70 cm. ----------------------------------!
      select case (iallom)
         case (0)
            !------------------------------------------------------------------------------!
            !    Original ED-2.1 (I don't know the source for this equation, though).      !
            !------------------------------------------------------------------------------!
            volume     = dbh2vol(hite,dbh,ipft)
            root_depth = b1Rd(ipft)  * volume ** b2Rd(ipft)

         case default
            !------------------------------------------------------------------------------!
            !    This is just a test allometry, that imposes root depth to be 0.5 m for    !
            ! plants that are 0.15-m tall, and 5.0 m for plants that are 35-m tall.        !
            !------------------------------------------------------------------------------!
            root_depth = b1Rd(ipft) * hite ** b2Rd(ipft)
      end select


      !------------------------------------------------------------------------------------!
      !     Root depth is the maximum root depth if the soil is that deep.  Find what is   !
      ! the deepest soil layer this root can go.                                           !
      !------------------------------------------------------------------------------------!
      dbh2krdepth = nzg
      do k=nzg,lsl+1,-1
         if (root_depth < slz(k)) dbh2krdepth = k-1
      end do
      !------------------------------------------------------------------------------------!

      return
   end function dbh2krdepth
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function finds the trunk height.  Currently this is based on the following    !
   ! reference, which is for a site in Bolivia:                                            !
   !                                                                                       !
   ! Poorter L., L. Bongers, F. Bongers, 2006: Architecture of 54 moist-forest tree        !
   !     species: traits, trade-offs, and functional groups. Ecology, 87, 1289-1301.       !
   !                                                                                       !
   ! Crown length is defined as the height of the tree minus the height of the first leaf  !
   !---------------------------------------------------------------------------------------!
   real function h2crownbh(height,ipft)
      use pft_coms, only : b1Cl & ! intent(in)
         , b2Cl ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: height
      integer, intent(in) :: ipft
      !----- Local variables. -------------------------------------------------------------!
      real                :: crown_length
      !------------------------------------------------------------------------------------!
      
      crown_length = b1Cl(ipft) * height ** b2Cl(ipft)
      h2crownbh    = max(0.05,height - crown_length)

      return
   end function h2crownbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the total above ground biomass (wood + leaves)              !
   !---------------------------------------------------------------------------------------!
   real function ed_biomass(cpatch,ico)
      use pft_coms,      only : agf_bs  !! intent(in)
      use ed_state_vars, only : patchtype     ! ! Structure

      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target :: cpatch
      integer, intent(in)     :: ico
      !------------------------------------------------------------------------------------!

      ed_biomass = cpatch%bleaf(ico) + cpatch%bsapwooda(ico)                               &
                   + (cpatch%bdead(ico) * agf_bs(cpatch%pft(ico)))

      return
   end function ed_biomass
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine estimates the tree area indices, namely leaf, branch(plus twigs), !
   ! and stem.  For the leaf area index (LAI), we use the specific leaf area (SLA), a      !
   ! constant.  The wood area index WAI is found using the model proposed by J�rvel�      !
   ! (2004) to find the specific projected area.                                           !
   !                                                                                       !
   ! G. Hormann, S. Irrgan, H. Jochheim, M. Lukes, H. Meesenburg, J. Muller, B. Scheler,   !
   !    J. Scherzer, G. Schuler, B. Schultze, B. Strohbach, F. Suckow, M. Wegehenkel, and  !
   !    G. Wessolek.   Wasserhaushalt von waldokosystemen: methodenleitfaden zur           !
   !    bestimmung der wasserhaushaltskomponenten auf level II-flachen. Technical note,    !
   !    Bundesministerium fur Verbraucherschutz, Ernahrung und Landwirtschaft (BMVEL),     !
   !    Bonn, Germany, 2003. URL http://www.wasklim.de/download/Methodenband.pdf.          !
   !                                                                                       !
   ! P. C. Olivas, S. F. Oberbauer, D. B. Clark, D. A. Clark, M. G. Ryan, J. J. O'Brien,   !
   !    and H. Ordonez.  Comparison of direct and indirect methods for assessing leaf area !
   !    index across a tropical rain forest landscape.  Agric. For. Meteorol., 177,        !
   !    110--116, 2013.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine area_indices(cpatch, ico)
      use ed_state_vars, only : patchtype     ! ! Structure
      use pft_coms     , only : dbh_crit        & ! intent(in)
                              , b1WAI           & ! intent(in)
                              , b2WAI           ! ! intent(in)
      use rk4_coms     , only : ibranch_thermo  ! ! intent(in)
      use ed_misc_coms , only : iallom          ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target :: cpatch
      integer, intent(in)     :: ico
      !------------------------------------------------------------------------------------!

      !----- Local variables --------------------------------------------------------------!
      real                    :: bleaf
      real                    :: nplant
      real                    :: hite
      real                    :: dbh
      integer                 :: ipft
      !------------------------------------------------------------------------------------!
      
      bleaf  = cpatch%bleaf(ico)
      nplant = cpatch%nplant(ico)
      hite   = cpatch%hite(ico)
      ipft   = cpatch%pft(ico)
      dbh    = cpatch%dbh(ico)

      !----- First, we compute the LAI ----------------------------------------------------!
      cpatch%lai(ico) = bleaf * nplant * cpatch%sla(ico)

      !----- Find the crown area. ---------------------------------------------------------!
      cpatch%crown_area(ico) = min(1.0, nplant * dbh2ca(dbh,hite,cpatch%sla(ico),ipft))

      !------------------------------------------------------------------------------------!
      !     Here we check whether we need to compute the branch, stem, and effective       !
      ! branch area indices.  These are only needed when branch thermodynamics is used,    !
      ! otherwise, simply assign zeroes to them.                                           !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (0) 
         !----- Ignore branches and trunk. ------------------------------------------------!
         cpatch%wai(ico)  = 0.
         !---------------------------------------------------------------------------------!

      case (1,2)
         !---------------------------------------------------------------------------------!
         !     Decide the WAI according to the allometry.                                  !
         !---------------------------------------------------------------------------------!
         select case (iallom)
         case (3)
            !------------------------------------------------------------------------------!
            !     Assume a simple extrapolation based on Olivas et al. (2013).  WAI is     !
            ! always 11% of the potential LAI.                                             !
            !------------------------------------------------------------------------------!
            cpatch%wai(ico) = 0.11 * nplant * cpatch%sla(ico) * size2bl(dbh,hite,ipft)
            !------------------------------------------------------------------------------!
         case default
            !----- Solve branches using the equations from Hormann et al. (2003) ----------!
            cpatch%wai(ico) = nplant * b1WAI(ipft) * min(dbh,dbh_crit(ipft)) ** b2WAI(ipft)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine area_indices
   !=======================================================================================!
   !=======================================================================================!
end module allometry
!==========================================================================================!
!==========================================================================================!


