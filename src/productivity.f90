! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! contacts :: David Montenegro Lapola <lapoladm ( at ) gmail.com>
! Author: JP Darela
! This program is based on the work of those that gave us the INPE-CPTEC-PVM2 model

module productivity
  implicit none
  private

  public :: prod


contains

  subroutine prod(dt,light_limit,catm,temp,ts,p0,w,ipar,rh,emax,&
       & cl1_prod,ca1_prod,cf1_prod,beta_leaf,beta_awood,beta_froot,&
       & wmax,psisoil,psi50,klmax,krcmax,psixylem,kxylem,knorm,&
       & ph,ar,nppa,laia,f5,vpd,rm,rg,rc,wue,c_defcit,vm_out,sla,e)

    use types
    use global_par
    use photo_par
    use photo
    use water

!Input
!-----
    real(r_8),dimension(ntraits),intent(in) :: dt     !PLS data
    real(r_4), intent(in) :: temp, ts                 !Mean monthly temperature (oC)
    real(r_4), intent(in) :: p0                       !Mean surface pressure (hPa)
    real(r_8), intent(in) :: w                        !Soil moisture kg m-2
    real(r_4), intent(in) :: ipar                     !Incident photosynthetic active radiation (w/m2)
    real(r_4), intent(in) :: rh,emax                  !Relative humidity/MAXIMUM EVAPOTRANSPIRATION
    real(r_8), dimension(3),intent(in) :: cl1_prod    !Total carbon in each cohort of leaves (kg/m2)
    real(r_8), intent(in) :: catm, cf1_prod, ca1_prod !Carbon in plant tissues (kg/m2)      
    real(r_8), intent(in) :: beta_leaf                !npp allocation to carbon pools (kg/m2/day)
    real(r_8), intent(in) :: beta_awood
    real(r_8), intent(in) :: beta_froot
    real(r_8), intent(in) :: wmax
    real(r_8), intent(in) :: psisoil
    logical(l_1), intent(in) :: light_limit           !True for no ligth limitation

!     Output
!     ------
    real(r_8), intent(out) :: psi50                   !xylem water potential when the plant loses 50% of their maximum xylem conductance (MPa)
    real(r_8), intent(out) :: klmax                   !Maximum xylem conductivity per unit leaf area (kgm-1s-1MPa-1)
    real(r_8), intent(out) :: krcmax                  !Maximum xylem conductance per unit leaf area (molm-2s-1Mpa-1)
    real(r_8), intent(out) :: psixylem                !Xylem water potential (MPa)
    real(r_8), intent(out) :: kxylem                  !Xylem conductance (molm-2s-1MPa-1)
    real(r_8), intent(out) :: knorm                   !Returns normalized xylem conductance (dimensionless)
    real(r_4), intent(out) :: ph                      !Canopy gross photosynthesis (kgC/m2/yr)
    real(r_4), intent(out) :: rc                      !Stomatal resistence (not scaled to canopy!) (s/m)
    real(r_8), intent(out) :: laia                    !Autotrophic respiration (kgC/m2/yr)
    real(r_4), intent(out) :: ar                      !Leaf area index (m2 leaf/m2 area)
    real(r_4), intent(out) :: nppa                    !Net primary productivity (kgC/m2/yr)
    real(r_4), intent(out) :: vpd
    real(r_8), intent(out) :: f5                      !Water stress response modifier (unitless)
    real(r_4), intent(out) :: rm                      !autothrophic respiration (kgC/m2/day)
    real(r_4), intent(out) :: rg
    real(r_4), intent(out) :: wue
    real(r_4), intent(out) :: c_defcit                !Carbon deficit gm-2 if it is positive, aresp was greater than npp + sto2(1)
    real(r_8), intent(out) :: sla
    real(r_8), intent(out) :: e                       !transpiration (molm2s)
    real(r_8), intent(out) :: vm_out


!     Internal
!     --------

    real(r_8) :: tleaf    !leaf/wood turnover time (yr) 
    real(r_8) :: awood            
    real(r_8) :: g1
    real(r_8) :: c4
    real(r_8) :: wd
    real(r_8) :: sla_var

    real(r_8) :: n2cl
    real(r_8) :: n2cl_resp
    real(r_8) :: n2cw_resp
    real(r_8) :: n2cf_resp
    real(r_8) :: p2cl
    integer(i_4) :: c4_int
    real(r_8) :: jl_out

    real(r_8), dimension(3) :: f1      !Leaf level gross photosynthesis (molCO2/m2/s)
    real(r_8) :: f1a                   !auxiliar_f1
    real(r_8), dimension(3) :: umol_penalties = (/-0.4, 1.0, 0.6/) !Penalization in photosynthesis for each cohort, defined by Wu et al (2016) and Albert et al (2018)
    real(r_8), dimension(3) :: leaf_age
    real(r_8), dimension(3) :: penalization_by_age
    real(r_8) :: age_crit
    real(r_8) :: cl_total              !Carbon sum of all the cohots (kg/m2)
    real(r_4) :: rc_pot, rc_aux, e_pot
    integer(i_4) :: i

    !Hydraulic parameters
    real(r_8) :: diameter
    real(r_8) :: height1
    !real(r_8) :: psi50
    !real(r_8) :: klmax
    !real(r_8) :: krcmax
    !real(r_8) :: psixylem
    !real(r_8) :: kxylem
    !real(r_8) :: knorm

! Getting pls parameters

    g1  = dt(1)
    tleaf = dt(3)
    awood = dt(7)
    c4  = dt(9)
    n2cl = dt(10)
    n2cl_resp = n2cl
    n2cw_resp = dt(11)
    n2cf_resp = dt(12)
    p2cl = dt(13)
    wd = dt(18)
    sla_var = dt(19)

    !Simulation of leaf demography
    !Obtain critical age
    age_crit = (tleaf / 3.0) * 2.0

    !Obtain leaf age (a) - middle age of each cohort
    leaf_age(1) = (tleaf * (1.0/12.0))
    leaf_age(2) = (tleaf * (1.0/2.0))
    leaf_age(3) = (tleaf * (5.0/6.0))

    do i = 1, 3
        penalization_by_age(i) = leaf_age_factor(umol_penalties(i), age_crit, leaf_age(i))
    enddo

    !print*,'fa jovem',penalization_by_age(1)
    !print*,'fa madura',penalization_by_age(2)
    !print*,'fa velha',penalization_by_age(3)

    !Obtain total carbon of the leaf cohorts
    cl_total = sum(cl1_prod)

    n2cl = real(n2cl * (cl_total * 1e3), r_4) ! N in leaf g m-2
    p2cl = real(p2cl * (cl_total * 1e3), r_4) ! P in leaf g m-2

    c4_int = idnint(c4)


!     ==============
!     Photosynthesis
!     ==============
!   rate (molCO2/m2/s)

    call photosynthesis_rate(catm,temp,p0,ipar,light_limit,c4_int,tleaf,n2cl,&
         & p2cl,cl1_prod(:),f1a,vm_out,jl_out)

    ! VPD
    !========
    vpd = vapor_p_defcit(temp,rh)

    !Stomatal resistence
    !===================
    rc_pot = canopy_resistence_pot(vpd, f1a, g1, catm) ! Potential RCM leaf level - s m-1

    e_pot = transpiration(rc_pot, p0, vpd, 1)
    !print*, 'e_pot:',e_pot

    !==========================
    !  Hydraulic without stress
    !==========================

    diameter = diameter_pls(wd,ca1_prod)
    height1 = height_pls(diameter)
    !print*,'height:',height1,'cawood:',ca1_prod

    !   P50
    !=========
    psi50 = psi_fifty(wd,ca1_prod)
    !print*,'P50',psi50, 'wd',wd

    ! Klmax
    !=========
    klmax = conductivity_xylemleaf(wd,jl_out* 1e6,ca1_prod)
    !print*,'klmax',klmax

    ! Krcmax
    !=========
    krcmax = conductance_xylemax(klmax,height1,ca1_prod)   
    !print*,'krcmax',krcmax,'klmax',klmax

    ! Psixylem
    !=========
    psixylem = xylem_waterpotential(psisoil,height1,ca1_prod)
    !print*,'psixylem',psixylem, 'psisoil',psisoil

    ! k xylem
    !=========
    kxylem = xylem_conductance(krcmax,psixylem,psi50,ca1_prod)
    !print*,'KXYLEM:',kxylem, 'P50:',psi50, 'KRCMAX:',krcmax

    ! k xylem
    !=========
    knorm = conductance_normalized(krcmax,kxylem,ca1_prod)
    !print*,'knorm:',knorm,'krcmax:',krcmax,'kxylem:',kxylem

    !Water stress response modifier (dimensionless)
    !----------------------------------------------
    f5 =  water_stress_modifier(w, cf1_prod, rc_pot, emax, wmax, knorm, ca1_prod)
    !print*,'f5',f5


!     Photosysthesis minimum and maximum temperature
!     ----------------------------------------------
    if ((temp.ge.-10.0).and.(temp.le.50.0)) then
        do i = 1,3
            f1(i) = f1a * f5 * penalization_by_age(i) ! water stress factor and factor age ! Ancient floating-point underflow spring (from CPTEC-PVM2)
        enddo
    else
        f1 = 0.0      !Temperature above/below photosynthesis windown
    endif

    rc_aux = canopy_resistence_real(vpd, f1(:), g1, catm)  ! RCM leaf level -!s m-1

    wue = water_ue(f1(:), rc_aux, p0, vpd)

    ! Calcula a transpiração em mm/s
    e = transpiration(rc_aux, p0, vpd, 2)
    !print*, 'e_real:',e

    ! Leaf area index (m2/m2)
    ! recalcula rc e escalona para dossel
    ! laia = 0.2D0 * dexp((2.5D0 * f1)/p25)
    sla = spec_leaf_area(tleaf)  ! m2 g-1  ! Convertions made in leaf_area_index &  gross_ph + calls therein

    laia = leaf_area_index(cl1_prod(:), sla_var)
    !print*,'lai:',laia,'cl1_prod(1):',cl1_prod(1),'cl1_prod(2):',cl1_prod(2),'cl1_prod(3):',cl1_prod(3),'sla:',sla_var

    rc = rc_aux !* real(laia,kind=r_4) ! RCM -!s m-1 ! CANOPY SCALING --

!     Canopy gross photosynthesis (kgC/m2/yr)
!     =======================================x

    ph = gross_ph(f1(:),cl1_prod(:), sla_var)        ! kg m-2 year-1

!     Autothrophic respiration
!     ========================
!     Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
    rm = m_resp(temp,ts,cl1_prod(:),cf1_prod,ca1_prod &
         &,n2cl_resp,n2cw_resp,n2cf_resp,awood)

! c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
! c     2003; Levis et al. 2004)
    rg = g_resp(beta_leaf,beta_awood, beta_froot,awood)

    if (rg.lt.0) then
       rg = 0.0
    endif

!     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
!     Respiration minimum and maximum temperature
!     -------------------------------------------
    if ((temp.ge.-10.0).and.(temp.le.50.0)) then
       ar = rm + rg
    else
       ar = 0.0               !Temperature above/below respiration windown
    endif
!     Net primary productivity(kgC/m2/yr)
!     ====================================
    nppa = ph - ar
! this operation affects the model mass balance
! If ar is bigger than ph, what is the source or respired C?

    if(ar .gt. ph) then
       c_defcit = ((ar - ph) * 2.73791) ! tranform kg m-2 year-1 in  g m-2 day-1
       nppa = 0.0
    else
       c_defcit = 0.0
    endif

  end subroutine prod

end module productivity
