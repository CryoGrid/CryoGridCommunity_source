function [psim, psic] = GetPsiRSL(vegetation, za, hc, dt, beta,zeta,obu, PrSc) 
% za, hc, dt, beta, zeta, obu, PrSc
%(za, hc, dt, beta, zeta, obu, PrSc, psim, psic)


%   !-----------------------------------------------------------------------
%   subroutine GetPsiRSL (za, hc, dt, beta, zeta, obu, PrSc, psim, psic)
%     !
%     ! !DESCRIPTION:
%     ! Calculate the stability functions psi for momentum and scalars. The
%     ! returned function values (psim, psic) are the Monin-Obukhov psi functions
%     ! and additionally include the roughness sublayer psihat functions.
%     ! These are evaluated between the height za and at the canopy height hc.
%     !
%     ! !USES:
%     use clm_varctl, only : turb_type
%     use clm_varcon, only : vkc
%     use clm_varcon, only : dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
%     !
%     ! !ARGUMENTS:
%     implicit none
%     real(r8), intent(in)  :: za        ! Atmospheric height (m)
%     real(r8), intent(in)  :: hc        ! Canopy height (m)
%     real(r8), intent(in)  :: dt        ! Height below canopy top (dt = ztop - zdisp)
%     real(r8), intent(in)  :: beta      ! Value of u*/u(h) at canopy top
%     real(r8), intent(in)  :: zeta      ! Monin-Obukhov stability parameter (za-d)/L
%     real(r8), intent(in)  :: obu       ! Obukhov length (m)
%     real(r8), intent(in)  :: PrSc      ! Turbulent Prandtl (Schmidt) number at canopy top
%     real(r8), intent(out) :: psim      ! psi function for momentum including RSL influence
%     real(r8), intent(out) :: psic      ! psi function for scalars  including RSL influence
%     !
%     ! !LOCAL VARIABLES:
%     real(r8) :: phim                   ! Monin-Obukhov phi function for momentum at canopy top
%     real(r8) :: phic                   ! Monin-Obukhov phi function for scalars at canopy top
%     real(r8) :: c1                     ! RSL magnitude multiplier
%     real(r8) :: psihat1                ! Multiply by c1 to get the RSL psihat function evaluated at za
%     real(r8) :: psihat2                ! Multiply by c1 to get the RSL psihat function evaluated at hc
%     !---------------------------------------------------------------------
% 
%     ! In the RSL theory, c1 and c2 represent the scaled magnitude and
%     ! height over which the RSL theory modifies MOST via the psihat functions:
%     !
%     !
%     !     z
%     !     ^
%     !     |                                    . ___
%     !     |                                    .  ^
%     !     |                                   .   |
%     !     |                                  .    |
%     !     |                                 .     |
%     !     |                               .       |
%     !     |                             .        c2
%     !     |                          .            |
%     !     |                       .               |
%     !     |                  .                    |
%     !     |           .                           |
%     !     |   .                                 _\/_
%     !     -------------------------------------------> u
%     !         |<-------------- c1 --------------->|
% 
%     ! Evaluate the roughness sublayer psihat function for momentum at
%     ! the height za and at the canopy height hc. Values for psihat are obtained
%     ! from a look-up table. Here, heights are above the canopy for compatibility
%     ! with the supplied look-up table. These heights are also scaled by dt = hc-d
%     ! so that the look-up table uses (za-hc)/dt and (hc-hc)/dt. Also the term
%     ! dt in the integration of psihat is scaled by the Obukhov length L (dt/obu).
%     ! This means that the returned psihat value needs to be scaled (multiplied) by
%     ! c1 before it fully represents psihat as it appears in the RSL equations.

%     select case (turb_type)
%        case (1, 3)

vkc = 0.4;
c2 = 0.5;


%[dtLgridM,zdtgridM,psigridM,dtLgridH,zdtgridH,psigridH] = LookupPsihatINI;
%replaced by outside lookup table build to save time
dtLgridM = vegetation.PsiLookup.dtLgridM;
zdtgridM = vegetation.PsiLookup.zdtgridM;
psigridM = vegetation.PsiLookup.psigridM;
dtLgridH = vegetation.PsiLookup.dtLgridH;
zdtgridH = vegetation.PsiLookup.zdtgridH;
psigridH = vegetation.PsiLookup.psigridH;


[psihat1] = LookupPsihat ((za-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM);
[psihat2] = LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM);

%           call LookupPsihat ((za-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM, psihat1)
%           call LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM, psihat2)
%        case (2)
%           psihat1 = 0._r8
%           psihat2 = 0._r8
%     end select

phim = phi_m_monin_obukhov(dt/obu);
c1 = (1. - vkc / (2. * beta * phim)) * exp(0.5*c2);
psim = -psi_m_monin_obukhov(zeta) + psi_m_monin_obukhov(dt/obu) + c1*(psihat1 - psihat2) + vkc / beta;

%     ! Now do the same for heat

%     select case (turb_type)
%        case (1, 3)
[psihat1] = LookupPsihat ((za-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH);
[psihat2] = LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH);
%           call LookupPsihat ((za-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH, psihat1)
%           call LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH, psihat2)
%        case (2)
%           psihat1 = 0._r8
%           psihat2 = 0._r8
%     end select

phic = phi_c_monin_obukhov(dt/obu);
c1 = (1. - PrSc*vkc / (2. * beta * phic)) * exp(0.5*c2);
psic = -psi_c_monin_obukhov(zeta) + psi_c_monin_obukhov(dt/obu) + c1*(psihat1 - psihat2);

end

