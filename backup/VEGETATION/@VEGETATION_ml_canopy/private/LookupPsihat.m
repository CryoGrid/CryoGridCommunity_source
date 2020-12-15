%   !-----------------------------------------------------------------------
function  [psihat] = LookupPsihat (zdt, dtL, zdtgrid, dtLgrid, psigrid)
%     !
%     ! !DESCRIPTION:
%     ! Determines the RSL function psihat as provided through a look-up table
%     ! for input values of zdt and dtL. Linearly interpolates between values
%     ! supplied on the look-up table grid defined by zdtgrid, dtLgrid, psigrid.
%     !
%     ! NOTE: The psihat presented in Harman and Finnigan (2007,2008) and Harman (2012)
%     ! has been re-written in non-dimensional form such that it now appears as:
%     !
%     ! psihat(z) = c1 * A(z/(beta^2*Lc),(beta^2*Lc)/L)
%     !
%     ! This routine gets the value of A from a look-up table. Noting that dt=beta^2*Lc,
%     ! this routine therefore requires values of z/dt and dt/L. In addition, this means
%     ! that the returned psihat value needs to be scaled (multiplied) by c1 before it fully
%     ! represents psihat as it appears in the RSL equations.
%     !
%     ! !USES:
%     use clm_varcon, only : nZ, nL
%     !
%     ! !ARGUMENTS:
%     implicit none
%     real(r8), intent(in) :: zdt            ! Height (above canopy) normalized by dt
%     real(r8), intent(in) :: dtL            ! dt/L (displacement height/Obukhov length)
%     real(r8), intent(in) :: zdtgrid(nZ,1)  ! Grid of zdt on which psihat is given
%     real(r8), intent(in) :: dtLgrid(1,nL)  ! Grid of dtL on which psihat is given
%     real(r8), intent(in) :: psigrid(nZ,nL) ! Grid of psihat values
%     real(r8), intent(out):: psihat         ! Value of psihat
%     !
%     !LOCAL VARIABLES
%     integer  :: ii, jj                     ! Looping indices
%     integer  :: L1, L2, Z1, Z2             ! Grid indices for psihat sought
%     real(r8) :: wL1, wL2, wZ1, wZ2         ! Weights for averaging
%     !---------------------------------------------------------------------
% 
%     ! Find indices and weights for dtL values which bracket the specified dtL

    L1 = 0; 
    L2 = 0;
    if (dtL <= dtLgrid(1))
       L1 = 1;
       L2 = 1;
       wL1 = 0.5;
       wL2 = 0.5;
    elseif (dtL > dtLgrid(end))
       L1 = dtLgrid(end);
       L2 = dtLgrid(end);
       wL1 = 0.5;
       wL2 = 0.5;
    else
       for jj = 1:length(dtLgrid)-1
          if ((dtL <= dtLgrid(jj+1)) && (dtL > dtLgrid(jj)))
             L1 = jj;
             L2 = jj + 1;
             wL1 = (dtLgrid(L2) - dtL) / (dtLgrid(L2) - dtLgrid(L1));
             wL2 = 1. - wL1;
          end
       end
    end

    if (L1 == 0 || L2 == 0)
       disp(' ERROR: CanopyTurbulenceMod: LookupPsihat error, indices L1 and L2 not found');
    end

%     ! Find indices and weights for zdt values which bracket the specified zdt

    Z1 = 0; 
    Z2 = 0;
    if ( zdt > zdtgrid(1,1))
       Z1 = 1;
       Z2 = 1;
       wZ1 = 0.5;
       wZ2 = 0.5;
    elseif (zdt < zdtgrid(end,1))
       Z1 = length(zdtgrid);
       Z2 = length(zdtgrid);
       wZ1 = 0.5;
       wZ2 = 0.5;
    else
       for ii = 1:length(zdtgrid)-1
          if ( (zdt >= zdtgrid(ii+1,1)) && (zdt < zdtgrid(ii,1)) )
             Z1 = ii;
             Z2 = ii + 1;
             wZ1 = (zdt - zdtgrid(ii+1,1)) / (zdtgrid(ii,1) - zdtgrid(ii+1,1));
             wZ2 = 1. - wZ1;
          end
       end
    end

    if (Z1 == 0 || Z2 == 0)
       error(' ERROR: CanopyTurbulenceMod: LookupPsihat error, indices Z1 and Z2 not found');
    end

%     ! Calculate psihat as a weighted average of the values of psihat on the grid
    psihat = wZ1*wL1*psigrid(Z1,L1) + wZ2*wL1*psigrid(Z2,L1) + wZ1*wL2*psigrid(Z1,L2) + wZ2*wL2*psigrid(Z2,L2);
    
    
    end
           
    
    
    
    
    