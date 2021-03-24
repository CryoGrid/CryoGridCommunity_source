function [beta] = GetBeta (beta_neutral, LcL)
    %
    % DESCRIPTION:
    % Calculate beta = u* / u(h) for current Obukhov length
    %
    % ARGUMENTS:
    % Input:
    
   
    
%     real(r8), intent(in)  :: beta_neutral    % Neutral value for beta = u*/u(h), ~0.3-0.5
%     real(r8), intent(in)  :: LcL             % Canopy density scale (Lc) / Obukhov length (obu)
%     real(r8), intent(out) :: beta            % Value of u*/u(h) at canopy top
    %
    % LOCAL VARIABLES:
    %aa, bb, cc, dd, qq, rr       % Terms for quadratic or cubic solution
    %LcL_val                      % LcL with limits applied
    %---------------------------------------------------------------------

    LcL_val = LcL;
    
    if (LcL_val <= 0 )  

       % Unstable case: quadratic equation for beta^2 at LcL_val

       bb = 16 .* LcL_val.* beta_neutral^4     ;
       beta = sqrt( 0.5.*(-bb + sqrt(bb^2 + 4 .* beta_neutral^4)) );

    else

       % Stable case: cubic equation for beta at LcL_val

       aa = 5 .* LcL_val;
       bb = 0 ;
       cc = 1 ;
       dd = -beta_neutral;
       qq = (2.*bb^3 - 9.*aa*bb*cc + 27.*(aa^2)*dd)^2 - 4.*(bb^2 - 3.*aa*cc)^3;
       qq = sqrt(qq);
       rr = 0.5 .* (qq + 2.*bb^3 - 9.*aa*bb*cc + 27.*(aa^2)*dd);
       rr = rr^(1 /3 );
       beta = -(bb+rr)/(3.*aa) - (bb^2 - 3.*aa*cc)/(3.*aa*rr)    ;

    end

    beta = min(0.50 , max(beta,0.20));

end