
%     !
%     ! !DESCRIPTION:
%     ! Initialize the look-up tables needed to calculate the RSL psihat functions
%     !
%     ! !USES:
%     use clm_varcon, only : nZ, nL, dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
%     use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit
%     !
%     ! !ARGUMENTS:
%     implicit none
%     !
%     !LOCAL VARIABLES
%     character(len=256) :: fin          ! File name
%     integer :: nin                     ! Fortran unit number
%     integer :: ii, jj                  ! Looping indices
%     !----------------------------------------

function [dtLgridM,zdtgridM,psigridM,dtLgridH,zdtgridH,psigridH] = LookupPsihatINI()
table_M = dlmread('psihatM.dat');
dtLgridM = table_M(1,2:end);
zdtgridM = table_M(2:end,1);
psigridM = table_M(2:end,2:end);

table_H = dlmread('psihatH.dat');
dtLgridH = table_H(1,2:end);
zdtgridH = table_H(2:end,1);
psigridH = table_H(2:end,2:end);

 
% % Examples
% % Read and display the file fgetl.m one line at a time:
% fid = fopen('psihatM.dat');
% tline = fgetl(fid);
% while ischar(tline)
%     disp(tline)
%     tline = fgetl(fid);
% end
% fclose(fid);
% 
% 
%     fin = '../canopy/psihatM.dat'
%     nin = shr_file_getUnit()
%     open (unit=nin, file=trim(fin), action="read", status="old", err=10)
%     read (nin,*,err=11) dtLgridM(1,1),(dtLgridM(1,jj), jj = 1,nL)
%     do ii = 1  nZ
%       read (nin,*,err=12) zdtgridM(ii,1),(psigridM(ii,jj), jj = 1,nL)
%     end do
%     close (nin)
%     call shr_file_freeUnit (nin)
%     
%     fin = '../canopy/psihatH.dat'
%     nin = shr_file_getUnit()
%     open (unit=nin, file=trim(fin), action="read", status="old", err=20)
%     read (nin,*,err=21) dtLgridH(1,1),(dtLgridH(1,jj), jj = 1,nL)
%     do ii = 1,nZ
%       read (nin,*,err=22) zdtgridH(ii,1),(psigridH(ii,jj), jj = 1,nL)
%     end do
%     close (nin)
%     call shr_file_freeUnit (nin)
% 
%     return

% 10  continue
%     call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error opening psihatM.dat')
% 
% 11  continue
%     call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error reading dtLgridM')
% 
% 12  continue
%     call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error reading psigridM')
% 
% 20  continue
%     call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error opening psihatH.dat')
% 
% 21  continue
%     call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error reading dtLgridH')
% 
% 22  continue
%     call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error reading psigridH')

    end
