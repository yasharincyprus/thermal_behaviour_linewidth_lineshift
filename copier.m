% ofstr = 'temp_fits_deb_fixed_RD3.m';
% offstr = 'temp_fits_df3.m';
% offfstr = 'ssqmin_cfit_df3.m';
% offffstr = 'cfit_output_df6.m';
% for i = 7:13
% fstr = strcat('temp_fits_deb_fixed_RD',num2str(i),'.m');
% copyfile(ofstr,fstr)
% ffstr = strcat('temp_fits_df',num2str(i),'.m');
% copyfile(offstr,ffstr)
% fffstr = strcat('ssqmin_cfit_df',num2str(i),'.m');
% copyfile(offfstr,fffstr)
% ffffstr = strcat('cfit_output_df',num2str(i),'.m');
% copyfile(offffstr,ffffstr)
% end

    
for i = 1:5
    try
        temp_shift
        pause
    catch
        clc
        disp('This programme is corrupted.')
        pause
        continue
    end
end