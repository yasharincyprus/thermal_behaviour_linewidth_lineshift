
pos_ph = [88,229,357,472,110.6,126.2,180.1,226.75,292,342.5,387.1];
pos_ph = sort(pos_ph');
len_pos_ph = numel(pos_ph);
fileID = fopen('errorlog.txt','w');

if mod(len_pos_ph,2) ~= 0
    comb_num = nchoosek(len_pos_ph,(len_pos_ph + 1)/2);
else
    comb_num = nchoosek(len_pos_ph,len_pos_ph/2);
end

comb_cell = cell(comb_num,len_pos_ph-1);
combos_sz = zeros(len_pos_ph-1,1);

for i = 2:len_pos_ph
    combos = combntns(pos_ph,i);
    sz = size(combos);
    combos_sz(i-1,1) = sz(1);
    for j = 1:sz(1)
        comb_cell{j,i-1} = combos(j,:);
    end
end

szcc = size(comb_cell);

for i = 1:szcc(2)
   try 
    switch i
        case 1
            for j = 1:combos_sz(1,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD2
                temp_fits_df2
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 2
            for j = 1:combos_sz(2,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD3
                temp_fits_df3
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 3
            for j = 1:combos_sz(3,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD4
                temp_fits_df4
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 4
            for j = 1:combos_sz(4,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end 
                temp_fits_deb_fixed_RD5
                temp_fits_df5
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 5
            for j = 1:combos_sz(5,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD6
                temp_fits_df6
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 6
            for j = 1:combos_sz(6,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD7
                temp_fits_df7
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 7
            for j = 1:combos_sz(7,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD8
                temp_fits_df8
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 8
            for j = 1:combos_sz(8,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD9
                temp_fits_df9
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 9
            for j = 1:combos_sz(9,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end 
                temp_fits_deb_fixed_RD10
                temp_fits_df10
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 10
            for j = 1:combos_sz(10,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD11
                temp_fits_df11
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 11
            for j = 1:combos_sz(11,1)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD12
                temp_fits_df12
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        case 12
            for j = 1:combos_sz(12)
                for k = 1:i+1
                    eval(['E' i '= comb_cell{' j ',' i '}(1,' k ')']);
                end
                temp_fits_deb_fixed_RD13
                temp_fits_df13
                comb_cell{j,i}(2,1) = alpha_final;
                comb_cell{j,i}(2,2) = debye;
                comb_cell{j,i}(2,3) = ssq_out_min;
                for l = 1:length(npar_n)
                    comb_cell{j,i}(2,l+3) = npar_n(l,1);
                end
            end
        otherwise
            disp('All possible direct phonon combinations have been created and tested successfully.')
    end
   catch
       msg = strcat('One of the programmes "temp_fits_deb_fixed_RD',num2str(i+1),'.m and "temp_fits_df',num2str(i+1),'.m throws an error.');
       fprintf(fileID,'\r\n%s',msg);
       continue
   end
end

T = cell2table(comb_cell(2:end,:));
writetable(T,'combs_params_linewidths.csv')
