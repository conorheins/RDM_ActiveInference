scenes = {'Scene A','Scene B','Scene C','Scene D'};
flipud_flags = [0 1];
fliplr_flags = [0 1];
diag_flags = [0 1];
trans_flags = [0 1];

for ii = 1:length(scenes)
    for jj = 1:length(flipud_flags)
        for kk = 1:length(fliplr_flags)
            for ll = 1:length(diag_flags)
                for pp = 1:length(trans_flags)
                    displayAmatrix(A,scenes{ii},flipud_flags(jj),fliplr_flags(kk),diag_flags(ll),trans_flags(pp));
                    if flipud_flags(jj)
                        vert_str = 'vertical flip';
                    else
                        vert_str = 'no vertical flip';
                    end
                    if fliplr_flags(kk)
                        lr_str = 'horizontal flip';
                    else
                        lr_str = 'no horizontal flip';
                    end
                    if diag_flags(ll)
                        diag_str = 'diagonalization';
                    else
                        diag_str = 'no diagonalization';
                    end
                    if trans_flags(pp)
                        trans_str = 'transpose';
                    else
                        trans_str = 'no transpose';
                    end
                
                    if ~flipud_flags(jj) && ~fliplr_flags(kk) && ~diag_flags(ll) && ~trans_flags(pp)
                        title(sprintf('%s, no spatial transformations',scenes{ii}))
                    else
                        title(sprintf('%s, %s, %s, %s, %s',scenes{ii},vert_str,lr_str,diag_str,trans_str));
                    end
                    filename = fullfile('SceneConfigs12062018',sprintf('%s_%d%d%d%d.png',scenes{ii},flipud_flags(jj),fliplr_flags(kk),diag_flags(ll),trans_flags(pp)));
                    saveas(gcf,filename);
                end
            end
        end
    end
end