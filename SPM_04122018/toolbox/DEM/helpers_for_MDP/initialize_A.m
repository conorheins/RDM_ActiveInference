% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(d);
for f = 1:Nf
    Ns(f) = numel(d{f});
end
No    = [7 9];
Ng    = numel(No);
for g = 1:Ng
    A{g} = zeros([No(g),Ns]);
end
for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)
        for f3 = 1:Ns(3)
            for f4 = 1:Ns(4)
                for f5 = 1:Ns(5)
                    for f6 = 1:Ns(6)
                    
                        if f1 == 1, a = {'90','0';'null','null'}; end % Scene A
                        if f1 == 2, a = {'0','270';'null','null'}; end % Scene B
                        if f1 == 3, a = {'270','180';'null','null'}; end % Scene C
                        if f1 == 4, a = {'180','90';'null','null'}; end % Scene D
                    
                        % flip cues according to hidden (invariants) states
                        %----------------------------------------------------------
                        % first compute the diagonal switch, then the
                        % horizontal/vertical flips
                        
                        % this corresponds to the transposition (moving
                        % upper right quadrant to lower left)
                        if f6 == 2
                            a{2} = a{3}; a{3} = 'null';
                        end
                        
                        % this corresponds to the 'diagonalization' (moving
                        % the upper right to the lower right)
                        if f5 == 2
                            a{4} = a{3}; a{3} = 'null';
                        end

                        if f3 == 2, a = flipud(a); end
                        if f4 == 2, a = fliplr(a); end
                    
                        % what: A{1} {'null','0','90','180','270','right,'wrong'}
                        %==========================================================
                        if f2 == 1
                            
                            % at fixation location
                            %----------------------------------------------------------
                            A{1}(1,f1,f2,f3,f4,f5,f6) = true;
                            
                        elseif f2 > 1 && f2 < 6
                            
                            
                            % saccade to cue location
                            %----------------------------------------------------------
                            A{1}(1,f1,f2,f3,f4,f5,f6) = strcmp(a{f2 - 1},'null');
                            A{1}(2,f1,f2,f3,f4,f5,f6) = strcmp(a{f2 - 1},'0');
                            A{1}(3,f1,f2,f3,f4,f5,f6) = strcmp(a{f2 - 1},'90');
                            A{1}(4,f1,f2,f3,f4,f5,f6) = strcmp(a{f2 - 1},'180');
                            A{1}(5,f1,f2,f3,f4,f5,f6) = strcmp(a{f2 - 1},'270');
                            
                        elseif f2 > 5
                            
                            % saccade choice location
                            %------------------------------------------------------
                            A{1}(6,f1,f2,f3,f4,f5,f6) = (f2 - 5) == f1;
                            A{1}(7,f1,f2,f3,f4,f5,f6) = (f2 - 5) ~= f1;
                            
                        end
                        
                        % where: A{2} {'start','1',...,'4','A','B','C',D'}
                        %----------------------------------------------------------
                        A{2}(f2,f1,f2,f3,f4,f5,f6) = 1;
                        
                    end
                end
            end
        end
    end
end

for g = 1:Ng
    A{g} = double(A{g});
end