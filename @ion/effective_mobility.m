function eff_mobility=effective_mobility(obj, pH, I)
    % This function takes the ion and a pH and an ionic strength. It calls the ionization
    % fraction function and uses this information to compute the
    % effective mobility of the ion. 
			
	%%%%%%%%%
	% This function uses ionic strength for ionization fraction 
	% calculation, but it does not use it for fi_mobility.
	%%%%%%%%%
			
	if ~exist('I', 'var')
		I=0;
    end
            
	if ~isempty(obj.fi_mobility_effective)
		fi_mobility=obj.fi_mobility_effective;
	else
		fi_mobility=obj.robinson_stokes_mobility(I);
	end
	
	eff_mobility=zeros(size(pH));
    for i = 1:length(pH)
        i_frac=ionization_fraction(obj, pH(i), I);
        eff_mobility(i)=sum(i_frac.*fi_mobility);
    end
			
end