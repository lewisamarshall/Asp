function m_cond=molar_conductivity(obj, pH, I)
    % This function takes the ion and the pH and computes the molar
    % conductivity. This function can take an ionic strength. 
	% Provides conducitivity in Siemens per meter per mole.
            
	if ~exist('I', 'var')
		I=0;
	end
	
	if ~isempty(obj.fi_mobility_effective)
		fi_mobility=obj.fi_mobility_effective;
	else
		fi_mobility=obj.robinson_stokes_mobility(I);
	end
			
    i_frac=ionization_fraction(obj, pH, I);
            
    m_cond=sum(obj.F*obj.z.*i_frac.*fi_mobility*obj.Lpm3);
            
end