function m_cond=molar_conductivity(obj, pH, I)
    % This function takes the ion and the pH and computes the molar
    % conductivity. This function can take an ionic strength. 
	% Provides conducitivity in Siemens per meter per mole.
			
	%%%%%%%%%%%%%%%%
	% Does not yet use the ionic strength corrected mobilities.
	%%%%%%%%%%%%%%%%
            
	if ~exist('I', 'var')
		I=0;
	end
			
    i_frac=ionization_fraction(obj, pH, I);
            
    m_cond=sum(obj.F*obj.z.*i_frac.*obj.fi_mobility*obj.Lpm3);
            
end